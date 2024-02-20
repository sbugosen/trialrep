import PowerModels
#import Ipopt
import JuMP
import DataFrames
import CSV
using JuMP
using PowerModels
#using Ipopt
using DataFrames
using CSV

function build_nf_opf(data::Dict{String,Any}, model=Model())
    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    
    @variable(model, ref[:bus][i]["vmin"]^2 <= w[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"]^2, start=1.001)
    @variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    @variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])

    @variable(model, -Inf <= p[(l,i,j) in ref[:arcs]] <= Inf)
    @variable(model, -Inf <= q[(l,i,j) in ref[:arcs]] <= Inf)

    from_idx = Dict(arc[1] => arc for arc in ref[:arcs_from_dc])

    #@objective(model, Min,
    #           sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen])
    #           #+sum(dcline["cost"][1]*p_dc[from_idx[i]]^2 + dcline["cost"][2]*p_dc[from_idx[i]] + dcline["cost"][3] for (i,dcline) in ref[:dcline])
    #          )

    @objective(model, Max, sum(pg[i] for (i,gen) in ref[:gen])+ sum(qg[i] for (i,gen) in ref[:gen]))

    for (i,bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Bus KCL
        @constraint(model,
            sum(p[a] for a in ref[:bus_arcs][i]) 
            #+sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) 
            == sum(pg[g] for g in ref[:bus_gens][i]) -
            sum(load["pd"] for load in bus_loads) - 
            sum(shunt["gs"] for shunt in bus_shunts)*w[i]
        )

        @constraint(model,
            sum(q[a] for a in ref[:bus_arcs][i])
            #+sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i])
            == sum(qg[g] for g in ref[:bus_gens][i]) -
            sum(load["qd"] for load in bus_loads)+   
            sum(shunt["bs"] for shunt in bus_shunts)*w[i]
        )
    end

    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])
        bp_idx = (branch["f_bus"], branch["t_bus"])

        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]

        w_fr = w[branch["f_bus"]]
        w_to = w[branch["t_bus"]]

        # Line Flow
        r, x = ref[:branch][i]["br_r"], ref[:branch][i]["br_x"] 
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]^2
        
        # Real Line Losses Constraint (lower bound)
        @constraint(model, p_fr + p_to >= 0)
        @constraint(model, q_fr + q_to >= 0)

        resis_fr = (g+g_fr)*tm / ((g+g_fr)^2+(b+b_fr)^2)
        reac_fr = -(b+b_fr)*tm / ((g+g_fr)^2+(b+b_fr)^2)

        @constraint(model, (-resis_fr*q_fr + reac_fr*p_fr) <= tan(branch["angmax"])*(w_fr - resis_fr*p_fr + reac_fr*q_fr))
        @constraint(model, (-resis_fr*q_fr + reac_fr*p_fr) >= tan(branch["angmin"])*(w_fr - resis_fr*p_fr + reac_fr*q_fr))
        
        resis_to = (g+g_to)*tm / ((g+g_to)^2+(b+b_to)^2)
        reac_to = -(b+b_to)*tm / ((g+g_to)^2+(b+b_to)^2)

        @constraint(model, (-resis_to*q_to + reac_to*p_to) <= tan(branch["angmax"])*(w_to - resis_to*p_to + reac_to*q_to))
        @constraint(model, (-resis_to*q_to + reac_to*p_to) >= tan(branch["angmin"])*(w_to - resis_to*p_to + reac_to*q_to))

    end
    
    for (k,branch) in ref[:branch]
        i,j = branch["f_bus"], branch["t_bus"]
        f_idx = (k, branch["f_bus"], branch["t_bus"])
        t_idx = (k, branch["t_bus"], branch["f_bus"])
        bp_idx = (branch["f_bus"], branch["t_bus"])

        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]

        w_fr = w[branch["f_bus"]]
        w_to = w[branch["t_bus"]]

        # Line Flow
        #r, x = ref[:branch][i]["br_r"], ref[:branch][i]["br_x"] 
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]^2
        
        buspair = ref[:buspairs][i,j]
        vfub = buspair["vm_fr_max"]
        vflb = buspair["vm_fr_min"]
        vtub = buspair["vm_to_max"]
        vtlb = buspair["vm_to_min"]
        tdub = buspair["angmax"]
        tdlb = buspair["angmin"]
        
        resis_to = (g+g_to)*tm / ((g+g_to)^2+(b+b_to)^2)
        reac_to = -(b+b_to)*tm / ((g+g_to)^2+(b+b_to)^2)
        
        resis_fr = (g+g_fr)*tm / ((g+g_fr)^2+(b+b_fr)^2)
        reac_fr = -(b+b_fr)*tm / ((g+g_fr)^2+(b+b_fr)^2)

        phi = (tdub + tdlb)/2
        d   = (tdub - tdlb)/2

        sf = vflb + vfub
        st = vtlb + vtub
        
        A1 = -resis_fr*(
                (vfub*vtub * cos(d) * (vflb*vtlb - vfub*vtub) + vtub*cos(d)*st*w_fr + vfub*cos(d)*sf*w_to) / (sf * st) - w_fr
        )
 
        B1 = -resis_to*(
                (vtub*vfub * cos(d) * (vtlb*vflb - vtub*vfub) + vfub*cos(d)*sf*w_to + vtub*cos(d)*st*w_fr) / (st * sf) - w_to
        )
 
        A2 = -resis_fr*(
                (vflb*vtlb * cos(d) * (vfub*vtub - vflb*vtlb) + vtlb*cos(d)*st*w_fr + vflb*cos(d)*sf*w_to) / (sf * st) - w_fr
        )
 
        B2 = -resis_to*(
                (vtlb*vflb * cos(d) * (vtub*vfub - vtlb*vflb) + vflb*cos(d)*sf*w_to + vtlb*cos(d)*st*w_fr) / (st * sf) - w_to
        )

        @constraint(model, p_fr + p_to <= A1 + B1)
        @constraint(model, p_fr + p_to <= A2 + B2)

        # Voltage Product Relaxation Upperbound
        #@constraint(model, sf*st*(cos(phi)*((w[i] - resis_to*p_to + reac_to*q_to)) + sin(phi)*(-resis_to*q_to + reac_to*p_to)) - vtub*cos(d)*st*w[i] - vfub*cos(d)*sf*w[j] >=  vfub*vtub*cos(d)*(vflb*vtlb - vfub*vtub))
        #@constraint(model, sf*st*(cos(phi)*((w[i] - resis_to*p_to + reac_to*q_to)) + sin(phi)*(-resis_to*q_to + reac_to*p_to)) - vtlb*cos(d)*st*w[i] - vflb*cos(d)*sf*w[j] >= -vflb*vtlb*cos(d)*(vflb*vtlb - vfub*vtub))
        ##
        #@constraint(model, sf*st*(cos(phi)*((w[j] - resis_fr*p_fr + reac_fr*q_fr)) + sin(phi)*(-resis_fr*q_fr + reac_fr*p_fr)) - vfub*cos(d)*sf*w[j] - vtub*cos(d)*st*w[i] >=  vtub*vfub*cos(d)*(vtlb*vflb - vtub*vfub))
        #@constraint(model, sf*st*(cos(phi)*((w[j] - resis_fr*p_fr + reac_fr*q_fr)) + sin(phi)*(-resis_fr*q_fr + reac_fr*p_fr)) - vflb*cos(d)*sf*w[j] - vtlb*cos(d)*st*w[i] >= -vtlb*vflb*cos(d)*(vtlb*vflb - vtub*vfub))

    end

    return model
end

import AmplNLWriter
using AmplNLWriter
#data = PowerModels.parse_file("pglib_opf_case73_ieee_rts.m")
#model = Model(() -> AmplNLWriter.Optimizer("ipopt"))
#build_cp_opf(data, model)
#optimize!(model)


obj_vals = Vector{Float64}()
term_st  = Vector()
times = Vector{Float64}()
files2 = joinpath.(abspath("pglib-opf/sad"), readdir("pglib-opf/sad"))

for file in files2
    data = PowerModels.parse_file(file)
    model = Model(() -> AmplNLWriter.Optimizer("ipopt"))
    build_nf_opf(data, model)
    optimize!(model)
    time = JuMP.solve_time(model)
    push!(times, time)
    push!(obj_vals, (JuMP.objective_value(model)))
    push!(term_st, (JuMP.termination_status(model)))
end

df = DataFrame(Datasets = readdir("pglib-opf/sad"), Objective = obj_vals, Time = times, Termination = term_st)
CSV.write("nf_max_sad_lnc_FINAL.csv",df)

#println(has_duals(model))
#
#println("Active constraints:")
#for con in all_constraints(model, include_variable_in_set_constraints = true)
#    println(con," ","dual:", dual(con))
#end
