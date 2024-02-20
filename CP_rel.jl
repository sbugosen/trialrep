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

function build_cp_opf(data::Dict{String,Any}, model=Model())
    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    @variable(model, ref[:bus][i]["vmin"]^2 <= w[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"]^2, start=1.001)

    @variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    @variable(model, ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"])
    
    #@objective(model, Min,
    #    sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen])
    #    #+sum(dcline["cost"][1]*p_dc[from_idx[i]]^2 + dcline["cost"][2]*p_dc[from_idx[i]] + dcline["cost"][3] for (i,dcline) in ref[:dcline])
    #)

    @objective(model, Max, sum(pg[i] for (i,gen) in ref[:gen]) + sum(qg[i] for (i,gen) in ref[:gen]))    

    @variable(model, -Inf <= p[(l,i,j) in ref[:arcs]] <= Inf)
    @variable(model, -Inf <= q[(l,i,j) in ref[:arcs]] <= Inf)    
    
    bus_loads = [try ref[:load][l] catch error; nothing end for l = 1:length(ref[:bus_loads])]
    generators = reduce(vcat, [value for (key,value) in ref[:bus_gens] if !isempty(value)])
    bus_shunts = [try ref[:shunt][s] catch error; nothing end for s = 1:length(ref[:bus_shunts])]
    bus_shunts_filtered = filter(x -> x !== nothing, bus_shunts)
    shunts_gs = [shunt["gs"] for shunt in bus_shunts_filtered]
    shunts_bs = [shunt["bs"] for shunt in bus_shunts_filtered]
    buses = reduce(vcat, [key for (key,value) in ref[:bus_shunts] if !isempty(value)], init = [])
    
    A1 = Dict()
    B1 = Dict()
    A2 = Dict()
    B2 = Dict()

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

        A1[i,j] = -resis_fr*(
                   (vfub*vtub * cos(d) * (vflb*vtlb - vfub*vtub) + vtub*cos(d)*st*w_fr + vfub*cos(d)*sf*w_to) / (sf * st)
                   - w_fr
        )

        B1[i,j] = -resis_to*(
                   (vtub*vfub * cos(d) * (vtlb*vflb - vtub*vfub) + vfub*cos(d)*sf*w_to + vtub*cos(d)*st*w_fr) / (st * sf)
                   - w_to
        )

        A2[i,j] = -resis_fr*(
                   (vflb*vtlb * cos(d) * (vfub*vtub - vflb*vtlb) + vtlb*cos(d)*st*w_fr + vflb*cos(d)*sf*w_to) / (sf * st)
                   - w_fr
        )

        B2[i,j] = -resis_to*(
                   (vtlb*vflb * cos(d) * (vtub*vfub - vtlb*vflb) + vflb*cos(d)*sf*w_to + vtlb*cos(d)*st*w_fr) / (st * sf)
                   - w_to
        )

    end
   
    @constraint(model,
                sum(pg[g] for g in generators) -
                sum(load["pd"] for load in bus_loads if load !== nothing) - 
                sum(shunts_gs[i]*w[buses[i]] for i in 1:length(shunts_gs)) 
                <= sum(A1[branch["f_bus"],branch["t_bus"]] for (k,branch) in ref[:branch]) 
                + sum(B1[branch["f_bus"],branch["t_bus"]] for (k,branch) in ref[:branch])
    )

    @constraint(model,
                sum(pg[g] for g in generators) -
                sum(load["pd"] for load in bus_loads if load !== nothing) -
                sum(shunts_gs[i]*w[buses[i]] for i in 1:length(shunts_gs))
                <= sum(A2[branch["f_bus"],branch["t_bus"]] for (k,branch) in ref[:branch]) 
                + sum(B2[branch["f_bus"],branch["t_bus"]] for (k,branch) in ref[:branch])
    )

    @constraint(model,
                sum(pg[g] for g in generators) -
                sum(load["pd"] for load in bus_loads if load !== nothing) - 
                sum(shunts_gs[i]*w[buses[i]] for i in 1:length(shunts_gs)) >= 0 
               )

    @constraint(model,
                sum(qg[g] for g in generators) -
                sum(load["qd"] for load in bus_loads if load !== nothing) + 
                sum(shunts_bs[i]*w[buses[i]] for i in 1:length(shunts_bs)) >= 0 
               )
    
   return model
end

import AmplNLWriter
using AmplNLWriter
#data = PowerModels.parse_file("/Users/sbugosen/Repositories/trialrep/pglib-opf/sad/pglib_opf_case10192_epigrids__sad.m")
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
    build_cp_opf(data, model)
    optimize!(model)
    time = JuMP.solve_time(model)
    push!(times, time)
    push!(obj_vals, (JuMP.objective_value(model)))
    push!(term_st, (JuMP.termination_status(model)))
end

df = DataFrame(Datasets = readdir("pglib-opf/sad"), Objective = obj_vals, Time = times, Termination = term_st)
CSV.write("cp_max_sad_lnc_FINAL_shunts.csv",df)
 

