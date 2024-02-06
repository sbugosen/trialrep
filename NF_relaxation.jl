#import Pkg; Pkg.add("PowerModels")
import PowerModels
import Ipopt
import JuMP
using JuMP
using PowerModels
using Ipopt

#data = PowerModels.parse_file("../goc-c3-datasets/pglib/pglib-opf/pglib_opf_case118_ieee.m")
data = PowerModels.parse_file("pglib_opf_case118_ieee.m")
#data = PowerModels.parse_file("pglib_opf_case73_ieee_rts.m")
#data = PowerModels.parse_file("pglib_opf_case162_ieee_dtc.m")
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
    for arc in ref[:arcs]
        (l,i,j) = arc
        branch = ref[:branch][l]
        if haskey(branch, "rate_a")
            JuMP.set_lower_bound(p[arc], -branch["rate_a"])
            JuMP.set_upper_bound(p[arc],  branch["rate_a"])
            JuMP.set_lower_bound(q[arc], -branch["rate_a"])
            JuMP.set_upper_bound(q[arc],  branch["rate_a"])
        end
    end
    @variable(model, p_dc[a in ref[:arcs_dc]])
    @variable(model, q_dc[a in ref[:arcs_dc]])

    for (l,dcline) in ref[:dcline]
        f_idx = (l, dcline["f_bus"], dcline["t_bus"])
        t_idx = (l, dcline["t_bus"], dcline["f_bus"])

        JuMP.set_lower_bound(p_dc[f_idx], dcline["pminf"])
        JuMP.set_upper_bound(p_dc[f_idx], dcline["pmaxf"])
        JuMP.set_lower_bound(q_dc[f_idx], dcline["qminf"])
        JuMP.set_upper_bound(q_dc[f_idx], dcline["qmaxf"])

        JuMP.set_lower_bound(p_dc[t_idx], dcline["pmint"])
        JuMP.set_upper_bound(p_dc[t_idx], dcline["pmaxt"])
        JuMP.set_lower_bound(q_dc[f_idx], dcline["qmint"])
        JuMP.set_upper_bound(q_dc[f_idx], dcline["qmaxt"])
    end

    from_idx = Dict(arc[1] => arc for arc in ref[:arcs_from_dc])

    @objective(model, Min,
        sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen]) 
       + sum(dcline["cost"][1]*p_dc[from_idx[i]]^2 + dcline["cost"][2]*p_dc[from_idx[i]] + dcline["cost"][3] for (i,dcline) in ref[:dcline])
    )
    
    for (i,bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Bus KCL
        @constraint(model,
            sum(p[a] for a in ref[:bus_arcs][i]) +
            sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) -
            sum(load["pd"] for load in bus_loads) -
            sum(shunt["gs"] for shunt in bus_shunts)*w[i]
        )

        @constraint(model,
            sum(q[a] for a in ref[:bus_arcs][i]) +
            sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==
            sum(qg[g] for g in ref[:bus_gens][i]) -
            sum(load["qd"] for load in bus_loads) +
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
        
        # Imaginary Line Losses COnstraint (lower bound)
        #@constraint(model, q_fr + q_to >= (w_to / tm) + w_fr)

        # Phase Angle Difference Limit
        #@constraint(model, (-r*q_fr + x*p_fr) <= tan(branch["angmax"])*(w_fr - r*p_fr + x*q_fr))
        #@constraint(model, (-r*q_fr + x*p_fr) >= tan(branch["angmin"])*(w_fr - r*p_fr + x*q_fr))
        
        #resis_fr = (g+g_fr)*tm / ((g+g_fr)^2+(b+b_fr)^2)
        #reac_fr = -(b+b_fr)*tm / ((g+g_fr)^2+(b+b_fr)^2)

        #@constraint(model, (-resis_fr*q_fr + reac_fr*p_fr) <= tan(branch["angmax"])*(w_fr - resis_fr*p_fr + reac_fr*q_fr))
        #@constraint(model, (-resis_fr*q_fr + reac_fr*p_fr) >= tan(branch["angmin"])*(w_fr - resis_fr*p_fr + reac_fr*q_fr))
        
        resis_to = (g+g_to)*tm / ((g+g_to)^2+(b+b_to)^2)
        reac_to = -(b+b_to)*tm / ((g+g_to)^2+(b+b_to)^2)

        @constraint(model, (-resis_to*q_to + reac_to*p_to) <= tan(branch["angmax"])*(w_to - resis_to*p_to + reac_to*q_to))
        @constraint(model, (-resis_to*q_to + reac_to*p_to) >= tan(branch["angmin"])*(w_to - resis_to*p_to + reac_to*q_to))

        #@constraint(model, (-r*q_to + x*p_to) <= tan(branch["angmax"])*(w_to - r*p_to + x*q_to))
        #@constraint(model, (-r*q_to + x*p_to) >= tan(branch["angmin"])*(w_to - r*p_to + x*q_to))
        
    end
    
    #for (bp, buspair) in ref[:buspairs], (k,branch) in ref[:branch]
    #    i,j = bp
    #    f_idx = (k, branch["f_bus"], branch["t_bus"])
    #    t_idx = (k, branch["t_bus"], branch["f_bus"])
    #    bp_idx = (branch["f_bus"], branch["t_bus"])

    #    p_fr = p[f_idx]
    #    q_fr = q[f_idx]
    #    p_to = p[t_idx]
    #    q_to = q[t_idx]

    #    w_fr = w[branch["f_bus"]]
    #    w_to = w[branch["t_bus"]]

    #    # Line Flow
    #    g, b = PowerModels.calc_branch_y(branch)
    #    tr, ti = PowerModels.calc_branch_t(branch)
    #    g_fr = branch["g_fr"]
    #    b_fr = branch["b_fr"]
    #    g_to = branch["g_to"]
    #    b_to = branch["b_to"]
    #    tm = branch["tap"]^2

    #    vfub = buspair["vm_fr_max"]
    #    vflb = buspair["vm_fr_min"]
    #    vtub = buspair["vm_to_max"]
    #    vtlb = buspair["vm_to_min"]
    #    tdub = buspair["angmax"]
    #    tdlb = buspair["angmin"]
    #    
    #    resis_to = (g+g_to)*tm / ((g+g_to)^2+(b+b_to)^2)
    #    reac_to = -(b+b_to)*tm / ((g+g_to)^2+(b+b_to)^2)

    #    phi = (tdub + tdlb)/2
    #    d   = (tdub - tdlb)/2

    #    sf = vflb + vfub
    #    st = vtlb + vtub

    #    # Voltage Product Relaxation Upperbound
    #    @constraint(model, sf*st*(cos(phi)*((w[i] - resis_to*p_to + reac_to*q_to)) + sin(phi)*(-resis_to*q_to + reac_to*p_to)) - vtub*cos(d)*st*w[i] - vfub*cos(d)*sf*w[j] >=  vfub*vtub*cos(d)*(vflb*vtlb - vfub*vtub))
    #    @constraint(model, sf*st*(cos(phi)*((w[i] - resis_to*p_to + reac_to*q_to)) + sin(phi)*(-resis_to*q_to + reac_to*p_to)) - vtlb*cos(d)*st*w[i] - vflb*cos(d)*sf*w[j] >= -vflb*vtlb*cos(d)*(vflb*vtlb - vfub*vtub))
    #end

    for (i,dcline) in ref[:dcline]
        # DC Line Flow Constraint
        f_idx = (i, dcline["f_bus"], dcline["t_bus"])
        t_idx = (i, dcline["t_bus"], dcline["f_bus"])
        @constraint(model, (1-dcline["loss1"])*p_dc[f_idx] + (p_dc[t_idx] - dcline["loss0"]) == 0)
    end

    return model
end

model = Model(Ipopt.Optimizer)
build_nf_opf(data, model)
optimize!(model)

