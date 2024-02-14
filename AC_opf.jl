import PowerModels
import Ipopt
import JuMP
import DataFrames
import CSV
using JuMP
using PowerModels
using Ipopt
using DataFrames
using CSV

function build_ac_opf(data::Dict{String,Any}, model=Model())
    @assert !haskey(data, "multinetwork")
    @assert !haskey(data, "conductors")

    PowerModels.standardize_cost_terms!(data, order=2)
    PowerModels.calc_thermal_limits!(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    @variable(model, va[i in keys(ref[:bus])])
    @variable(model, ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], start=1.0)

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

    #@variable(model, p_dc[a in ref[:arcs_dc]])
    #@variable(model, q_dc[a in ref[:arcs_dc]])

    #for (l,dcline) in ref[:dcline]
    #    f_idx = (l, dcline["f_bus"], dcline["t_bus"])
    #    t_idx = (l, dcline["t_bus"], dcline["f_bus"])

    #    JuMP.set_lower_bound(p_dc[f_idx], dcline["pminf"])
    #    JuMP.set_upper_bound(p_dc[f_idx], dcline["pmaxf"])
    #    JuMP.set_lower_bound(q_dc[f_idx], dcline["qminf"])
    #    JuMP.set_upper_bound(q_dc[f_idx], dcline["qmaxf"])

    #    JuMP.set_lower_bound(p_dc[t_idx], dcline["pmint"])
    #    JuMP.set_upper_bound(p_dc[t_idx], dcline["pmaxt"])
    #    JuMP.set_lower_bound(q_dc[t_idx], dcline["qmint"])
    #    JuMP.set_upper_bound(q_dc[t_idx], dcline["qmaxt"])
    #end

    from_idx = Dict(arc[1] => arc for arc in ref[:arcs_from_dc])
    @objective(model, Max, sum(pg[i] for (i,gen) in ref[:gen])+ sum(qg[i] for (i,gen) in ref[:gen]))
    #@objective(model, Min,
    #    sum(gen["cost"][1]*pg[i]^2 + gen["cost"][2]*pg[i] + gen["cost"][3] for (i,gen) in ref[:gen])
    #    #+sum(dcline["cost"][1]*p_dc[from_idx[i]]^2 + dcline["cost"][2]*p_dc[from_idx[i]] + dcline["cost"][3] for (i,dcline) in ref[:dcline])
    #)

    for (i,bus) in ref[:ref_buses]
        # Refrence Bus
        @constraint(model, va[i] == 0)
    end

    for (i,bus) in ref[:bus]
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Bus KCL
        @constraint(model,
            sum(p[a] for a in ref[:bus_arcs][i]) 
            #+ sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) 
            == sum(pg[g] for g in ref[:bus_gens][i]) -
            sum(load["pd"] for load in bus_loads) -
            sum(shunt["gs"] for shunt in bus_shunts)*vm[i]^2
        )

        @constraint(model,
            sum(q[a] for a in ref[:bus_arcs][i])
            #+ sum(q_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i])
            == sum(qg[g] for g in ref[:bus_gens][i]) -
            sum(load["qd"] for load in bus_loads) +
            sum(shunt["bs"] for shunt in bus_shunts)*vm[i]^2
        )
    end

    for (i,branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        t_idx = (i, branch["t_bus"], branch["f_bus"])

        p_fr = p[f_idx]
        q_fr = q[f_idx]
        p_to = p[t_idx]
        q_to = q[t_idx]

        vm_fr = vm[branch["f_bus"]]
        vm_to = vm[branch["t_bus"]]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        # Line Flow
        g, b = PowerModels.calc_branch_y(branch)
        tr, ti = PowerModels.calc_branch_t(branch)
        g_fr = branch["g_fr"]
        b_fr = branch["b_fr"]
        g_to = branch["g_to"]
        b_to = branch["b_to"]
        tm = branch["tap"]^2

        # AC Line Flow Constraints
        @constraint(model, p_fr ==  (g+g_fr)/tm*vm_fr^2 + (-g*tr+b*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-b*tr-g*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )
        @constraint(model, q_fr == -(b+b_fr)/tm*vm_fr^2 - (-b*tr-g*ti)/tm*(vm_fr*vm_to*cos(va_fr-va_to)) + (-g*tr+b*ti)/tm*(vm_fr*vm_to*sin(va_fr-va_to)) )

        @constraint(model, p_to ==  (g+g_to)*vm_to^2 + (-g*tr-b*ti)/tm*(vm_to*vm_fr*cos(va_to-va_fr)) + (-b*tr+g*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )
        @constraint(model, q_to == -(b+b_to)*vm_to^2 - (-b*tr+g*ti)/tm*(vm_to*vm_fr*cos(va_fr-va_to)) + (-g*tr-b*ti)/tm*(vm_to*vm_fr*sin(va_to-va_fr)) )

        # Phase Angle Difference Limit
        @constraint(model, va_fr - va_to <= branch["angmax"])
        @constraint(model, va_fr - va_to >= branch["angmin"])

        # Apparent Power Limit, From and To
        if haskey(branch, "rate_a")
            @constraint(model, p[f_idx]^2 + q[f_idx]^2 <= branch["rate_a"]^2)
            @constraint(model, p[t_idx]^2 + q[t_idx]^2 <= branch["rate_a"]^2)
        end
    end

    #for (i,dcline) in ref[:dcline]
    #    # DC Line Flow Constraint
    #    f_idx = (i, dcline["f_bus"], dcline["t_bus"])
    #    t_idx = (i, dcline["t_bus"], dcline["f_bus"])

    #    @constraint(model, (1-dcline["loss1"])*p_dc[f_idx] + (p_dc[t_idx] - dcline["loss0"]) == 0)
    #end

    return model
end

obj_vals = Vector{Float64}()
term_st  = Vector()
times = Vector{Float64}()
files2 = joinpath.(abspath("sad"), readdir("sad"))

for file in files2
    data = PowerModels.parse_file(file)
    model = Model(Ipopt.Optimizer)
    JuMP.set_optimizer_attribute(model, "linear_solver", "ma27")
    build_ac_opf(data, model)
    optimize!(model)
    time = JuMP.solve_time(model)
    push!(times, time)
    push!(obj_vals, (JuMP.objective_value(model)))
    push!(term_st, (JuMP.termination_status(model)))
end

df = DataFrame(Datasets = readdir("sad"), Objective = obj_vals, Time = times, Termination = term_st)
CSV.write("acopf_max_sad",df)

