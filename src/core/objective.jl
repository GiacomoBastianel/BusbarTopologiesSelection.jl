function calc_gen_cost_topology(pm::_PM.AbstractPowerModel; nw::Int=_PM.nw_id_default)
    cost = JuMP.AffExpr(0.0)
    for (g_id,g) in pm.ref[:it][:pm][:nw][nw][:gen]
        if length(g["cost"]) ≥ 2
            JuMP.add_to_expression!(cost, g["cost"][end-1], pm.var[:it][:pm][:nw][nw][:pg][g_id])
        end
    end
    return cost
end

function objective_min_gen_cost_topology(pm::_PM.AbstractPowerModel)
    nws = collect(1:length(pm.data["nw"])) # Vector{String}
    #nws = sort(parse.(Int, collect(keys(pm.data["nw"]))))
    println("nws in objective: ", nws)
    JuMP.@objective(pm.model, Min, sum(calc_gen_cost_topology(pm; nw=n) for n in nws))
    return
end
