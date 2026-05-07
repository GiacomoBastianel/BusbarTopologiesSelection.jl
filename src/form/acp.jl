function constraint_power_balance_choose_topology(pm::_PM.AbstractACPModel, n::Int, i::Int, bus_arcs, bus_arcs_dc, bus_arcs_sw, bus_gens, bus_storage, bus_pd, bus_qd, bus_gs, bus_bs)
    vm   = _PM.var(pm, n, :vm, i)
    p    = _PM.var(pm, n,    :p) 
    q    = _PM.var(pm, n,    :q) 
    pg   = _PM.var(pm, n,   :pg) 
    qg   = _PM.var(pm, n,   :qg) 
    ps   = _PM.var(pm, n,   :ps) 
    qs   = _PM.var(pm, n,   :qs) 
    psw  = _PM.var(pm, n,  :psw) 
    qsw  = _PM.var(pm, n,  :qsw) 
    p_dc = _PM.var(pm, n, :p_dc) 
    q_dc = _PM.var(pm, n, :q_dc) 
    z_c  = _PM.var(pm, n, :z_c)


    cstr_p = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        + sum(p_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(psw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(pg[g] for g in bus_gens)
        - sum(ps[s] for s in bus_storage)
        - z_c*sum(pd for (i,pd) in bus_pd)
        - sum(gs for (i,gs) in bus_gs)*vm^2
    )

    cstr_q = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        + sum(q_dc[a_dc] for a_dc in bus_arcs_dc)
        + sum(qsw[a_sw] for a_sw in bus_arcs_sw)
        ==
        sum(qg[g] for g in bus_gens)
        - sum(qs[s] for s in bus_storage)
        - z_c*sum(qd for (i,qd) in bus_qd)
        + sum(bs for (i,bs) in bus_bs)*vm^2
    )

    if _IM.report_duals(pm)
        sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end