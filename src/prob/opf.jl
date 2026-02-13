export solve_opf_choose_topology

function solve_opf_choose_topology(data::Dict{String,Any}, model_type::Type, optimizer; kwargs...)
    return _PM.solve_model(data, model_type, optimizer, build_opf_choose_topology; multinetwork = true, solution_processors=[add_z_c_to_solution!], kwargs...)
end

function build_opf_choose_topology(pm::_PM.AbstractPowerModel)
    for (n, network) in _PM.nws(pm)
        println("network: ", n)
        _PM.variable_bus_voltage(pm, nw=n)
        _PM.variable_gen_power(pm, nw=n)
        _PM.variable_branch_power(pm, nw=n)
        _PM.variable_storage_power(pm, nw=n)
        _PM.variable_dcline_power(pm, nw=n)
        _PM.variable_switch_power(pm, nw=n)
    
        variable_chosen_topology(pm, nw=n)

        _PM.constraint_model_voltage(pm, nw=n)

        for i in _PM.ids(pm, :ref_buses, nw=n)
            _PM.constraint_theta_ref(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :bus, nw=n)
            constraint_power_balance_choose_topology(pm, i, nw=n)
        end

        for i in _PM.ids(pm, :branch, nw=n)
            _PM.constraint_ohms_yt_from(pm, i, nw=n)
            _PM.constraint_ohms_yt_to(pm, i, nw=n)
            _PM.constraint_voltage_angle_difference(pm, i, nw=n) #angle difference across transformer and reactor - useful for LPAC if available?
            _PM.constraint_thermal_limit_from(pm, i, nw=n)
            _PM.constraint_thermal_limit_to(pm, i, nw=n)
        end
        for i in _PM.ids(pm, :dcline, nw=n)
            _PM.constraint_dcline_power_losses(pm, i, nw=n)
        end
    end
    sos_binary_constraint(pm)
    objective_min_gen_cost_topology(pm)
end 
