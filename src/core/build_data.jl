function make_multinetwork_time_series_scenarios(
    sn_data::Dict{String,Any},n_scenarios,n_hours,hours,hour_simulation,time_series::Dict{String,Any};
    global_keys = ["dim","name","per_unit","source_type","source_version"],
    check_dim::Bool = true,
    )

    mn_data = Dict{String,Any}("nw"=>Dict{String,Any}())
    _FP._add_mn_global_values!(mn_data, sn_data, global_keys)
    #template_nw = _make_template_nw(sn_data, global_keys)
    for hour in 1:length(hours)
        for scenario_idx in 1:n_scenarios
            n = (hour - 1)*n_scenarios + scenario_idx
            mn_data["nw"]["$n"] = deepcopy(sn_data)#_build_nw(template_nw, sn_data, time_series_idx; share_data = true)
            delete!(mn_data["nw"]["$n"],"dim")
            add_hour_scenario_probability(mn_data,hour,scenario_idx,n,time_series)
            for (g_id,g) in mn_data["nw"]["$n"]["gen"]
                mn_data["nw"]["$n"]["gen"][g_id]["pmax"] = time_series["gen"][g_id]["$hour_simulation"]["$scenario_idx"]["pmax_hourly"]
            end
            for (l_id,l) in mn_data["nw"]["$n"]["load"]
                mn_data["nw"]["$n"]["load"][l_id]["pd"] = time_series["load"][l_id]["$hour_simulation"]["$scenario_idx"]["pd"]
            end
        end
    end
    mn_data["scenarios"] = n_scenarios
    mn_data["hours"] = n_hours
    return mn_data
end