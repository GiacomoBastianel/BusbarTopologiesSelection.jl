
function split_every_hour(data,n_hours,wind_series,load_series,load_multiplier)
    result = Dict{String,Any}()
    hourly_grid = deepcopy(data)
    for h in 1:n_hours
        println("Processing hour $h")
        hourly_grid["gen"]["30"]["pmax"] = data["gen"]["30"]["pmax"] * wind_series[h]
        for (l_id,l) in hourly_grid["load"]
            hourly_grid["load"][l_id]["pd"] = data["load"][l_id]["pd"] * load_multiplier * load_series[h]
        end
        result_switches_lpac  = _PMTA.run_acdcsw_AC_grid(hourly_grid,LPACCPowerModel,gurobi)
        result["$h"] = result_switches_lpac
    end
    return result
end

function split_opf(data_opf,n_hours,wind_series,load_series,formulation,load_multiplier)
    result_opf = Dict{String,Any}()
    hourly_grid_opf = deepcopy(data_opf)
    for h in 1:n_hours
        hourly_grid_opf["gen"]["30"]["pmax"] = data_opf["gen"]["30"]["pmax"] * wind_series[h]
        for (l_id,l) in hourly_grid_opf["load"]
            hourly_grid_opf["load"][l_id]["pd"] = data_opf["load"][l_id]["pd"] * load_multiplier * load_series[h]
        end
        result_opf_hour = _PM.solve_opf(hourly_grid_opf,formulation,ipopt)
        result_opf["$h"] = result_opf_hour
    end
    return result_opf
end