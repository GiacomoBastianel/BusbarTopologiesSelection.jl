using PowerModelsTopologicalActionsII; const _PMTA = PowerModelsTopologicalActionsII
using PowerModels; const _PM = PowerModels
using PowerModelsACDC; const _PMACDC = PowerModelsACDC
using JuMP, Ipopt, JSON, HiGHS
using Gurobi
using PowerPlots
using CSV
using DataFrames
using StatsBase
using Plots
using BusbarTopologiesSelection; const _BTS = BusbarTopologiesSelection

#test_case_folder = joinpath(dirname(@__DIR__),"test_cases","RTS_GMLC_data","RTS_data")
test_case_folder = joinpath(dirname(@__DIR__),"test_cases")

ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0)
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"MIPGap" => 5e-4, "time_limit" => 180)
highs = JuMP.optimizer_with_attributes(HiGHS.Optimizer)

#test_case = _PM.parse_file(joinpath(test_case_folder,"FormattedData","MATPOWER","RTS_GMLC.m"))
test_case = _PM.parse_file(joinpath(test_case_folder,"pglib_opf_case118_ieee.m"))

function split_every_hour(data,n_hours,wind_series,load_series,load_multiplier)
    result = Dict{String,Any}()
    hourly_grid = deepcopy(data)
    for h in 1:n_hours
        hourly_grid["gen"]["30"]["pmax"] = data["gen"]["30"]["pmax"] * wind_series[h]
        for (l_id,l) in hourly_grid["load"]
            if parse(Int64,l_id) < 100 
                hourly_grid["load"][l_id]["pd"] = data["load"][l_id]["pd"] * load_multiplier * load_series[h]
            end
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
            if parse(Int64,l_id) < 100 
                hourly_grid_opf["load"][l_id]["pd"] = data_opf["load"][l_id]["pd"] * load_multiplier * load_series[h]
            end
        end
        result_opf_hour = _PM.solve_opf(hourly_grid_opf,formulation,ipopt)
        result_opf["$h"] = result_opf_hour
    end
    return result_opf
end

###########
rep = CSV.read(joinpath(test_case_folder,"Daily_K_Means_clustering.csv"),DataFrame)
load = rep.demand
wind_cf = rep.wind_cf
n_timesteps = length(load)


wind_series = CSV.read(joinpath(test_case_folder,"RTS_GMLC_data","RTS_Data","timeseries_data_files","WIND","DAY_AHEAD_wind.csv"),DataFrame)
wind_69 = wind_series[:,7]
maximum(wind_69)
cap_factor_69 = wind_69 ./ maximum(wind_69)

load_series = CSV.read(joinpath(test_case_folder,"RTS_GMLC_data","RTS_Data","timeseries_data_files","Load","DAY_AHEAD_regional_Load.csv"),DataFrame)
load_series_2 = load_series[:,6]
maximum(load_series_2)
cap_factor_load = load_series_2 ./ maximum(load_series_2)
maximum(cap_factor_load)

time_series = DataFrame(
    demand = cap_factor_load,      # MW
    wind_cf = cap_factor_69     # ∈ [0,1]
)


Plots.scatter(time_series.demand, time_series.wind_cf,label = "Original times series",xlabel = "Demand [p.u.]",ylabel = "Wind capacity factor [-]",grid=:none,xlims = (0.25,1))
Plots.scatter!(rep.demand, rep.wind_cf,color = :red,label = "K-means clustering")

results_figures_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results/Figures"
savefig(joinpath(results_figures_folder,"Time_series_k_means_clustering.png"))
savefig(joinpath(results_figures_folder,"Time_series_k_means_clustering.pdf"))
savefig(joinpath(results_figures_folder,"Time_series_k_means_clustering.svg"))



test_case = _PM.parse_file(joinpath(test_case_folder,"pglib_opf_case118_ieee.m"))
test_case_original = _PM.parse_file(joinpath(test_case_folder,"pglib_opf_case118_ieee.m"))

# Adding load
test_case["load"]["100"] = deepcopy(test_case["load"]["99"])
test_case["load"]["100"]["source_id"][2] = 69
test_case["load"]["100"]["load_bus"] = 69
test_case["load"]["100"]["pd"] = deepcopy(test_case["load"]["97"]["pd"])

test_case_original["load"]["100"] = deepcopy(test_case["load"]["99"])
test_case_original["load"]["100"]["source_id"][2] = 69
test_case_original["load"]["100"]["load_bus"] = 69
test_case_original["load"]["100"]["pd"] = deepcopy(test_case["load"]["97"]["pd"])

1.84/sum(test_case_original["load"][l]["pd"] for l in keys(test_case_original["load"]))


for (b_id,b) in test_case["bus"]
    b["vmax"] = 1.1
    b["vmin"] = 0.9
end
for (b_id,b) in test_case_original["bus"]
    b["vmax"] = 1.1
    b["vmin"] = 0.9
end

function add_VOLL_generators(data)
    first_l = maximum(parse.(Int, keys(data["gen"])))
    count = 0
    for (b_id,b) in data["bus"]
        count += 1
        l = first_l + count
        data["gen"]["$l"] = deepcopy(data["gen"]["26"])
        data["gen"]["$l"]["gen_bus"] = parse(Int64,b_id) 
        data["gen"]["$l"]["pmax"] = 99.99
        data["gen"]["$l"]["source_id"][2] = deepcopy(l)
        data["gen"]["$l"]["index"] = l 
        data["gen"]["$l"]["type"] = "VOLL"
        data["gen"]["$l"]["cost"][1] = 10000
        data["gen"]["$l"]["zone"] = b["zone"]
        println("Added VOLL gen $l at bus $b_id")
    end
end
add_VOLL_generators(test_case)

test_case_opf = deepcopy(test_case)
splitted_bus_ac = [49,46]
test_case_updated_split = deepcopy(test_case)



# Preparing data
test_case_updated_split,  switches_couples,  extremes_ZILs  = _PMTA.AC_busbars_split(test_case,splitted_bus_ac)

result_bs_congested = split_every_hour(test_case_updated_split,n_timesteps,wind_cf,load,2)
result_opf_congested = split_opf(test_case_opf,n_timesteps,wind_cf,load,LPACCPowerModel,2)
result_opf_ac_congested = split_opf(test_case_opf,n_timesteps,wind_cf,load,ACPPowerModel,2)


results_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results"

json_result_bs_congested = JSON.json(result_bs_congested)
open(joinpath(results_folder,"result_bs_congested_49_46_data_center.json"),"w") do f 
    write(f, json_result_bs_congested) 
end 

json_result_opf_congested = JSON.json(result_opf_congested)
open(joinpath(results_folder,"result_opf_congested_49_46_data_center.json"),"w") do f 
    write(f, json_result_opf_congested) 
end 

json_result_opf_ac_congested = JSON.json(result_opf_ac_congested)
open(joinpath(results_folder,"result_opf_ac_congested_49_46_data_center.json"),"w") do f 
    write(f, json_result_opf_ac_congested) 
end 

#######################################

obj_bs_congested = [result_bs_congested["$h"]["objective"] for h in 1:n_timesteps]
sol_bs_congested = [result_bs_congested["$h"]["primal_status"] for h in 1:n_timesteps]
countmap(sol_bs_congested)

obj_opf_congested = [result_opf_congested["$h"]["objective"] for h in 1:n_timesteps]
sol_opf_congested = [result_opf_congested["$h"]["primal_status"] for h in 1:n_timesteps]
countmap(sol_opf_congested)


diff_congested = obj_opf_congested - obj_bs_congested
Plots.scatter(diff_congested)


test_case_bs_check = deepcopy(test_case_updated_split)
test_case_bs_check_auxiliary = deepcopy(test_case_updated_split)
_PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs_congested["44"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case_original)
test_case_bs_check["gen"]["30"]["pmax"] = test_case["gen"]["30"]["pmax"] * wind_cf[44]
for (l_id,l) in test_case_bs_check["load"]
    test_case_bs_check["load"][l_id]["pd"] = test_case["load"][l_id]["pd"] * 2 * load[44]
end
results_ac_check_single = _PM.solve_opf(test_case_bs_check,ACPPowerModel,ipopt)
results_lpac_check_single = _PM.solve_opf(test_case_bs_check,LPACCPowerModel,ipopt)

##################

conf_switches = []
for h in 1:n_timesteps
    if result_bs_congested["$h"]["primal_status"] == FEASIBLE_POINT
        count_negative = 0
        conf_switches_hour = []
        for sw_id in 1:length(test_case_updated_split["switch"])
            if result_bs_congested["$h"]["solution"]["switch"]["$sw_id"]["status"] == -0.0
                count_negative += 1
                push!(conf_switches_hour,0.0)
            else
                push!(conf_switches_hour,result_bs_congested["$h"]["solution"]["switch"]["$sw_id"]["status"])
            end
        end
        println("Hour $h, $count_negative")
        push!(conf_switches,conf_switches_hour)
    end
end
map_conf = countmap(conf_switches)
positions_map = Dict{Any, Vector{Int}}()
for (idx, elem) in enumerate(conf_switches)
    if haskey(positions_map, elem)
        push!(positions_map[elem], idx)
    else
        positions_map[elem] = [idx]
    end
end



count_and_positions = Dict(k => (countmap(conf_switches)[k], positions_map[k]) for k in keys(countmap(conf_switches)))
println(count_and_positions)
sorted_countmap = Dict(k => v for (k, v) in sort(collect(countmap(conf_switches)), by = x -> x[1]))
println(sorted_countmap)



confs = Dict{String,Any}()
count_ = 0
for k in keys(map_conf)
    count_ += 1
    confs["$(count_)"] = Dict{String,Any}()
    confs["$(count_)"]["configuration"] = k  
    confs["$(count_)"]["occurrences"] = map_conf[k]
    confs["$(count_)"]["timesteps"] = positions_map[k]
end
confs["30"]["configuration"]
times = [[k,confs["$(k)"]["occurrences"]] for k in keys(confs)]
sort_times = sort(times, by = x -> x[2], rev = true)



#confs_selected = [89,118,214,180] -> no 49_46
confs_selected = [74, 162, 230, 39]
for i in 1:n_timesteps
    sw_i = []
    for sw_id in 1:length(test_case_updated_split["switch"])
        if result_bs_congested["$i"]["primal_status"] == FEASIBLE_POINT
            if result_bs_congested["$i"]["solution"]["switch"]["$sw_id"]["status"] == -0.0
                push!(sw_i,0.0)
            else
                push!(sw_i,result_bs_congested["$i"]["solution"]["switch"]["$sw_id"]["status"])
            end
        end
    end
    for j in 1:length(sort_times)
        if sw_i == confs["$(j)"]["configuration"] && j in confs_selected
            println("Timestep $i matches configuration $j with occurrences $(confs["$(j)"]["occurrences"])")
        end
    end
end
timeseries_selected = [1,92,260,41]

function try_different_topologies(test_case,test_case_split,wind_cf,load)
    results_different_topologies = Dict{String,Any}()
    for t in 1:n_timesteps
        results_different_topologies["$t"] = Dict{String,Any}()
        results_different_topologies["$t"]["AC_OPF"] = Dict{String,Any}()
        results_different_topologies["$t"]["AC_OPF"]["objective"] = deepcopy(result_opf_ac_congested["$t"]["objective"])
        results_different_topologies["$t"]["AC_OPF"]["primal_status"] = deepcopy(result_opf_ac_congested["$t"]["primal_status"])
        results_different_topologies["$t"]["LPAC_OPF"] = Dict{String,Any}()
        results_different_topologies["$t"]["LPAC_OPF"]["objective"] = deepcopy(result_opf_congested["$t"]["objective"])
        results_different_topologies["$t"]["LPAC_OPF"]["primal_status"] = deepcopy(result_opf_congested["$t"]["primal_status"])
        obj_ac_confs = []
        push!(obj_ac_confs,result_opf_ac_congested["$t"]["objective"])
        obj_lpac_confs = []
        push!(obj_ac_confs,result_opf_congested["$t"]["objective"])
        for conf in confs_selected
            test_case_split_conf = deepcopy(test_case_split)
            test_case_bs_check = deepcopy(test_case_split_conf)
            test_case_bs_check_auxiliary = deepcopy(test_case_split_conf)
            results_different_topologies["$t"]["configuration_$(conf)"] = Dict{String,Any}()
            #if conf == 89
            #    _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs_congested["44"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case)
            #    println("89")
            #elseif conf == 118
            #    _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs_congested["3"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case)
            #    println("118")
            #elseif conf == 214
            #    _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs_congested["32"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case)
            #    println("214")
            #elseif conf == 180
            #    _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs_congested["26"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case)
            #    println("180")
            #end
            # Data center
            if conf == 74
                _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs_congested["1"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case)
                println("74")
            elseif conf == 162
                _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs_congested["92"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case)
                println("162")
            elseif conf == 230
                _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs_congested["260"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case)
                println("230")
            elseif conf == 39
                _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs_congested["41"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case)
                println("139")
            end
            test_case_bs_check["gen"]["30"]["pmax"] = test_case["gen"]["30"]["pmax"] * wind_cf[t]
            for (l_id,l) in test_case_bs_check["load"]
                if parse(Int64,l_id) < 100 
                    test_case_bs_check["load"][l_id]["pd"] = test_case["load"][l_id]["pd"] * 2 * load[t]
                end
            end
            results_ac_check_single = deepcopy(_PM.solve_opf(test_case_bs_check,ACPPowerModel,ipopt))
            println("AC OPF objective: ",results_ac_check_single["objective"])
            push!(obj_ac_confs,results_ac_check_single["objective"])
            results_lpac_check_single = deepcopy(_PM.solve_opf(test_case_bs_check,LPACCPowerModel,ipopt))
            println("LPAC OPF objective: ",results_lpac_check_single["objective"])
            push!(obj_lpac_confs,results_lpac_check_single["objective"])
            results_different_topologies["$t"]["configuration_$(conf)"]["AC_OPF"] = Dict{String,Any}()
            results_different_topologies["$t"]["configuration_$(conf)"]["AC_OPF"]["objective"] = deepcopy(results_ac_check_single["objective"])
            results_different_topologies["$t"]["configuration_$(conf)"]["AC_OPF"]["primal_status"] = deepcopy(results_ac_check_single["primal_status"])
            results_different_topologies["$t"]["configuration_$(conf)"]["LPAC_OPF"] = Dict{String,Any}()
            results_different_topologies["$t"]["configuration_$(conf)"]["LPAC_OPF"]["objective"] = deepcopy(results_lpac_check_single["objective"])
            results_different_topologies["$t"]["configuration_$(conf)"]["LPAC_OPF"]["primal_status"] = deepcopy(results_lpac_check_single["primal_status"])
        end
    end
    return results_different_topologies
end
# This to be corrected, now it is fine
result_diff_top = try_different_topologies(test_case_original,test_case_updated_split,wind_cf,load)

json_result_diff_top = JSON.json(result_diff_top)
open(joinpath(results_folder,"result_diff_top_49_46_data_center.json"),"w") do f 
    write(f, json_result_diff_top) 
end 

for t in 1:50
    println("Timestep $t:")
    println(" Original AC OPF obj: $(result_opf_ac_congested["$t"]["objective"]), primal status: $(result_opf_ac_congested["$t"]["primal_status"])")
    for conf in confs_selected
        println("  Configuration $conf: AC OPF obj: $(result_diff_top["$t"]["configuration_$(conf)"]["AC_OPF"]["objective"]), primal status: $(result_diff_top["$t"]["configuration_$(conf)"]["AC_OPF"]["primal_status"])")
    end
end

#########

function split_every_hour_og_load(data,n_hours,wind_series,load_series)
    result = Dict{String,Any}()
    hourly_grid = deepcopy(data)
    for h in 1:n_hours
        hourly_grid["gen"]["30"]["pmax"] = data["gen"]["30"]["pmax"] * wind_series[h]
        for (l_id,l) in hourly_grid["load"]
            if parse(Int64,l_id) < 100 
                hourly_grid["load"][l_id]["pd"] = data["load"][l_id]["pd"] * load_series[h]
            end
        end
        result_switches_lpac  = _PMTA.run_acdcsw_AC_grid(hourly_grid,LPACCPowerModel,gurobi)
        result["$h"] = result_switches_lpac
    end
    return result
end

function split_opf_og_load(data_opf,n_hours,wind_series,load_series,formulation)
    result_opf = Dict{String,Any}()
    hourly_grid_opf = deepcopy(data_opf)
    for h in 1:n_hours
        hourly_grid_opf["gen"]["30"]["pmax"] = data_opf["gen"]["30"]["pmax"] * wind_series[h]
        for (l_id,l) in hourly_grid_opf["load"]
            if parse(Int64,l_id) < 100 
                hourly_grid_opf["load"][l_id]["pd"] = data_opf["load"][l_id]["pd"] * load_series[h]
            end
        end
        result_opf_hour = _PM.solve_opf(hourly_grid_opf,formulation,ipopt)
        result_opf["$h"] = result_opf_hour
    end
    return result_opf
end

result_bs = split_every_hour_og_load(test_case_updated_split,n_timesteps,wind_cf,load)
result_opf = split_opf_og_load(test_case_opf,n_timesteps,wind_cf,load,LPACCPowerModel)
result_opf_ac = split_opf_og_load(test_case_opf,n_timesteps,wind_cf,load,ACPPowerModel)


test_case_bs_check = deepcopy(test_case_updated_split)
test_case_bs_check_auxiliary = deepcopy(test_case_updated_split)
_PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs["44"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case_original)
test_case_bs_check["gen"]["30"]["pmax"] = test_case["gen"]["30"]["pmax"] * wind_cf[44]
for (l_id,l) in test_case_bs_check["load"]
    test_case_bs_check["load"][l_id]["pd"] = test_case["load"][l_id]["pd"] * load[44]
end
results_ac_check_single = _PM.solve_opf(test_case_bs_check,ACPPowerModel,ipopt)
results_lpac_check_single = _PM.solve_opf(test_case_bs_check,LPACCPowerModel,ipopt)

json_result_bs = JSON.json(result_bs)
open(joinpath(results_folder,"result_bs_49_46_data_center.json"),"w") do f 
    write(f, json_result_bs) 
end 

json_result_opf = JSON.json(result_opf)
open(joinpath(results_folder,"result_opf_49_46_data_center.json"),"w") do f 
    write(f, json_result_opf) 
end 

json_result_opf_ac = JSON.json(result_opf_ac)
open(joinpath(results_folder,"result_opf_ac_49_46_data_center.json"),"w") do f 
    write(f, json_result_opf_ac) 
end 


obj_bs_original = [result_bs["$h"]["objective"] for h in 1:n_timesteps]
sol_bs_original = [result_bs["$h"]["primal_status"] for h in 1:n_timesteps]
countmap(sol_bs_original)

obj_opf_original = [result_opf["$h"]["objective"] for h in 1:n_timesteps]
sol_opf_original = [result_opf["$h"]["primal_status"] for h in 1:n_timesteps]
countmap(sol_opf_original)

##################

conf_switches_original = []
for h in 1:n_timesteps
    if result_bs["$h"]["primal_status"] == FEASIBLE_POINT
        count_negative = 0
        conf_switches_hour = []
        for sw_id in 1:length(test_case_updated_split["switch"])
            if result_bs["$h"]["solution"]["switch"]["$sw_id"]["status"] == -0.0
                count_negative += 1
                push!(conf_switches_hour,0.0)
            else
                push!(conf_switches_hour,result_bs["$h"]["solution"]["switch"]["$sw_id"]["status"])
            end
        end
        println("Hour $h, $count_negative")
        push!(conf_switches_original,conf_switches_hour)
    end
end
map_conf_original = countmap(conf_switches_original)
positions_map_original = Dict{Any, Vector{Int}}()
for (idx, elem) in enumerate(conf_switches_original)
    if haskey(positions_map_original, elem)
        push!(positions_map_original[elem], idx)
    else
        positions_map_original[elem] = [idx]
    end
end

confs_original = Dict{String,Any}()
count_ = 0
for k in keys(map_conf_original)
    count_ += 1
    confs_original["$(count_)"] = Dict{String,Any}()
    confs_original["$(count_)"]["configuration"] = k  
    confs_original["$(count_)"]["occurrences"] = map_conf_original[k]
end
times_original = [[k,confs_original["$(k)"]["occurrences"]] for k in keys(confs_original)]
sort_times_original = sort(times_original, by = x -> x[2], rev = true)


confs_selected_original = [18,82,70,46]
for i in 1:n_timesteps
    sw_i = []
    for sw_id in 1:length(test_case_updated_split["switch"])
        if result_bs["$i"]["primal_status"] == FEASIBLE_POINT
            if result_bs["$i"]["solution"]["switch"]["$sw_id"]["status"] == -0.0
                push!(sw_i,0.0)
            else
                push!(sw_i,result_bs["$i"]["solution"]["switch"]["$sw_id"]["status"])
            end
        end
    end
    for j in 1:length(sort_times_original)
        if sw_i == confs_original["$(j)"]["configuration"] && j in confs_selected_original
            println("Timestep $i matches configuration $j with occurrences $(confs["$(j)"]["occurrences"])")
        end
    end
end

timeseries_selected_original = [6,54,63,2]


test_case_bs_check = deepcopy(test_case_updated_split)
test_case_bs_check_auxiliary = deepcopy(test_case_updated_split)
_PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs["3"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case_original)
test_case_bs_check["gen"]["30"]["pmax"] = test_case["gen"]["30"]["pmax"] * wind_cf[3]
for (l_id,l) in test_case_bs_check["load"]
    test_case_bs_check["load"][l_id]["pd"] = test_case["load"][l_id]["pd"] * load[3]
end
results_ac_check_single = _PM.solve_opf(test_case_bs_check,ACPPowerModel,ipopt)
results_lpac_check_single = _PM.solve_opf(test_case_bs_check,LPACCPowerModel,ipopt)



function try_different_topologies_original(test_case,test_case_split,wind_cf,load)
    results_different_topologies = Dict{String,Any}()
    for t in 1:n_timesteps
        results_different_topologies["$t"] = Dict{String,Any}()
        results_different_topologies["$t"]["AC_OPF"] = Dict{String,Any}()
        results_different_topologies["$t"]["AC_OPF"]["objective"] = deepcopy(result_opf_ac["$t"]["objective"])
        results_different_topologies["$t"]["AC_OPF"]["primal_status"] = deepcopy(result_opf_ac["$t"]["primal_status"])
        results_different_topologies["$t"]["LPAC_OPF"] = Dict{String,Any}()
        results_different_topologies["$t"]["LPAC_OPF"]["objective"] = deepcopy(result_opf["$t"]["objective"])
        results_different_topologies["$t"]["LPAC_OPF"]["primal_status"] = deepcopy(result_opf["$t"]["primal_status"])
        obj_ac_confs = []
        push!(obj_ac_confs,result_opf_ac["$t"]["objective"])
        obj_lpac_confs = []
        push!(obj_ac_confs,result_opf["$t"]["objective"])
        for conf in confs_selected_original
            test_case_split_conf = deepcopy(test_case_split)
            test_case_bs_check = deepcopy(test_case_split_conf)
            test_case_bs_check_auxiliary = deepcopy(test_case_split_conf)
            results_different_topologies["$t"]["configuration_$(conf)"] = Dict{String,Any}()
            #if conf == 29
            #    _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs["3"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case)
            #    println("3")
            #elseif conf == 10
            #    _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs["4"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case)
            #    println("4")
            #elseif conf == 73
            #    _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs["2"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case)
            #    println("2")
            #elseif conf == 63
            #    _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs["1"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case)
            #    println("1")
            #end
            if conf == 18
                _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs["6"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case)
                println("18")
            elseif conf == 82
                _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs["54"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case)
                println("82")
            elseif conf == 70
                _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs["63"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case)
                println("70")
            elseif conf == 46
                _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs["2"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,test_case)
                println("46")
            end
            test_case_bs_check["gen"]["30"]["pmax"] = test_case["gen"]["30"]["pmax"] * wind_cf[t]
            for (l_id,l) in test_case_bs_check["load"]
                if parse(Int64,l_id) < 100         
                    test_case_bs_check["load"][l_id]["pd"] = test_case["load"][l_id]["pd"] * load[t]
                end
            end
            results_ac_check_single = deepcopy(_PM.solve_opf(test_case_bs_check,ACPPowerModel,ipopt))
            println("AC OPF objective: ",results_ac_check_single["objective"])
            push!(obj_ac_confs,results_ac_check_single["objective"])
            results_lpac_check_single = deepcopy(_PM.solve_opf(test_case_bs_check,LPACCPowerModel,ipopt))
            println("LPAC OPF objective: ",results_lpac_check_single["objective"])
            push!(obj_lpac_confs,results_lpac_check_single["objective"])
            results_different_topologies["$t"]["configuration_$(conf)"]["AC_OPF"] = Dict{String,Any}()
            results_different_topologies["$t"]["configuration_$(conf)"]["AC_OPF"]["objective"] = deepcopy(results_ac_check_single["objective"])
            results_different_topologies["$t"]["configuration_$(conf)"]["AC_OPF"]["primal_status"] = deepcopy(results_ac_check_single["primal_status"])
            results_different_topologies["$t"]["configuration_$(conf)"]["LPAC_OPF"] = Dict{String,Any}()
            results_different_topologies["$t"]["configuration_$(conf)"]["LPAC_OPF"]["objective"] = deepcopy(results_lpac_check_single["objective"])
            results_different_topologies["$t"]["configuration_$(conf)"]["LPAC_OPF"]["primal_status"] = deepcopy(results_lpac_check_single["primal_status"])
        end
    end
    return results_different_topologies
end
result_diff_top_original = try_different_topologies_original(test_case_original,test_case_updated_split,wind_cf,load)


json_result_diff_top_original = JSON.json(result_diff_top_original)
open(joinpath(results_folder,"result_diff_top_original_49_46_data_center.json"),"w") do f 
    write(f, json_result_diff_top_original) 
end 

for t in 1:20
    println("Timestep $t:")
    println(" Original AC OPF obj: $(result_opf_ac["$t"]["objective"]), primal status: $(result_opf_ac["$t"]["primal_status"])")
    for conf in confs_selected_original
        println("  Configuration $conf: AC OPF obj: $(result_diff_top_original["$t"]["configuration_$(conf)"]["AC_OPF"]["objective"]), primal status: $(result_diff_top_original["$t"]["configuration_$(conf)"]["AC_OPF"]["primal_status"])")
    end
end


count_opf_original = 0
count_bs_original = 0
for t in 1:n_timesteps
    opfs = []
    for conf in confs_selected_original
        push!(opfs,result_diff_top_original["$t"]["configuration_$(conf)"]["AC_OPF"]["objective"])
    end
    min_opf = minimum(opfs)
    if result_opf_ac["$t"]["objective"] < min_opf
        count_opf_original += 1
    else
        count_bs_original += 1
    end
end



count_opf_congested = 0
count_bs_congested = 0
for t in 1:n_timesteps
    opfs = []
    for conf in confs_selected
        push!(opfs,result_diff_top["$t"]["configuration_$(conf)"]["AC_OPF"]["objective"])
    end
    min_opf = minimum(opfs)
    if result_opf_ac_congested["$t"]["objective"] < min_opf
        count_opf_congested += 1
    else
        count_bs_congested += 1
    end
end


count_18 = 0
count_82 = 0
count_70 = 0
count_46 = 0
count_opf = 0

for t in 1:n_timesteps
    opfs = []
    push!(opfs,[result_opf_ac["$t"]["objective"],"OPF"])
    for conf in confs_selected_original
        push!(opfs,[result_diff_top_original["$t"]["configuration_$(conf)"]["AC_OPF"]["objective"],"$conf"])
    end
    #println(opfs)
    opfs_sorted = sort(opfs, by = x -> first(x))
    min_conf = first(opfs_sorted)[2]
    #println("min_conf is $min_conf")
    if min_conf == "18"
        count_18 += 1
    elseif min_conf == "82"
        count_82 += 1
    elseif min_conf == "70"
        count_70 += 1
    elseif min_conf == "46"
        count_46 += 1
    elseif min_conf == "OPF"
        count_opf += 1
    end
end

test_case_original["branch"]["$(test_case_updated_split["switch"]["3"]["original"])"]["f_bus"]

function print_switch_configuration(conf,confs_original,test_case_updated_split,test_case_original)
    for sw_id in 1:length(test_case_updated_split["switch"])
        if haskey(test_case_updated_split["switch"]["$sw_id"],"auxiliary")
            if test_case_updated_split["switch"]["$sw_id"]["auxiliary"] == "branch"
                println("Switch $sw_id between bus $(test_case_updated_split["switch"]["$sw_id"]["f_bus"]) and bus $(test_case_updated_split["switch"]["$sw_id"]["t_bus"]), auxiliary $(test_case_updated_split["switch"]["$sw_id"]["auxiliary"]), from $(test_case_original["branch"]["$(test_case_updated_split["switch"]["$sw_id"]["original"])"]["f_bus"]),to $(test_case_original["branch"]["$(test_case_updated_split["switch"]["$sw_id"]["original"])"]["t_bus"]),status $(confs_original["$conf"]["configuration"][sw_id])")
            else
                println("Switch $sw_id between bus $(test_case_updated_split["switch"]["$sw_id"]["f_bus"]) and bus $(test_case_updated_split["switch"]["$sw_id"]["t_bus"]), auxiliary $(test_case_updated_split["switch"]["$sw_id"]["auxiliary"]), $(test_case_updated_split["switch"]["$sw_id"]["original"]),status $(confs_original["$conf"]["configuration"][sw_id])")        
            end
        else
            println("Switch $sw_id, status $(confs_original["$conf"]["configuration"][sw_id])") 
        end
    end 
end   

print_switch_configuration(18,confs_original,test_case_updated_split,test_case_original)
print_switch_configuration(82,confs_original,test_case_updated_split,test_case_original)
print_switch_configuration(70,confs_original,test_case_updated_split,test_case_original)
print_switch_configuration(46,confs_original,test_case_updated_split,test_case_original)

######################
# Printing the





count_74 = 0
count_162 = 0
count_230 = 0
count_39 = 0
count_opf = 0

for t in 1:n_timesteps
    opfs = []
    push!(opfs,[result_opf_ac_congested["$t"]["objective"],"OPF"])
    for conf in confs_selected
        push!(opfs,[result_diff_top["$t"]["configuration_$(conf)"]["AC_OPF"]["objective"],"$conf"])
    end
    #println(opfs)
    opfs_sorted = sort(opfs, by = x -> first(x))
    min_conf = first(opfs_sorted)[2]
    #println("min_conf is $min_conf")
    if min_conf == "74"
        count_74 += 1
    elseif min_conf == "162"
        count_162 += 1
    elseif min_conf == "230"
        count_230 += 1
    elseif min_conf == "39"
        count_39 += 1
    elseif min_conf == "OPF"
        count_opf += 1
    end
end


function print_switch_configuration(conf,confs_original,test_case_updated_split,test_case_original)
    println("Configuration $conf:")
    for sw_id in 1:length(test_case_updated_split["switch"])
        if haskey(test_case_updated_split["switch"]["$sw_id"],"auxiliary")
            if test_case_updated_split["switch"]["$sw_id"]["auxiliary"] == "branch"
                println("Switch $sw_id between bus $(test_case_updated_split["switch"]["$sw_id"]["f_bus"]) and bus $(test_case_updated_split["switch"]["$sw_id"]["t_bus"]), auxiliary $(test_case_updated_split["switch"]["$sw_id"]["auxiliary"]), from $(test_case_original["branch"]["$(test_case_updated_split["switch"]["$sw_id"]["original"])"]["f_bus"]),to $(test_case_original["branch"]["$(test_case_updated_split["switch"]["$sw_id"]["original"])"]["t_bus"]),status $(confs_original["$conf"]["configuration"][sw_id])")
            else
                println("Switch $sw_id between bus $(test_case_updated_split["switch"]["$sw_id"]["f_bus"]) and bus $(test_case_updated_split["switch"]["$sw_id"]["t_bus"]), auxiliary $(test_case_updated_split["switch"]["$sw_id"]["auxiliary"]), $(test_case_updated_split["switch"]["$sw_id"]["original"]),status $(confs_original["$conf"]["configuration"][sw_id])")        
            end
        else
            println("Switch $sw_id, status $(confs_original["$conf"]["configuration"][sw_id])") 
        end
    end 
end   

print_switch_configuration(74 ,confs,test_case_updated_split,test_case_original)
print_switch_configuration(162,confs,test_case_updated_split,test_case_original)
print_switch_configuration(230,confs,test_case_updated_split,test_case_original)
print_switch_configuration(39 ,confs,test_case_updated_split,test_case_original)






obj_opf_congested
obj_opf_original

obj_bs_congested
obj_bs_original

scatter(obj_bs_congested)
scatter!(obj_bs_original)

scatter(obj_opf_congested)
scatter!(obj_opf_original)

hourly_curt_opf_congested = []
hourly_curt_bs_congested = []
for t in 1:n_timesteps
    if result_bs_congested["$t"]["primal_status"] == FEASIBLE_POINT
        curt_opf = 0
        curt_bs = 0
        for (g_id,g) in test_case_opf["gen"]
            if haskey(g,"type") && g["type"] == "VOLL"
                curt_bs += result_bs_congested["$t"]["solution"]["gen"]["$g_id"]["pg"]
                curt_opf += result_opf_congested["$t"]["solution"]["gen"]["$g_id"]["pg"]
            end
        end
        push!(hourly_curt_opf_congested,curt_opf)
        push!(hourly_curt_bs_congested,curt_bs)
    end
end
scatter(hourly_curt_opf_congested)
scatter!(hourly_curt_bs_congested)

sum(hourly_curt_bs_congested)/sum(hourly_curt_opf_congested)


hourly_curt_opf_69 = []
hourly_curt_bs_69 = []
for t in 1:n_timesteps
    if result_bs_congested["$t"]["primal_status"] == FEASIBLE_POINT
        curt_opf_69 = 0
        curt_bs_69 = 0
        for (g_id,g) in test_case_opf["gen"]
            if haskey(g,"type") && g["type"] == "VOLL" && g["gen_bus"] == 69
                curt_bs_69 += result_bs_congested["$t"]["solution"]["gen"]["$g_id"]["pg"]
                curt_opf_69 += result_opf_congested["$t"]["solution"]["gen"]["$g_id"]["pg"]
            end
        end
        push!(hourly_curt_opf_69,curt_opf_69)
        push!(hourly_curt_bs_69,curt_bs_69)
    end
end
sum(hourly_curt_opf_69)
sum(hourly_curt_bs_69 )