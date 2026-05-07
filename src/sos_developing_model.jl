using PowerModelsTopologicalActionsII; const _PMTA = PowerModelsTopologicalActionsII
using PowerModels; const _PM = PowerModels
using PowerModelsACDC; const _PMACDC = PowerModelsACDC
using JuMP, Ipopt, JSON, HiGHS
using Gurobi
using PowerPlots
using CSV
using DataFrames
using StatsBase
using InfrastructureModels
using Plots
using BusbarTopologiesSelection; const _BTS = BusbarTopologiesSelection
using Juniper

################
test_case_folder = joinpath(dirname(@__DIR__),"test_cases")
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0)
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"MIPGap" => 5e-4, "time_limit" => 180)
juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt, "mip_solver" => gurobi, "time_limit" => 36000)

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

################
test_case_1 = _PM.parse_file(joinpath(test_case_folder,"pglib_opf_case118_ieee.m"))
test_case_original_1 = _PM.parse_file(joinpath(test_case_folder,"pglib_opf_case118_ieee.m"))

test_case_2 = _PM.parse_file(joinpath(test_case_folder,"pglib_opf_case118_ieee.m"))
test_case_original_2 = _PM.parse_file(joinpath(test_case_folder,"pglib_opf_case118_ieee.m"))

for (b_id,b) in test_case_1["bus"]
    b["vmax"] = 1.1
    b["vmin"] = 0.9
end
for (b_id,b) in test_case_original_1["bus"]
    b["vmax"] = 1.1
    b["vmin"] = 0.9
end

for (b_id,b) in test_case_2["bus"]
    b["vmax"] = 1.1
    b["vmin"] = 0.9
end
for (b_id,b) in test_case_original_2["bus"]
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

add_VOLL_generators(test_case_1)
add_VOLL_generators(test_case_2)
test_case_1["gen"]["138"]["cost"][1] = 9000
test_case_2["gen"]["138"]["cost"][1] = 9000

#=
test_case_1["load"]["100"] = deepcopy(test_case_1["load"]["99"])
test_case_1["load"]["100"]["source_id"][2] = 69
test_case_1["load"]["100"]["load_bus"] = 69
test_case_1["load"]["100"]["pd"] = deepcopy(test_case_1["load"]["97"]["pd"])

test_case_2["load"]["100"] = deepcopy(test_case_2["load"]["99"])
test_case_2["load"]["100"]["source_id"][2] = 69
test_case_2["load"]["100"]["load_bus"] = 69
test_case_2["load"]["100"]["pd"] = deepcopy(test_case_2["load"]["97"]["pd"])

test_case_original_1["load"]["100"] = deepcopy(test_case_original_1["load"]["99"])
test_case_original_1["load"]["100"]["source_id"][2] = 69
test_case_original_1["load"]["100"]["load_bus"] = 69
test_case_original_1["load"]["100"]["pd"] = deepcopy(test_case_original_1["load"]["97"]["pd"])

test_case_original_2["load"]["100"] = deepcopy(test_case_original_2["load"]["99"])
test_case_original_2["load"]["100"]["source_id"][2] = 69
test_case_original_2["load"]["100"]["load_bus"] = 69
test_case_original_2["load"]["100"]["pd"] = deepcopy(test_case_original_2["load"]["97"]["pd"])
=#

test_case_opf_1 = deepcopy(test_case_1)
test_case_opf_2 = deepcopy(test_case_2)
splitted_bus_ac = [69,24]
name_file_1 = "69_24_standard"
name_file_2 = "69_24_congested"

test_case_updated_split_1_result = deepcopy(test_case_1)
test_case_updated_split_2_result = deepcopy(test_case_2)
test_case_1_result = deepcopy(test_case_1)
test_case_2_result = deepcopy(test_case_2)

test_case_updated_split_1_result,  switches_couples_1_result,  extremes_ZILs_1_result  = _PMTA.AC_busbars_split(test_case_1_result,splitted_bus_ac)
test_case_updated_split_2_result,  switches_couples_2_result,  extremes_ZILs_2_result  = _PMTA.AC_busbars_split(test_case_2_result,splitted_bus_ac)

###
# Upload results
results_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results"

bs_congested_1 = JSON.parsefile(joinpath(results_folder, "result_bs_$(name_file_1).json"))
bs_congested_2 = JSON.parsefile(joinpath(results_folder, "result_bs_$(name_file_2).json"))

#################
function sort_configurations(grid_bs,result_bs,n_timesteps)
    conf_switches = []
    for h in 1:n_timesteps
        if result_bs["$h"]["primal_status"] == "FEASIBLE_POINT"
            count_negative = 0
            conf_switches_hour = []
            for sw_id in 1:length(grid_bs["switch"])
                if result_bs["$h"]["solution"]["switch"]["$sw_id"]["status"] == -0.0
                    count_negative += 1
                    push!(conf_switches_hour,0.0)
                else
                    push!(conf_switches_hour,result_bs["$h"]["solution"]["switch"]["$sw_id"]["status"])
                end
            end
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

    confs = Dict{String,Any}()
    count_ = 0 
    for k in keys(map_conf)
        count_ += 1
        confs["$(count_)"] = Dict{String,Any}()
        confs["$(count_)"]["configuration"] = k  
        confs["$(count_)"]["occurrences"] = map_conf[k]
        confs["$(count_)"]["timesteps"] = positions_map[k]
    end
    times = [[k,confs["$(k)"]["occurrences"]] for k in keys(confs)]
    sort_times = sort(times, by = x -> x[2], rev = true)
    return confs, sort_times
end
n_timesteps = 365
dict_confs_1, sorted_confs_1 = sort_configurations(test_case_updated_split_1_result,bs_congested_1,n_timesteps)
dict_confs_2, sorted_confs_2 = sort_configurations(test_case_updated_split_2_result,bs_congested_2,n_timesteps)

n_confs_selected_validation_4 = 4

# Og
confs_selected_1_validation_4 = [parse(Int, sorted_confs_1[i][1]) for i in 1:n_confs_selected_validation_4]
timeseries_selected_1_validation_4 = [first(dict_confs_1["$(confs_selected_1_validation_4[i])"]["timesteps"]) for i in 1:n_confs_selected_validation_4]
conf_and_timeseries_1_validation_4 = [[confs_selected_1_validation_4[i],timeseries_selected_1_validation_4[i]] for i in 1:n_confs_selected_validation_4]


# Congested
confs_selected_2_validation_4 = [parse(Int, sorted_confs_2[i][1]) for i in 1:n_confs_selected_validation_4]
timeseries_selected_2_validation_4 = [first(dict_confs_2["$(confs_selected_2_validation_4[i])"]["timesteps"]) for i in 1:n_confs_selected_validation_4]
conf_and_timeseries_2_validation_4 = [[confs_selected_2_validation_4[i],timeseries_selected_2_validation_4[i]] for i in 1:n_confs_selected_validation_4]

################

results_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results/Validation"

function batch_opf_validation_choosing_topology(first_hour,last_hour,size_batch,load_time_series,wind_time_series,optimizer,formulation,results_folder,grid_original,result_bs,load_multiplier,conf_and_timeseries,name_file,splitted_bus_ac)
    for sample in Int64(first_hour/size_batch):Int64(last_hour/size_batch)
        dict = Dict{String,Any}()
        start_idx = size_batch*(sample-1) + 1
        end_idx = size_batch*(sample-1) + 24
        s = Dict("output" => Dict("duals" => false))
        for t in start_idx:end_idx
            println("===================================")
            println("Let's go with timestep $t")
            println("===================================")
            dict["$t"] = Dict{String,Any}()    
            test_case = deepcopy(grid_original)       
            test_case_updated_split = deepcopy(grid_original)
            test_case_updated_split,  switches_couples,  extremes_ZILs  = _PMTA.AC_busbars_split(test_case,splitted_bus_ac)    
            mn_test_case = _PM.replicate(test_case_updated_split, length(conf_and_timeseries))
            mn_test_case_auxiliary = _PM.replicate(test_case_updated_split, length(conf_and_timeseries))
            for conf in 1:length(conf_and_timeseries)    
                _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs["$(conf_and_timeseries[conf][2])"],mn_test_case_auxiliary["nw"]["$conf"],mn_test_case["nw"]["$conf"],switches_couples,extremes_ZILs,grid_original)
                mn_test_case["nw"]["$conf"]["gen"]["30"]["pmax"] = grid_original["gen"]["30"]["pmax"] * wind_time_series[t]
                for (l_id,l) in mn_test_case["nw"]["$conf"]["load"]
                    if parse(Int64,l_id) < 100 
                        mn_test_case["nw"]["$conf"]["load"][l_id]["pd"] = grid_original["load"][l_id]["pd"] * load_multiplier * load_time_series[t]
                    end
                end
            end
            mn_opf_choose_topology = _BTS.solve_opf_choose_topology(mn_test_case, formulation, optimizer)
            dict["$t"] = deepcopy(mn_opf_choose_topology)
        end
        json_opf = JSON.json(dict)        
        open(joinpath(results_folder,"Choose_topology_$(formulation)_$(name_file)_$(start_idx)_$(end_idx)_$(length(conf_and_timeseries))_configurations.json"),"w") do f 
            write(f, json_opf) 
        end
    end
end

batch_opf_validation_choosing_topology(24,720*4,24,time_series.demand,time_series.wind_cf,gurobi,LPACCPowerModel,results_folder,test_case_1,bs_congested_1,1,conf_and_timeseries_1_validation_4,name_file_1,splitted_bus_ac)
batch_opf_validation_choosing_topology(24,720*4,24,time_series.demand,time_series.wind_cf,gurobi,LPACCPowerModel,results_folder,test_case_2,bs_congested_2,2,conf_and_timeseries_2_validation_4,name_file_2,splitted_bus_ac)

optimized_topology =JSON.parsefile(joinpath(dirname(results_folder),"Results_batch","Result_BS_Congested_hours_1_24.json"))
chosen_topology = JSON.parsefile(joinpath(results_folder,"Choose_topology_$(name_file_2)_1_24_4_configurations.json"))



function upload_results_choosing_topology(first_hour,last_hour,size_batch,formulation,results_folder,name_file,conf_and_timeseries)
    results_dict = Dict{String,Any}()
    for sample in Int64(first_hour/size_batch):Int64(last_hour/size_batch)
        dict = Dict{String,Any}()
        start_idx = size_batch*(sample-1) + 1
        end_idx = size_batch*(sample-1) + 24
        file = joinpath(results_folder,"Choose_topology_$(formulation)_$(name_file)_$(start_idx)_$(end_idx)_$(length(conf_and_timeseries))_configurations.json")
        json_opf = JSON.parsefile(file) 
        for t in start_idx:end_idx
            results_dict["$t"] = deepcopy(json_opf["$t"])
        end       
    end
    return results_dict
end

results_standard = upload_results_choosing_topology(24,168,24,LPACCPowerModel,results_folder,name_file_1,conf_and_timeseries_1_validation_4)
results_congested = upload_results_choosing_topology(24,168,24,LPACCPowerModel,results_folder,name_file_2,conf_and_timeseries_2_validation_4)

function count_chosen_topology(results,conf_and_timeseries)
    count_chosen_topology = Dict{String,Int}()
    for i in 1:length(conf_and_timeseries)
        count_chosen_topology["$i"] = 0
    end
    for t in keys(results)
        for i in 1:length(conf_and_timeseries)
            if results["$t"]["solution"]["solution"]["nw"]["$i"]["z_c"] == 1.0
                count_chosen_topology["$i"] += 1
            end
        end
    end
    return count_chosen_topology
end
count_chosen_topology_standard = count_chosen_topology(results_standard,conf_and_timeseries_1_validation_4)
count_chosen_topology_congested = count_chosen_topology(results_congested,conf_and_timeseries_2_validation_4)

#############
BuS_congested_full_year = JSON.parsefile(joinpath(results_folder,"BuS_49_46_congested_24_8784.json"))
BuS_standard_full_year = JSON.parsefile(joinpath(results_folder,"BuS_49_46_standard_24_8784.json")) 
#=
obj_opt_congested = [BuS_congested_full_year["$t"]["objective"] for t in 1:168]
obj_opt_standard = [BuS_standard_full_year["$t"]["objective"] for t in 1:168]

obj_chosen_congested = [results_congested["$t"]["objective"] for t in 1:168]
obj_chosen_standard = [results_standard["$t"]["objective"] for t in 1:168]

time_congested = [results_congested["$t"]["solve_time"] for t in 1:168]
time_standard = [results_standard["$t"]["solve_time"] for t in 1:168]


obj_opt_congested_diff = [obj_chosen_congested[i] - obj_opt_congested[i] for i in 1:length(obj_opt_congested)]
obj_opt_standard_diff = [obj_chosen_standard[i] - obj_opt_standard[i] for i in 1:length(obj_opt_standard)]
=#

function upload_results_BuS(first_hour,last_hour,size_batch,results_folder,name_file)
    results_dict = Dict{String,Any}()
    for sample in Int64(first_hour):Int64(last_hour/size_batch)
        start_idx = size_batch*(sample-1) + 1
        end_idx = size_batch*(sample-1) + 24
        println("Uploading results from hour $start_idx to hour $end_idx")
        file = joinpath(results_folder,"BuS_$(name_file)_$(start_idx)_$(end_idx)_BuS.json")
        json_bus = JSON.parsefile(file) 
        for t in start_idx:end_idx
            results_dict["$t"] = deepcopy(json_bus["$t"])
        end       
    end
    return results_dict
end
results_BuS_folder = joinpath(results_folder)
last_hour_BuS = 168

BuS_results = upload_results_BuS(1,168,24,results_BuS_folder,name_file_1)

time_original = [BuS_results["$t"]["solve_time"] for t in 1:last_hour_BuS]
obj_original = [BuS_results["$t"]["objective"] for t in 1:last_hour_BuS]

time_standard = [results_standard["$t"]["solve_time"] for t in 1:last_hour_BuS]
obj_standard = [results_standard["$t"]["objective"] for t in 1:last_hour_BuS]

obj_diff_standard = [obj_standard[i] - obj_original[i] for i in 1:60]#length(obj_original)]
solve_time_diff_standard = [time_standard[i] - time_original[i] for i in 1:60]#length(time_original)]

for t in 1:168
    if isnothing(obj_original[t])
        println("Hour $t has no solution")
    end
end
obj_original[36]


BuS_results["61"]

BuS_results_original["36"]

for t in 1:168
    println("Hour $t, value ",obj_standard[t])
end