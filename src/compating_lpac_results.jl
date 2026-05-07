using PowerModelsTopologicalActionsII; const _PMTA = PowerModelsTopologicalActionsII
using PowerModels; const _PM = PowerModels
using JuMP, Ipopt, JSON, HiGHS
using BusbarTopologiesSelection; const _BTS = BusbarTopologiesSelection

################
test_case_folder = joinpath(dirname(@__DIR__),"test_cases")

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
results_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results/Validation"
n_configurations = 4

function upload_results_choosing_topology(first_hour,last_hour,size_batch,formulation,results_folder,name_file,conf_and_timeseries)
    results_dict = Dict{String,Any}()
    for sample in Int64(first_hour/size_batch):Int64(last_hour/size_batch)
        dict = Dict{String,Any}()
        start_idx = size_batch*(sample-1) + 1
        end_idx = size_batch*(sample-1) + 24
        file = joinpath(results_folder,"Choose_topology_$(formulation)_$(name_file)_$(start_idx)_$(end_idx)_$(conf_and_timeseries)_configurations.json")
        json_opf = JSON.parsefile(file) 
        for t in start_idx:end_idx
            results_dict["$t"] = deepcopy(json_opf["$t"])
        end       
    end
    return results_dict
end

results_standard = upload_results_choosing_topology(24,2448,24,LPACCPowerModel,results_folder,name_file_1,n_configurations)
#results_congested = upload_results_choosing_topology(24,168,24,LPACCPowerModel,results_folder,name_file_2,conf_and_timeseries_2_validation_4)

function count_chosen_topology(results,conf_and_timeseries)
    count_chosen_topology = Dict{String,Int}()
    for i in 1:conf_and_timeseries
        count_chosen_topology["$i"] = 0
    end
    for t in keys(results)
        for i in 1:conf_and_timeseries
            if results["$t"]["solution"]["solution"]["nw"]["$i"]["z_c"] == 1.0
                count_chosen_topology["$i"] += 1
            end
        end
    end
    return count_chosen_topology
end
count_chosen_topology_standard = count_chosen_topology(results_standard,n_configurations)
#count_chosen_topology_congested = count_chosen_topology(results_congested,conf_and_timeseries_2_validation_4)

#############
#BuS_congested_full_year = JSON.parsefile(joinpath(results_folder,"BuS_49_46_congested_24_8784.json"))
#BuS_standard_full_year = JSON.parsefile(joinpath(results_folder,"BuS_49_46_standard_24_8784.json")) 
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
last_hour_BuS = 720

BuS_results = upload_results_BuS(1,last_hour_BuS,24,results_folder,name_file_1)
OPF_results = JSON.parsefile(joinpath(results_folder,"OPF_$(name_file_1)_24_8784.json"))

time_original = [BuS_results["$t"]["solve_time"] for t in 1:last_hour_BuS]
obj_original = [BuS_results["$t"]["objective"] for t in 1:last_hour_BuS]
obj_original_minus_bc = [BuS_results["$t"]["objective"] - sum([test_case_updated_split_1_result["switch"]["$sw_id"]["cost"] * (1 - BuS_results["$t"]["solution"]["switch"]["$sw_id"]["status"]) for sw_id in 1:2]) for t in 1:last_hour_BuS]

time_standard = [results_standard["$t"]["solve_time"] for t in 1:last_hour_BuS]
obj_standard = [results_standard["$t"]["objective"] for t in 1:last_hour_BuS]
obj_opf = [OPF_results["$t"]["objective"] for t in 1:last_hour_BuS]

obj_diff_standard = [obj_standard[i] - obj_original_minus_bc[i] for i in 1:last_hour_BuS]
obj_diff_opf = [obj_opf[i] - obj_original_minus_bc[i] for i in 1:last_hour_BuS]
solve_time_diff_standard = [time_standard[i] - time_original[i] for i in 1:last_hour_BuS]

scatter(obj_diff_standard)
scatter!(obj_diff_opf)

for t in 1:last_hour_BuS
    if isnothing(obj_original[t])
        println("Hour $t has no solution")
    end
end
obj_original[61]
obj_diff_standard[61]

BuS_results["61"]

BuS_results_original["36"]

for t in 1:168
    println("Hour $t, value ",obj_standard[t])
end