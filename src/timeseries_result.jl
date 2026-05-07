using PowerModelsTopologicalActionsII; const _PMTA = PowerModelsTopologicalActionsII
using PowerModels; const _PM = PowerModels
using JuMP, Ipopt, JSON, HiGHS
using BusbarTopologiesSelection; const _BTS = BusbarTopologiesSelection
using CSV, DataFrames

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
splitted_bus_ac = [49,46]
name_file_1 = "49_46_standard"
name_file_2 = "49_46_congested"

test_case_updated_split_1_result = deepcopy(test_case_1)
test_case_updated_split_2_result = deepcopy(test_case_2)
test_case_1_result = deepcopy(test_case_1)
test_case_2_result = deepcopy(test_case_2)

test_case_updated_split_1_result,  switches_couples_1_result,  extremes_ZILs_1_result  = _PMTA.AC_busbars_split(test_case_1_result,splitted_bus_ac)
test_case_updated_split_2_result,  switches_couples_2_result,  extremes_ZILs_2_result  = _PMTA.AC_busbars_split(test_case_2_result,splitted_bus_ac)

results_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results/Validation/Full_results"
#result_chosen_topology_1 = JSON.parsefile(joinpath(results_folder,"Choose_topology_$(name_file_1)_4_configurations.json"))
result_chosen_topology_1 = JSON.parsefile(joinpath(dirname(dirname(results_folder)),"result_opf_ac_$(name_file_1).json"))

result_chosen_topology_1["1"]


BE_grid = JSON.parsefile(joinpath("/Users/giacomobastianel/Downloads/Belgian_grid_data/BE_grid_with_energy_island.json"))