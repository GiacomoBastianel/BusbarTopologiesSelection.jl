### Script analyzing results
# Upload packages
using PowerModelsTopologicalActionsII; const _PMTA = PowerModelsTopologicalActionsII
using PowerModels; const _PM = PowerModels
using PowerModelsACDC; const _PMACDC = PowerModelsACDC
using JuMP, Ipopt, JSON, HiGHS, Gurobi
using PowerPlots, CSV, DataFrames, StatsBase, Plots

################
test_case_folder = joinpath(dirname(@__DIR__),"test_cases")
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0)

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

test_case_opf_1 = deepcopy(test_case_1)
test_case_opf_2 = deepcopy(test_case_2)
splitted_bus_ac = [49,46]
name_file_1 = "49_46_data_center"
name_file_2 = "49_46_congested_data_center"


test_case_updated_split_1_result = deepcopy(test_case_1)
test_case_updated_split_2_result = deepcopy(test_case_2)
test_case_1_result = deepcopy(test_case_1)
test_case_2_result = deepcopy(test_case_2)

test_case_updated_split_1_result,  switches_couples_1_result,  extremes_ZILs_1_result  = _PMTA.AC_busbars_split(test_case_1_result,splitted_bus_ac)
test_case_updated_split_2_result,  switches_couples_2_result,  extremes_ZILs_2_result  = _PMTA.AC_busbars_split(test_case_2_result,splitted_bus_ac)
###
# Upload results
results_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results"

bs_congested_1 = JSON.parsefile(joinpath(results_folder,        "result_bs_49_46_data_center.json"))
#opf_congested_1 = JSON.parsefile(joinpath(results_folder,      "result_opf_49_46_data_center.json"))
#opf_ac_congested_1 = JSON.parsefile(joinpath(results_folder,"result_opf_ac_49_46_data_center.json"))

bs_congested_2 = JSON.parsefile(joinpath(results_folder,        "result_bs_congested_49_46_data_center.json"))
#opf_congested_2 = JSON.parsefile(joinpath(results_folder,      "result_opf_congested_49_46_data_center.json"))
#opf_ac_congested_2 = JSON.parsefile(joinpath(results_folder,"result_opf_ac_congested_49_46_data_center.json"))

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

#################
# Selecting configurations and timesteps
n_confs_selected_validation_4 = 4
n_confs_selected_validation_1 = 1

# Og
confs_selected_1_validation_4 = [parse(Int, sorted_confs_1[i][1]) for i in 1:n_confs_selected_validation_4]
timeseries_selected_1_validation_4 = [first(dict_confs_1["$(confs_selected_1_validation_4[i])"]["timesteps"]) for i in 1:n_confs_selected_validation_4]
conf_and_timeseries_1_validation_4 = [[confs_selected_1_validation_4[i],timeseries_selected_1_validation_4[i]] for i in 1:n_confs_selected_validation_4]
#result_congested_comparison_1 = JSON.parsefile(joinpath(results_folder,"comparison_results_$(name_file_1).json"))

confs_selected_1_validation_1 = [parse(Int, sorted_confs_1[i][1]) for i in 1:n_confs_selected_validation_1]
timeseries_selected_1_validation_1 = [first(dict_confs_1["$(confs_selected_1_validation_1[i])"]["timesteps"]) for i in 1:n_confs_selected_validation_1]
conf_and_timeseries_1_validation_1 = [[confs_selected_1_validation_1[i],timeseries_selected_1_validation_1[i]] for i in 1:n_confs_selected_validation_1]


# Congested
confs_selected_2_validation_4 = [parse(Int, sorted_confs_2[i][1]) for i in 1:n_confs_selected_validation_4]
timeseries_selected_2_validation_4 = [first(dict_confs_2["$(confs_selected_2_validation_4[i])"]["timesteps"]) for i in 1:n_confs_selected_validation_4]
conf_and_timeseries_2_validation_4 = [[confs_selected_2_validation_4[i],timeseries_selected_2_validation_4[i]] for i in 1:n_confs_selected_validation_4]

confs_selected_2_validation_1 = [parse(Int, sorted_confs_2[i][1]) for i in 1:n_confs_selected_validation_1]
timeseries_selected_2_validation_1 = [first(dict_confs_2["$(confs_selected_2_validation_4[i])"]["timesteps"]) for i in 1:n_confs_selected_validation_1]
conf_and_timeseries_2_validation_1 = [[confs_selected_2_validation_1[i],timeseries_selected_2_validation_1[i]] for i in 1:n_confs_selected_validation_1]


#=
function count_configurations(n_timesteps,conf_and_timeseries,result_comparison,opf_ac)
    count_configurations = Dict{String,Any}()
    count_configurations["OPF"] = Dict{String,Any}()
    count_configurations["OPF"]["times"] = 0
    count_configurations["OPF"]["timesteps"] = []
    for i in 1:n_confs_selected
        conf = conf_and_timeseries[i][1]
        count_configurations["$conf"] = Dict{String,Any}()
        count_configurations["$conf"]["times"] = 0
        count_configurations["$conf"]["timesteps"] = []
    end

    for t in 1:n_timesteps
        opfs = []
        push!(opfs,[opf_ac["$t"]["objective"],"OPF"])
        for i in 1:length(conf_and_timeseries)
            conf = conf_and_timeseries[i][1]
            push!(opfs,[result_comparison["$t"]["$(conf)"]["AC_OPF"]["objective"],"$conf"])
        end
        opfs_sorted = sort(opfs, by = x -> first(x))
        min_conf = first(opfs_sorted)[2]
        count_configurations["$min_conf"]["times"] += 1
        push!(count_configurations["$min_conf"]["timesteps"],t)
    end
    return count_configurations
end
number_configurations_1 = count_configurations(n_timesteps,conf_and_timeseries_1,result_congested_comparison_1,opf_ac_congested_1)
number_configurations_2 = count_configurations(n_timesteps,conf_and_timeseries_2,result_congested_comparison_2,opf_ac_congested_2)
=#

###############

results_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results/Validation"

function upload_batch_opf_validation_results(first_hour,last_hour,size_batch,results_folder,name_file,type)
    dict = Dict{String,Any}()
    for sample in Int64(first_hour/size_batch):Int64(last_hour/size_batch)
        start_idx = size_batch*(sample-1) + 1
        end_idx = size_batch*(sample-1) + 24
        s = Dict("output" => Dict("duals" => false))
        dict_sample = JSON.parsefile(joinpath(results_folder,"Validation_$(name_file)_$(start_idx)_$(end_idx)_$(type).json"))
        for t in start_idx:end_idx
            println("===================================")
            println("Let's go with timestep $t")
            println("===================================")
            dict["$t"] = deepcopy(dict_sample["$t"])
        end
        empty!(dict_sample)
        GC.gc()  # Trigger garbage collection to free up
    end
    json_dict = JSON.json(dict)        
    open(joinpath(results_folder,"Validation_$(name_file)_$(Int64(first_hour/size_batch))_$(last_hour)_$(type).json"),"w") do f 
        write(f, json_dict) 
    end
    return dict
end

#dict = upload_batch_opf_validation_results(24,8784,24,results_folder,conf_and_timeseries_2_validation_4,name_file_2,"OPF")
#dict_opf = JSON.parsefile(joinpath(results_folder,"Validation_$(name_file_2)_1_8784_OPF.json"))

dict_4_configurations = upload_batch_opf_validation_results(24,8784,24,results_folder,name_file_1,"4_configurations")
dict_1_configurations = upload_batch_opf_validation_results(24,8784,24,results_folder,name_file_1,"1_configurations")
dict_1_OPF = upload_batch_opf_validation_results(24,8784,24,results_folder,name_file_1,"OPF")

dict_2_4_configurations = upload_batch_opf_validation_results(24,8784,24,results_folder,name_file_2,"4_configurations")
dict_2_1_configurations = upload_batch_opf_validation_results(24,8784,24,results_folder,name_file_2,"1_configurations")
dict_2_1_OPF = upload_batch_opf_validation_results(24,8784,24,results_folder,name_file_2,"OPF")
