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
#test_case_1["gen"]["138"]["cost"][1] = 9000
#test_case_2["gen"]["138"]["cost"][1] = 9000

test_case_opf_1 = deepcopy(test_case_1)
test_case_opf_2 = deepcopy(test_case_2)
splitted_bus_ac = [49,46]
name_file_1 = "49_46_standard"
name_file_2 = "49_46_congested"
#busbars = "49_46"

test_case_updated_split_1_result = deepcopy(test_case_1)
test_case_updated_split_2_result = deepcopy(test_case_2)
test_case_1_result = deepcopy(test_case_1)
test_case_2_result = deepcopy(test_case_2)

test_case_updated_split_1_result,  switches_couples_1_result,  extremes_ZILs_1_result  = _PMTA.AC_busbars_split(test_case_1_result,splitted_bus_ac)
test_case_updated_split_2_result,  switches_couples_2_result,  extremes_ZILs_2_result  = _PMTA.AC_busbars_split(test_case_2_result,splitted_bus_ac)
###
# Upload results
results_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results/Validation"

###############

function upload_batch_opf_choose_topology_results(first_hour,last_hour,size_batch,results_folder,name_file,type,n_configurations)
    dict = Dict{String,Any}()
    for sample in Int64(first_hour/size_batch):Int64(last_hour/size_batch)
        start_idx = size_batch*(sample-1) + 1
        end_idx = size_batch*(sample-1) + 24
        s = Dict("output" => Dict("duals" => false))
        dict_sample = JSON.parsefile(joinpath(results_folder,"Choose_topology_$(name_file)_$(start_idx)_$(end_idx)_$(type).json"))
        for t in start_idx:end_idx
            dict["$t"] = Dict{String,Any}()
            dict["$t"]["termination_status"] = dict_sample["$t"]["termination_status"]
            dict["$t"]["objective"] = dict_sample["$t"]["objective"]
            dict["$t"]["solution"] = Dict{String,Any}()
            for nw in 1:n_configurations
                if dict_sample["$t"]["solution"]["solution"]["nw"]["$nw"]["z_c"] > 0.99
                    dict["$t"]["solution"]["nw"] = "$nw"
                    dict["$t"]["solution"]["bus"] = deepcopy(dict_sample["$t"]["solution"]["nw"]["$nw"]["bus"])
                    dict["$t"]["solution"]["gen"] = deepcopy(dict_sample["$t"]["solution"]["nw"]["$nw"]["gen"])
                    #dict["$t"]["solution"]["branch"] = deepcopy(dict_sample["$t"]["solution"]["nw"]["$nw"]["branch"])
                end
            end
            println("===================================")
            println("Let's go with timestep $t")
            println("===================================")
        end
        empty!(dict_sample)
        GC.gc()  # Trigger garbage collection to free up
    end
    json_dict = JSON.json(dict)        
    open(joinpath(results_folder,"Choose_topology_$(name_file)_$(first_hour)_$(last_hour)_$(type).json"),"w") do f 
        write(f, json_dict) 
    end
    return dict
end

upload_batch_opf_choose_topology_results(24  ,2160,24,results_folder,name_file_1,"4_configurations",4)
upload_batch_opf_choose_topology_results(2160,4320,24,results_folder,name_file_1,"4_configurations",4)
upload_batch_opf_choose_topology_results(4320,6480,24,results_folder,name_file_1,"4_configurations",4)
upload_batch_opf_choose_topology_results(6480,8784,24,results_folder,name_file_1,"4_configurations",4)
upload_batch_opf_choose_topology_results(24,8784,24,results_folder,name_file_1,"4_configurations",4)

upload_batch_opf_choose_topology_results(24  ,2160,24,results_folder,name_file_2,"4_configurations",4)
upload_batch_opf_choose_topology_results(2160,4320,24,results_folder,name_file_2,"4_configurations",4)
upload_batch_opf_choose_topology_results(4320,6480,24,results_folder,name_file_2,"4_configurations",4)
upload_batch_opf_choose_topology_results(6480,8784,24,results_folder,name_file_2,"4_configurations",4)
upload_batch_opf_choose_topology_results(24,8784,24,results_folder,name_file_2,"4_configurations",4)


function upload_batch_opf_validation_results(first_hour,last_hour,size_batch,results_folder,name_file,type,simulation)
    dict = Dict{String,Any}()
    for sample in 1:366
        start_idx = size_batch*(sample-1) + 1
        end_idx = size_batch*(sample-1) + 24
        s = Dict("output" => Dict("duals" => false))
        dict_sample = JSON.parsefile(joinpath(results_folder,"Validation_$(name_file)_$(start_idx)_$(end_idx)_$(type).json"))
        for t in start_idx:end_idx
            dict["$t"] = Dict{String,Any}()
            dict["$t"]["termination_status"] = dict_sample["$t"]["termination_status"]
            dict["$t"]["objective"] = dict_sample["$t"]["objective"]
            println("===================================")
            println("Let's go with timestep $t")
            println("===================================")
        end
        empty!(dict_sample)
        GC.gc()  # Trigger garbage collection to free up
    end
    json_dict = JSON.json(dict)        
    open(joinpath(results_folder,"$(simulation)_$(name_file)_$(first_hour)_$(last_hour).json"),"w") do f 
        write(f, json_dict) 
    end
    return dict
end

function upload_batch_conf_validation_results_name(first_hour,last_hour,size_batch,results_folder,name_file,type)
    dict = Dict{String,Any}()
    for sample in 1:366
        start_idx = size_batch*(sample-1) + 1
        end_idx = size_batch*(sample-1) + 24
        s = Dict("output" => Dict("duals" => false))
        dict_sample = JSON.parsefile(joinpath(results_folder,"Validation_$(name_file)_$(start_idx)_$(end_idx)_$(type).json"))
        for t in start_idx:end_idx
            dict["$t"] = Dict{String,Any}()
            for top in keys(dict_sample["$t"])
                dict["$t"][top] = Dict{String,Any}()
                dict["$t"][top]["termination_status"] = dict_sample["$t"][top]["termination_status"]
                dict["$t"][top]["objective"] = dict_sample["$t"][top]["objective"]
            end
            println("===================================")
            println("Let's go with timestep $t")
            println("===================================")
        end
        empty!(dict_sample)
        GC.gc()  # Trigger garbage collection to free up
    end
    json_dict = JSON.json(dict)        
    open(joinpath(results_folder,"$(type)_standard_$(name_file)_$(first_hour)_$(last_hour).json"),"w") do f 
        write(f, json_dict) 
    end
    return dict
end


function upload_batch_conf_validation_results(first_hour,last_hour,size_batch,results_folder,name_file,type)
    dict = Dict{String,Any}()
    for sample in 1:366
        start_idx = size_batch*(sample-1) + 1
        end_idx = size_batch*(sample-1) + 24
        s = Dict("output" => Dict("duals" => false))
        dict_sample = JSON.parsefile(joinpath(results_folder,"Validation_$(name_file)_$(start_idx)_$(end_idx)_$(type).json"))
        for t in start_idx:end_idx
            dict["$t"] = Dict{String,Any}()
            for top in keys(dict_sample["$t"])
                dict["$t"][top] = Dict{String,Any}()
                dict["$t"][top]["termination_status"] = dict_sample["$t"][top]["termination_status"]
                dict["$t"][top]["objective"] = dict_sample["$t"][top]["objective"]
            end
            println("===================================")
            println("Let's go with timestep $t")
            println("===================================")
        end
        empty!(dict_sample)
        GC.gc()  # Trigger garbage collection to free up
    end
    json_dict = JSON.json(dict)        
    open(joinpath(results_folder,"$(type)_$(name_file)_$(first_hour)_$(last_hour)_no_137.json"),"w") do f 
        write(f, json_dict) 
    end
    return dict
end

upload_batch_opf_validation_results(24  ,8784,24,results_folder,name_file_1,"BuS","BuS")
upload_batch_opf_validation_results(24  ,8784,24,results_folder,name_file_2,"BuS","BuS")


upload_batch_opf_validation_results(24  ,8784,24,results_folder,name_file_1,"OPF")
upload_batch_opf_validation_results(24  ,8784,24,results_folder,name_file_2,"OPF")

upload_batch_conf_validation_results(24  ,8784,24,results_folder,name_file_1,"1_configurations_no_137")
upload_batch_conf_validation_results(24  ,8784,24,results_folder,name_file_1,"4_configurations_no_137")

upload_batch_conf_validation_results(24  ,8784,24,results_folder,name_file_2,"1_configurations")
upload_batch_conf_validation_results(24  ,8784,24,results_folder,name_file_2,"4_configurations")




file_merged = JSON.parsefile(joinpath(results_folder,"1_configurations_69_24_standard_data_center_24_8784.json"))
file_merged["8552"]["217"]
file_merged["8552"]["112"]
file_merged["8552"]["123"]
file_merged["8552"]["103"]
