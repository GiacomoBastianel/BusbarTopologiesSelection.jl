### Script analyzing results
# Upload packages
using PowerModelsTopologicalActionsII; const _PMTA = PowerModelsTopologicalActionsII
using PowerModels; const _PM = PowerModels
using PowerModelsACDC; const _PMACDC = PowerModelsACDC
using JuMP, Ipopt, JSON, HiGHS, Gurobi
using PowerPlots, CSV, DataFrames, StatsBase, Plots

# Load helper data
test_case_folder = joinpath(dirname(@__DIR__),"test_cases")
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0)

rep = CSV.read(joinpath(test_case_folder,"Daily_K_Means_clustering.csv"),DataFrame)
load = rep.demand
wind_cf = rep.wind_cf
n_timesteps = length(load)

wind_series_yearly = CSV.read(joinpath(test_case_folder,"RTS_GMLC_data","RTS_Data","timeseries_data_files","WIND","DAY_AHEAD_wind.csv"),DataFrame)
wind_69_yearly = wind_series_yearly[:,7]
maximum(wind_69_yearly)
cap_factor_69_yearly = wind_69_yearly ./ maximum(wind_69_yearly)

load_series_yearly = CSV.read(joinpath(test_case_folder,"RTS_GMLC_data","RTS_Data","timeseries_data_files","Load","DAY_AHEAD_regional_Load.csv"),DataFrame)
load_series_2_yearly = load_series_yearly[:,6]
maximum(load_series_2_yearly)
cap_factor_load_yearly = load_series_2_yearly ./ maximum(load_series_2_yearly)
maximum(cap_factor_load_yearly)

time_series = DataFrame(
    demand_yearly = cap_factor_load_yearly,      # MW
    wind_cf_yearly = cap_factor_69_yearly     # ∈ [0,1]
)
hours = length(time_series.demand_yearly)


test_case = _PM.parse_file(joinpath(test_case_folder,"pglib_opf_case118_ieee.m"))
test_case_original = _PM.parse_file(joinpath(test_case_folder,"pglib_opf_case118_ieee.m"))


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
splitted_bus_ac = [69,24]
#name_file = "69_24_standard"
name_file = "69_24_congested"
test_case_updated_split = deepcopy(test_case)
test_case_updated_split,  switches_couples,  extremes_ZILs  = _PMTA.AC_busbars_split(test_case,splitted_bus_ac)

###
# Upload results
results_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results"
load_multiplier = 2
#load_multiplier = 1

bs_congested = JSON.parsefile(joinpath(results_folder,        "result_bs_$(name_file).json"))
opf_congested = JSON.parsefile(joinpath(results_folder,      "result_opf_$(name_file).json"))
opf_ac_congested = JSON.parsefile(joinpath(results_folder,"result_opf_ac_$(name_file).json"))

obj_bs_congested     = [bs_congested["$h"]["objective"] for h in 1:n_timesteps] 
obj_opf_congested    = [opf_congested["$h"]["objective"] for h in 1:n_timesteps] 
obj_opf_ac_congested = [opf_ac_congested["$h"]["objective"] for h in 1:n_timesteps] 

# Standard
#conf_4 = JSON.parsefile(joinpath(results_folder,"Validation","Full_results","OPF_selection_4_conf_$(name_file)_1_8784.json"))
#conf_1 = JSON.parsefile(joinpath(results_folder,"Validation","Full_results","OPF_selection_1_conf_$(name_file)_1_8784.json"))

# Congested
conf_4   = JSON.parsefile(joinpath(results_folder,"Validation","Full_results","OPF_selection_4_conf_$(name_file)_1_8784.json"))
conf_4__ = JSON.parsefile(joinpath(results_folder,"Validation","Full_results","4_configurations_$(name_file)_24_8784.json"))
conf_1 = JSON.parsefile(joinpath(results_folder,"Validation","Full_results","OPF_selection_1_conf_$(name_file)_1_8784.json"))

conf_4["1"]  
conf_4__["1"]

# Common
opf_yearly = JSON.parsefile(joinpath(results_folder,"Validation","Full_results","OPF_$(name_file)_24_8784.json"))
#bus_yearly = JSON.parsefile(joinpath(results_folder,"Validation","Full_results","BuS_$(name_file)_24_8784.json"))


#################

function sort_topologies(grid_bs,result_bs,n_timesteps)
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
        confs["$(count_)"]["topology"] = k  
        confs["$(count_)"]["occurrences"] = map_conf[k]
        confs["$(count_)"]["timesteps"] = positions_map[k]
    end
    times = [[k,confs["$(k)"]["occurrences"]] for k in keys(confs)]
    sort_times = sort(times, by = x -> x[2], rev = true)
    return confs, sort_times
end

dict_confs, sorted_confs = sort_topologies(test_case_updated_split,bs_congested,n_timesteps)

#################
# Selecting Topologys and timesteps
n_confs_selected = 4

# Standard
#topologies_selected = [1,5,6,7]# -> standard conditions
#confs_selected = [parse(Int, sorted_confs[i][1]) for i in topologies_selected]

# Congested
confs_selected = [parse(Int, sorted_confs[i][1]) for i in 1:n_confs_selected]

# Common
timeseries_selected = [first(dict_confs["$(confs_selected[i])"]["timesteps"]) for i in 1:n_confs_selected]
conf_and_timeseries = [[confs_selected[i],timeseries_selected[i]] for i in 1:n_confs_selected]
conf_and_timeseries_1 = [[confs_selected[i],timeseries_selected[i]] for i in 1:1]

###################
#=
for i in 1:hours
    conf_vals = [("$c", conf_4["$i"]["$c"]["objective"]) for c in confs_selected]
    push!(conf_vals, ("OPF", opf_yearly["$i"]["AC_OPF"]["objective"]))
    min_val = minimum(v for (_, v) in conf_vals)
    min_confs = [c for (c, v) in conf_vals if v == min_val]
    if length(min_confs) > 1
        println("Timestep $i: tied minimum $min_val — $(min_confs)")
    end
end

function run_batch_opf_different_topologies(first_hour,last_hour,size_batch,results_folder,name_file,conf_and_timeseries,grid_original,grid_bs,result_bs,wind,load)
    for sample in Int64(first_hour):Int64(last_hour/size_batch)
        dict = Dict{String,Any}()
        start_idx = size_batch*(sample-1) + 1
        end_idx = size_batch*(sample-1) + 24
        s = Dict("output" => Dict("duals" => false))
        for t in start_idx:end_idx
            println("===================================")
            println("Let's go with timestep $t")
            println("===================================")
            dict["$t"] = Dict{String,Any}()
            for c in 1:length(conf_and_timeseries)
                dict["$t"]["$(conf_and_timeseries[c][1])"] = Dict{String,Any}()
                test_case_split_conf = deepcopy(grid_bs)
                test_case_bs_check = deepcopy(test_case_split_conf)
                test_case_bs_check_auxiliary = deepcopy(test_case_split_conf)
                _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs["$(conf_and_timeseries[c][2])"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,grid_original)
                test_case_bs_check["gen"]["30"]["pmax"] = grid_original["gen"]["30"]["pmax"] * wind[t]
                for (l_id,l) in test_case_bs_check["load"]
                    if parse(Int64,l_id) < 100 
                        test_case_bs_check["load"][l_id]["pd"] = grid_original["load"][l_id]["pd"] * load_multiplier * load[t]
                    end
                end
                results_ac_check_single = deepcopy(_PM.solve_opf(test_case_bs_check,ACPPowerModel,ipopt))
                dict["$t"]["$(conf_and_timeseries[c][1])"]["objective"] = deepcopy(results_ac_check_single["objective"])
                dict["$t"]["$(conf_and_timeseries[c][1])"]["termination_status"] = deepcopy(results_ac_check_single["termination_status"])
            end
        end
        json_dict = JSON.json(dict)        
        open(joinpath(results_folder,"Results_batch","OPF_selection_$(length(conf_and_timeseries))_conf_$(name_file)_$(start_idx)_$(end_idx).json"),"w") do f 
        write(f, json_dict) 
        end
    end
end

run_batch_opf_different_topologies(1,8784,24,results_folder,name_file,conf_and_timeseries,test_case_opf,test_case_updated_split,bs_congested,time_series.wind_cf_yearly,time_series.demand_yearly)
run_batch_opf_different_topologies(1,8784,24,results_folder,name_file,conf_and_timeseries_1,test_case_opf,test_case_updated_split,bs_congested,time_series.wind_cf_yearly,time_series.demand_yearly)



function upload_batch_opf_choose_topology_results(first_hour,last_hour,size_batch,results_folder,name_file,conf_and_timeseries)
    dict = Dict{String,Any}()
    for sample in Int64(first_hour):Int64(last_hour/size_batch)
        start_idx = size_batch*(sample-1) + 1
        end_idx = size_batch*(sample-1) + 24
        s = Dict("output" => Dict("duals" => false))
        dict_sample = JSON.parsefile(joinpath(results_folder,"OPF_selection_$(length(conf_and_timeseries))_conf_$(name_file)_$(start_idx)_$(end_idx).json"))
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
    open(joinpath(results_folder,"OPF_selection_$(length(conf_and_timeseries))_conf_$(name_file)_$(first_hour)_$(last_hour).json"),"w") do f 
        write(f, json_dict) 
    end
end
results_folder_validation = joinpath(results_folder,"Results_batch")
upload_batch_opf_choose_topology_results(1,8784,24,results_folder_validation,name_file,conf_and_timeseries)
upload_batch_opf_choose_topology_results(1,8784,24,results_folder_validation,name_file,conf_and_timeseries_1)
=#


#=
function run_batch_opf_different_topologies(first_hour,last_hour,size_batch,results_folder,name_file,conf_and_timeseries,grid_original,grid_bs,result_bs,wind,load)
    for sample in Int64(first_hour):Int64(last_hour/size_batch)
        dict = Dict{String,Any}()
        start_idx = size_batch*(sample-1) + 1
        end_idx = size_batch*(sample-1) + 24
        s = Dict("output" => Dict("duals" => false))
        for t in start_idx:end_idx
            println("===================================")
            println("Let's go with timestep $t")
            println("===================================")
            dict["$t"] = Dict{String,Any}()
            for c in 1:length(conf_and_timeseries)
                dict["$t"]["$(conf_and_timeseries[c][1])"] = Dict{String,Any}()
                test_case_split_conf = deepcopy(grid_bs)
                test_case_bs_check = deepcopy(test_case_split_conf)
                test_case_bs_check_auxiliary = deepcopy(test_case_split_conf)
                _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs["$(conf_and_timeseries[c][2])"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,grid_original)
                test_case_bs_check["gen"]["30"]["pmax"] = grid_original["gen"]["30"]["pmax"] * wind[t]
                for (l_id,l) in test_case_bs_check["load"]
                    if parse(Int64,l_id) < 100 
                        test_case_bs_check["load"][l_id]["pd"] = grid_original["load"][l_id]["pd"] * load_multiplier * load[t]
                    end
                end
                results_ac_check_single = deepcopy(_PM.solve_opf(test_case_bs_check,ACPPowerModel,ipopt))
                dict["$t"]["$(conf_and_timeseries[c][1])"]["objective"] = deepcopy(results_ac_check_single["objective"])
                dict["$t"]["$(conf_and_timeseries[c][1])"]["termination_status"] = deepcopy(results_ac_check_single["termination_status"])
            end
        end
        json_dict = JSON.json(dict)        
        open(joinpath(results_folder,"Results_batch","OPF_selection_$(name_file)_$(start_idx)_$(end_idx).json"),"w") do f 
        write(f, json_dict) 
        end
    end
end
run_batch_opf_different_topologies(1,8784,24,results_folder,name_file,conf_and_timeseries,test_case_opf,test_case_updated_split,bs_congested,time_series.wind_cf_yearly,time_series.demand_yearly)

function upload_batch_opf_choose_topology_results(first_hour,last_hour,size_batch,results_folder,name_file)
    dict = Dict{String,Any}()
    for sample in Int64(first_hour):Int64(last_hour/size_batch)
        start_idx = size_batch*(sample-1) + 1
        end_idx = size_batch*(sample-1) + 24
        s = Dict("output" => Dict("duals" => false))
        dict_sample = JSON.parsefile(joinpath(results_folder,"OPF_selection_$(name_file)_$(start_idx)_$(end_idx).json"))
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
    open(joinpath(results_folder,"OPF_selection_$(name_file)_$(first_hour)_$(last_hour).json"),"w") do f 
        write(f, json_dict) 
    end
end
results_folder_validation = joinpath(results_folder,"Results_batch")
upload_batch_opf_choose_topology_results(1,8784,24,results_folder_validation,name_file)

results_folder_batch = joinpath(results_folder,"Results_batch")
function batch_opf(grid,first_hour,last_hour,size_batch,load_series,wind_series,multiplier_load,results_folder,name_file)
    for sample in Int64(first_hour):Int64(last_hour/size_batch)
        result_opf = Dict{String,Any}()
        start_idx = size_batch*(sample-1) + 1
        end_idx = size_batch*(sample-1) + size_batch
        for h in start_idx:end_idx
            hourly_grid = deepcopy(grid)
            hourly_grid["gen"]["30"]["pmax"] = grid["gen"]["30"]["pmax"] * wind_series[h]
            for (l_id,l) in hourly_grid["load"]
                if parse(Int64,l_id) < 100 
                    hourly_grid["load"][l_id]["pd"] = grid["load"][l_id]["pd"] * multiplier_load * load_series[h]
                end
            end
            result_opf["$h"] = Dict{String,Any}()
            result_opf_lpac  = _PM.solve_opf(hourly_grid,LPACCPowerModel,ipopt)
            result_opf_ac  = _PM.solve_opf(hourly_grid,ACPPowerModel,ipopt)
            result_opf["$h"]["LPAC_OPF"] = result_opf_lpac
            result_opf["$h"]["AC_OPF"] = result_opf_ac
        end
        json_opf = JSON.json(result_opf)        
        open(joinpath(results_folder,"Result_OPF_$(name_file)_$(start_idx)_$(end_idx).json"),"w") do f 
            write(f, json_opf) 
        end
    end
end
batch_opf(test_case_opf,1,8784,24,time_series.demand_yearly,time_series.wind_cf_yearly,load_multiplier,results_folder_batch,name_file)


function upload_batch_opf_results(first_hour,last_hour,size_batch,results_folder,name_file)
    dict = Dict{String,Any}()
    for sample in Int64(first_hour):Int64(last_hour/size_batch)
        start_idx = size_batch*(sample-1) + 1
        end_idx = size_batch*(sample-1) + 24
        s = Dict("output" => Dict("duals" => false))
        dict_sample = JSON.parsefile(joinpath(results_folder,"Result_OPF_$(name_file)_$(start_idx)_$(end_idx).json"))
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
    open(joinpath(results_folder,"OPF_$(name_file)_$(first_hour)_$(last_hour).json"),"w") do f 
        write(f, json_dict) 
    end
end
results_folder_validation = joinpath(results_folder,"Results_batch")
upload_batch_opf_results(1,8784,24,results_folder_validation,name_file)

opf_try = JSON.parsefile(joinpath(results_folder_validation,"Result_OPF_69_24_congested_8761_8784.json"))
=#

###################
#=
function compare_different_topologies(grid_original,grid_bs,switches_couples,extremes_ZILs,n_timesteps,opf_ac,opf_lpac,conf_and_timeseries,result_bs,wind,load,load_multiplier)
    result_dict = Dict{String,Any}()
    for t in 1:n_timesteps
        println("Let's go with timestep $t")
        result_dict["$t"] = Dict{String,Any}()
        result_dict["$t"]["AC_OPF"] = Dict{String,Any}()
        result_dict["$t"]["AC_OPF"]["objective"] = deepcopy(opf_ac["$t"]["objective"])
        result_dict["$t"]["AC_OPF"]["primal_status"] = deepcopy(opf_ac["$t"]["primal_status"])
        result_dict["$t"]["LPAC_OPF"] = Dict{String,Any}()
        result_dict["$t"]["LPAC_OPF"]["objective"] = deepcopy(opf_lpac["$t"]["objective"])
        result_dict["$t"]["LPAC_OPF"]["primal_status"] = deepcopy(opf_lpac["$t"]["primal_status"])
        for i in 1:length(conf_and_timeseries)
            conf = conf_and_timeseries[i][1]
            println("conf is $conf") 
            timeseries = conf_and_timeseries[i][2]
            println("timeseries is $timeseries")
            test_case_split_conf = deepcopy(grid_bs)
            test_case_bs_check = deepcopy(test_case_split_conf)
            test_case_bs_check_auxiliary = deepcopy(test_case_split_conf)
            result_dict["$t"]["$(conf)"] = Dict{String,Any}()
            _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs["$timeseries"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,grid_original)
            test_case_bs_check["gen"]["30"]["pmax"] = grid_original["gen"]["30"]["pmax"] * wind[t]
            for (l_id,l) in test_case_bs_check["load"]
                if parse(Int64,l_id) < 100 
                    test_case_bs_check["load"][l_id]["pd"] = grid_original["load"][l_id]["pd"] * load_multiplier * load[t]
                end
            end
            results_ac_check_single = deepcopy(_PM.solve_opf(test_case_bs_check,ACPPowerModel,ipopt))
            results_lpac_check_single = deepcopy(_PM.solve_opf(test_case_bs_check,LPACCPowerModel,ipopt))
            result_dict["$t"]["$(conf)"]["AC_OPF"] = Dict{String,Any}()
            result_dict["$t"]["$(conf)"]["AC_OPF"]["objective"] = deepcopy(results_ac_check_single["objective"])
            result_dict["$t"]["$(conf)"]["AC_OPF"]["primal_status"] = deepcopy(results_ac_check_single["primal_status"])
            result_dict["$t"]["$(conf)"]["LPAC_OPF"] = Dict{String,Any}()
            result_dict["$t"]["$(conf)"]["LPAC_OPF"]["objective"] = deepcopy(results_lpac_check_single["objective"])
            result_dict["$t"]["$(conf)"]["LPAC_OPF"]["primal_status"] = deepcopy(results_lpac_check_single["primal_status"])
        end
    end
    return result_dict
end
result_congested_comparison = compare_different_topologies(test_case_original,test_case_updated_split,switches_couples,extremes_ZILs,n_timesteps,opf_ac_congested,opf_congested,conf_and_timeseries,bs_congested,wind_cf,load,load_multiplier)

json_result_congested_comparison = JSON.json(result_congested_comparison)
open(joinpath(results_folder,"comparison_results_$(name_file)_1_3_5_6.json"),"w") do f 
    write(f, json_result_congested_comparison) 
end 
=#
#result_congested_comparison = JSON.parsefile(joinpath(results_folder,"comparison_results_$(name_file)_1_3_5_6.json"))
#result_congested_comparison = JSON.parsefile(joinpath(results_folder,"comparison_results_$(name_file).json"))

#=
for i in 1:hours
    conf_vals = [("$c", conf_4["$i"]["$c"]["objective"]) for c in confs_selected]
    push!(conf_vals, ("OPF", opf_yearly["$i"]["objective"]))
    min_val = minimum(v for (_, v) in conf_vals)
    min_confs = [c for (c, v) in conf_vals if v == min_val]
    if length(min_confs) > 1
        println("Timestep $i: tied minimum $min_val — $(min_confs)")
    end
end

objs = Dict{String,Any}()
count_opf = 0
count_bus = 0

for i in 1:hours
    println("Timestep $i")
    objs["$i"] = []
    min = minimum([conf_4["$i"]["$c"]["objective"] for c in confs_selected])
    if min < opf_yearly["$i"]["objective"]
        push!(objs["$i"],[min,"BuS"])
        count_bus += 1
    else
        push!(objs["$i"],[opf_yearly["$i"]["objective"],"OPF"])
        count_opf += 1
    end
end
=#


###################
function count_topologies(n_timesteps,conf_and_timeseries,result_comparison,opf_ac,n_confs_selected)
    count_Topologys = Dict{String,Any}()
    count_Topologys["OPF"] = Dict{String,Any}()
    count_Topologys["OPF"]["times"] = 0
    count_Topologys["OPF"]["timesteps"] = []
    for i in 1:n_confs_selected
        conf = conf_and_timeseries[i][1]
        count_Topologys["$conf"] = Dict{String,Any}()
        count_Topologys["$conf"]["times"] = 0
        count_Topologys["$conf"]["timesteps"] = []
    end

    for t in 1:n_timesteps
        opfs = []
        push!(opfs,[opf_ac["$t"]["objective"],"OPF"])
        for i in 1:length(conf_and_timeseries)
            conf = conf_and_timeseries[i][1]
            push!(opfs,[result_comparison["$t"]["$(conf)"]["objective"],"$conf"])
        end
        opfs_sorted = sort(opfs, by = x -> first(x))
        min_conf = first(opfs_sorted)[2]
        count_Topologys["$min_conf"]["times"] += 1
        push!(count_Topologys["$min_conf"]["timesteps"],t)
    end
    return count_Topologys
end
number_topologies_yearly = count_topologies(hours,conf_and_timeseries,conf_4,opf_yearly,n_confs_selected)
number_topologies_1_yearly = count_topologies(hours,conf_and_timeseries_1,conf_1,opf_yearly,length(conf_and_timeseries_1))

timesteps_OPF= number_topologies_yearly["OPF"]["timesteps"]
timesteps_1 =  number_topologies_yearly["$(conf_and_timeseries[1][1])"]["timesteps"]
timesteps_2 =  number_topologies_yearly["$(conf_and_timeseries[2][1])"]["timesteps"]
timesteps_3 =  number_topologies_yearly["$(conf_and_timeseries[3][1])"]["timesteps"]
timesteps_4 =  number_topologies_yearly["$(conf_and_timeseries[4][1])"]["timesteps"]


wind_OPF = [time_series.wind_cf_yearly[t] for t in timesteps_OPF]
wind_1 =   [time_series.wind_cf_yearly[t] for t in timesteps_1]
wind_2 =   [time_series.wind_cf_yearly[t] for t in timesteps_2]
wind_3 =   [time_series.wind_cf_yearly[t] for t in timesteps_3]
wind_4 =   [time_series.wind_cf_yearly[t] for t in timesteps_4]

load_OPF = [time_series.demand_yearly[t] for t in timesteps_OPF]
load_1 = [time_series.demand_yearly[t] for t in timesteps_1]
load_2 = [time_series.demand_yearly[t] for t in timesteps_2]
load_3 = [time_series.demand_yearly[t] for t in timesteps_3]
load_4 = [time_series.demand_yearly[t] for t in timesteps_4]


opf_1 = sum(opf_yearly["$h"]["objective"] for h in timesteps_1)
bus_1 = sum(conf_4["$h"]["$(conf_and_timeseries[1][1])"]["objective"] for h in timesteps_1)
(1 - bus_1/opf_1)*100

opf_2 = sum(opf_yearly["$h"]["objective"] for h in timesteps_2)
bus_2 = sum(conf_4["$h"]["$(conf_and_timeseries[2][1])"]["objective"] for h in timesteps_2)
(1 - bus_2/opf_2)*100

opf_3 = sum(opf_yearly["$h"]["objective"] for h in timesteps_3)
bus_3 = sum(conf_4["$h"]["$(conf_and_timeseries[3][1])"]["objective"] for h in timesteps_3)
(1 - bus_3/opf_3)*100

opf_4 = sum(opf_yearly["$h"]["objective"] for h in timesteps_4)
bus_4 = sum(conf_4["$h"]["$(conf_and_timeseries[4][1])"]["objective"] for h in timesteps_4)
(1 - bus_4/opf_4)*100

# Add colors here
scatter(load_1, wind_1, label="Topology 1", color=:blue)
scatter!(load_2, wind_2, label="Topology 2", color=:green)
scatter!(load_3, wind_3, label="Topology 3", color=:orange)
scatter!(load_4, wind_4, label="Topology 4",color=:red)
scatter!(load_OPF,wind_OPF,label="Original",color=:gray,xlabel= "Demand (p.u.)",ylabel="Wind capacity factor [-]",legend=:topright,grid =:none,xticks=0.3:0.1:1.0)


results_figures_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results/Figures"
savefig(joinpath(results_figures_folder,"$(name_file)_distribution_$(hours)_timesteps.png"))
savefig(joinpath(results_figures_folder,"$(name_file)_distribution_$(hours)_timesteps.pdf"))
savefig(joinpath(results_figures_folder,"$(name_file)_distribution_$(hours)_timesteps.svg"))


##########
# 1 topology
timesteps_OPF_1= number_topologies_1_yearly["OPF"]["timesteps"]
timesteps_1_1 =  number_topologies_1_yearly["$(conf_and_timeseries_1[1][1])"]["timesteps"]

wind_OPF_1 = [time_series.wind_cf_yearly[t] for t in timesteps_OPF_1]
wind_1_1 =   [time_series.wind_cf_yearly[t] for t in timesteps_1_1]

load_OPF_1 = [time_series.demand_yearly[t] for t in timesteps_OPF_1]
load_1_1 = [time_series.demand_yearly[t] for t in timesteps_1_1]


scatter(load_1_1, wind_1_1, label="Topology 1", color=:blue)
scatter!(load_OPF_1,wind_OPF_1,label="Original",color=:gray,xlabel= "Demand (p.u.)",ylabel="Wind capacity factor [-]",legend=:topright,grid =:none,xticks=0.3:0.1:1.0)

results_figures_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results/Figures"
savefig(joinpath(results_figures_folder,"$(name_file)_distribution_$(n_timesteps)_timesteps_yearly_1_configuration.png"))
savefig(joinpath(results_figures_folder,"$(name_file)_distribution_$(n_timesteps)_timesteps_yearly_1_configuration.pdf"))
savefig(joinpath(results_figures_folder,"$(name_file)_distribution_$(n_timesteps)_timesteps_yearly_1_configuration.svg"))

###########
# To fill in the table
obj_1 = [conf_4["$t"]["$(conf_and_timeseries[1][1])"]["objective"] for t in timesteps_1]
obj_ac_1 = [opf_yearly["$t"]["objective"] for t in timesteps_1]

obj_2 = [conf_4["$t"]["$(conf_and_timeseries[2][1])"]["objective"] for t in timesteps_2]
obj_ac_2 = [opf_yearly["$t"]["objective"] for t in timesteps_2]
ts_2 = [opf_yearly["$t"]["termination_status"] for t in timesteps_2]

obj_3 = [conf_4["$t"]["$(conf_and_timeseries[3][1])"]["objective"] for t in timesteps_3]
obj_ac_3 = [opf_yearly["$t"]["objective"] for t in timesteps_3]

obj_4 = [conf_4["$t"]["$(conf_and_timeseries[4][1])"]["objective"] for t in timesteps_4]
obj_ac_4 = [opf_yearly["$t"]["objective"] for t in timesteps_4]
ts_4 = [opf_yearly["$t"]["termination_status"] for t in timesteps_4]

#(1 - (sum(obj_1)+sum(obj_2)+sum(obj_3)+sum(obj_4))/(sum(obj_ac_1)+sum(obj_ac_2)+sum(obj_ac_3)+sum(obj_ac_4)))*100
(1 - (sum(obj_3)+sum(obj_4))/(sum(obj_ac_3)+sum(obj_ac_4)))*100


remaining_timesteps = setdiff(1:hours, timesteps_1)
remaining_timesteps_ = setdiff(remaining_timesteps,timesteps_2)
remaining_timesteps__ = setdiff(remaining_timesteps_,timesteps_3)
remaining_timesteps___ = setdiff(remaining_timesteps__,timesteps_4)
obj_no_bus = [opf_yearly["$t"]["objective"] for t in remaining_timesteps___]

# Congested
(1 - (sum(obj_1)+sum(obj_2)+sum(obj_3)+sum(obj_4)+sum(obj_no_bus))/(sum(obj_ac_1)+sum(obj_ac_2)+sum(obj_ac_3)+sum(obj_ac_4)+sum(obj_no_bus)))*100
(sum(obj_1)+sum(obj_2)+sum(obj_3)+sum(obj_4)+sum(obj_no_bus))*100
(sum(obj_ac_1)+sum(obj_ac_2)+sum(obj_ac_3)+sum(obj_ac_4)+sum(obj_no_bus))*100

# Standard-> not for here
#(1 - (sum(obj_3)+sum(obj_4)+sum(obj_no_bus))/(sum(obj_ac_3)+sum(obj_ac_4)+sum(obj_no_bus)))*100
#(sum(obj_3)+sum(obj_4)+sum(obj_no_bus))*100
#(sum(obj_ac_3)+sum(obj_ac_4)+sum(obj_no_bus))*100

#(1 - (sum(obj_3)+sum(obj_4)+sum(obj_no_bus))/(sum(obj_ac_3)+sum(obj_ac_4)+sum(obj_no_bus)))*100
#fc_no_bus = [opf_ac_congested["$t"]["termination_status"] for t in 1:hours]
#fc_bs_cong = [bs_congested["$t"]["termination_status"] for t in 1:365]
#countmap(fc_bs_cong)

###################
remaining_timesteps_1 = setdiff(1:hours, timesteps_1_1)

opf_1_1 = sum(opf_yearly["$h"]["objective"] for h in timesteps_1_1)
bus_1_1 = sum(conf_1["$h"]["$(conf_and_timeseries_1[1][1])"]["objective"] for h in timesteps_1_1)
opf_1_1_no_bus = sum(opf_yearly["$h"]["objective"] for h in remaining_timesteps_1)


(1 - (bus_1_1+opf_1_1_no_bus)/(opf_1_1+opf_1_1_no_bus))*100
bus_1_1+opf_1_1_no_bus
opf_1_1+opf_1_1_no_bus


###############



