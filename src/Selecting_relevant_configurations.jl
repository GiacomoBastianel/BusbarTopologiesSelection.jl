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
n_timesteps = 365

################
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
#test_case["gen"]["138"]["cost"][1] = 9000

test_case["load"]["100"] = deepcopy(test_case["load"]["99"])
test_case["load"]["100"]["source_id"][2] = 69
test_case["load"]["100"]["load_bus"] = 69
test_case["load"]["100"]["pd"] = deepcopy(test_case["load"]["97"]["pd"])

test_case_opf = deepcopy(test_case)
splitted_bus_ac = [49,46]
name_file = "49_46_standard"


test_case_updated_split_result = deepcopy(test_case)
test_case_result = deepcopy(test_case)

test_case_updated_split_result,  switches_couples_result,  extremes_ZILs_result  = _PMTA.AC_busbars_split(test_case_result,splitted_bus_ac)

###
# Upload results
results_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results"

bs_congested = JSON.parsefile(joinpath(results_folder,"result_bs_$(name_file).json"))
opf_lpac_congested = JSON.parsefile(joinpath(results_folder, "result_opf_$(name_file).json"))

diff_ = []
diff_hourly = []
for t in 1:365 
    if bs_congested["$t"]["primal_status"] == "FEASIBLE_POINT"
        diff_t = opf_lpac_congested["$t"]["objective"] - bs_congested["$t"]["objective"]
        push!(diff_,[diff_t,t])
        push!(diff_hourly,diff_t)
    end
end
findmax(diff_)
scatter(diff_hourly)
sort(diff_,by = x -> x[1], rev = true)

bs_congested["66"]
bs_congested["342"]
bs_congested["20"]


[bs_congested["$t"]["termination_status"] for t in 1:n_timesteps]
countmap([bs_congested["$t"]["termination_status"] for t in 1:n_timesteps])


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
dict_confs, sorted_confs = sort_configurations(test_case_updated_split_result,bs_congested,n_timesteps)


n_confs_selected = 5
confs_selected = [parse(Int, sorted_confs[i][1]) for i in 1:n_confs_selected]
timeseries_selected = [first(dict_confs["$(confs_selected[i])"]["timesteps"]) for i in 1:n_confs_selected]
conf_and_timeseries = [[confs_selected[i],timeseries_selected[i]] for i in 1:n_confs_selected]


function print_switch_configuration(conf,confs_original,grid_bs,grid_original)
    for sw_id in 1:length(grid_bs["switch"])
        if haskey(grid_bs["switch"]["$sw_id"],"auxiliary")
            if grid_bs["switch"]["$sw_id"]["auxiliary"] == "branch"
                println("Switch $sw_id between bus $(grid_bs["switch"]["$sw_id"]["f_bus"]) and bus $(grid_bs["switch"]["$sw_id"]["t_bus"]), auxiliary $(grid_bs["switch"]["$sw_id"]["auxiliary"]), from $(grid_original["branch"]["$(grid_bs["switch"]["$sw_id"]["original"])"]["f_bus"]),to $(grid_original["branch"]["$(grid_bs["switch"]["$sw_id"]["original"])"]["t_bus"]),status $(confs_original["$conf"]["configuration"][sw_id])")
            else
                println("Switch $sw_id between bus $(grid_bs["switch"]["$sw_id"]["f_bus"]) and bus $(grid_bs["switch"]["$sw_id"]["t_bus"]), auxiliary $(grid_bs["switch"]["$sw_id"]["auxiliary"]), $(grid_bs["switch"]["$sw_id"]["original"]),status $(confs_original["$conf"]["configuration"][sw_id])")        
            end
        else
            println("Switch $sw_id, f_bus $(grid_bs["switch"]["$sw_id"]["f_bus"]), t_bus $(grid_bs["switch"]["$sw_id"]["t_bus"]) status $(confs_original["$conf"]["configuration"][sw_id])") 
        end
    end 
end   

print_switch_configuration(conf_and_timeseries[1][1],dict_confs,test_case_updated_split_result,test_case_original)
print_switch_configuration(conf_and_timeseries[2][1],dict_confs,test_case_updated_split_result,test_case_original)
print_switch_configuration(conf_and_timeseries[3][1],dict_confs,test_case_updated_split_result,test_case_original)
print_switch_configuration(conf_and_timeseries[4][1],dict_confs,test_case_updated_split_result,test_case_original)
print_switch_configuration(conf_and_timeseries[5][1],dict_confs,test_case_updated_split_result,test_case_original)


count_switches = Dict{String,Any}()
for t in 1:n_timesteps
    count_switches["$t"] = Dict{String,Any}()
    for busbar in eachindex(extremes_ZILs_result)
        for extreme in extremes_ZILs_result[busbar]
            count_switches["$t"]["$extreme"] = 0
        end
    end
    for (sw_id,sw) in test_case_updated_split_result["switch"]
        if haskey(sw,"auxiliary")
            bs_congested["$t"]["solution"]["switch"]["$sw_id"]["status"] >= 0.9 ? count_switches["$t"]["$(sw["t_bus"])"] += 1 : nothing
        end    
    end
end



timesteps_2 = dict_confs["$(conf_and_timeseries[2][1])"]["timesteps"]
timesteps_3 = dict_confs["$(conf_and_timeseries[3][1])"]["timesteps"]
timesteps_4 = dict_confs["$(conf_and_timeseries[4][1])"]["timesteps"]
timesteps_5 = dict_confs["$(conf_and_timeseries[5][1])"]["timesteps"]

diff_timesteps_2 = [diff_hourly[t] for t in timesteps_2]
diff_timesteps_3 = [diff_hourly[t] for t in timesteps_3]
diff_timesteps_4 = [diff_hourly[t] for t in timesteps_4]
diff_timesteps_5 = [diff_hourly[t] for t in timesteps_5]


count_switches["$(conf_and_timeseries[2][1])"]
count_switches["$(conf_and_timeseries[3][1])"]
count_switches["$(conf_and_timeseries[4][1])"]
count_switches["$(conf_and_timeseries[5][1])"]

#################
# Selecting configurations and timesteps
n_confs_selected_validation_4 = 4
n_confs_selected_validation_1 = 1

# Og
confs_selected_1_validation_4 = [parse(Int, sorted_confs[i][1]) for i in 2:n_confs_selected_validation_4+1]
timeseries_selected_1_validation_4 = [first(dict_confs["$(confs_selected_1_validation_4[i])"]["timesteps"]) for i in 1:n_confs_selected_validation_4]
conf_and_timeseries_1_validation_4 = [[confs_selected_1_validation_4[i],timeseries_selected_1_validation_4[i]] for i in 1:n_confs_selected_validation_4]
#result_congested_comparison_1 = JSON.parsefile(joinpath(results_folder,"comparison_results_$(name_file_1).json"))

confs_selected_1_validation_1 = [parse(Int, sorted_confs[i][1]) for i in n_confs_selected_validation_1:n_confs_selected_validation_1]
timeseries_selected_1_validation_1 = [first(dict_confs["$(confs_selected_1_validation_1[i])"]["timesteps"]) for i in n_confs_selected_validation_1:n_confs_selected_validation_1]
conf_and_timeseries_1_validation_1 = [[confs_selected_1_validation_1[i],timeseries_selected_1_validation_1[i]] for i in n_confs_selected_validation_1:n_confs_selected_validation_1]







###############

results_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results/Validation"

function batch_opf_validation_opf(first_hour,last_hour,size_batch,load_time_series,wind_time_series,optimizer,formulation,results_folder,grid_original,load_multiplier,name_file)
    for sample in Int64(first_hour/size_batch):Int64(last_hour/size_batch)
        dict = Dict{String,Any}()
        start_idx = size_batch*(sample-1) + 1
        end_idx = size_batch*(sample-1) + 24
        s = Dict("output" => Dict("duals" => false))
        for t in start_idx:end_idx
            println("===================================")
            println("Let's go with timestep $t")
            println("===================================")
            test_case = deepcopy(grid_original)  
            test_case["gen"]["30"]["pmax"] = grid_original["gen"]["30"]["pmax"] * wind_time_series[t]
            for (l_id,l) in test_case["load"]
                if parse(Int64,l_id) < 100 
                    test_case["load"][l_id]["pd"] = grid_original["load"][l_id]["pd"] * load_multiplier * load_time_series[t]
                end
            end
            dict["$t"] = deepcopy(_PM.solve_opf(test_case,formulation,optimizer))
        end
        json_opf = JSON.json(dict)        
        open(joinpath(results_folder,"Validation_$(name_file)_$(start_idx)_$(end_idx)_OPF.json"),"w") do f 
            write(f, json_opf) 
        end
    end
end

try_batch_og = batch_opf_validation_opf_optimized_topology(24,8784,24,time_series.demand,time_series.wind_cf,ipopt,ACPPowerModel,results_folder,test_case_1,bs_congested_1,1,conf_and_timeseries_1_validation_4,name_file_1,splitted_bus_ac)


try_batch_opf = batch_opf_validation_opf(24,8784,24,time_series.demand,time_series.wind_cf,ipopt,ACPPowerModel,results_folder,test_case_opf_1,1,name_file_1)





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
    open(joinpath(results_folder,"Validation_$(name_file)_$(Int64(first_hour/size_batch))_$(last_hour)_$(type)_no_173.json"),"w") do f 
        write(f, json_dict) 
    end
    return dict
end

upload_batch_opf_validation_results(24,8784,24,results_folder,name_file,"4_configurations")
