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

test_case_1 = _PM.parse_file(joinpath(test_case_folder,"pglib_opf_case118_ieee.m"))
test_case_original_1 = _PM.parse_file(joinpath(test_case_folder,"pglib_opf_case118_ieee.m"))

test_case_2 = _PM.parse_file(joinpath(test_case_folder,"pglib_opf_case118_ieee.m"))
test_case_original_2 = _PM.parse_file(joinpath(test_case_folder,"pglib_opf_case118_ieee.m"))

# Adding load

test_case_2["load"]["100"] = deepcopy(test_case["load"]["99"])
test_case_2["load"]["100"]["source_id"][2] = 69
test_case_2["load"]["100"]["load_bus"] = 69
test_case_2["load"]["100"]["pd"] = deepcopy(test_case["load"]["97"]["pd"])

test_case_original_2["load"]["100"] = deepcopy(test_case["load"]["99"])
test_case_original_2["load"]["100"]["source_id"][2] = 69
test_case_original_2["load"]["100"]["load_bus"] = 69
test_case_original_2["load"]["100"]["pd"] = deepcopy(test_case["load"]["97"]["pd"])


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

test_case_opf_1 = deepcopy(test_case_1)
test_case_opf_2 = deepcopy(test_case_2)
splitted_bus_ac = [49,46]
name_file_1 = "49_46_congested"
name_file_2 = "49_46_data_center_congested"

test_case_updated_split_1 = deepcopy(test_case_1)
test_case_updated_split_2 = deepcopy(test_case_2)

test_case_updated_split_1,  switches_couples_1,  extremes_ZILs_1  = _PMTA.AC_busbars_split(test_case_1,splitted_bus_ac)
test_case_updated_split_2,  switches_couples_2,  extremes_ZILs_2  = _PMTA.AC_busbars_split(test_case_2,splitted_bus_ac)

###
# Upload results
results_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results"
load_multiplier = 2

bs_congested_1 = JSON.parsefile(joinpath(results_folder,        "result_bs_congested_49_46.json"))
opf_congested_1 = JSON.parsefile(joinpath(results_folder,      "result_opf_congested_49_46.json"))
opf_ac_congested_1 = JSON.parsefile(joinpath(results_folder,"result_opf_ac_congested_49_46.json"))

bs_congested_2 = JSON.parsefile(joinpath(results_folder,        "result_bs_congested_49_46_data_center.json"))
opf_congested_2 = JSON.parsefile(joinpath(results_folder,      "result_opf_congested_49_46_data_center.json"))
opf_ac_congested_2 = JSON.parsefile(joinpath(results_folder,"result_opf_ac_congested_49_46_data_center.json"))


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

dict_confs_1, sorted_confs_1 = sort_configurations(test_case_updated_split_1,bs_congested_1,n_timesteps)
dict_confs_2, sorted_confs_2 = sort_configurations(test_case_updated_split_2,bs_congested_2,n_timesteps)

#################
# Selecting configurations and timesteps
n_confs_selected = 4
confs_selected_1 = [parse(Int, sorted_confs_1[i][1]) for i in 1:n_confs_selected]
timeseries_selected_1 = [first(dict_confs_1["$(confs_selected[i])"]["timesteps"]) for i in 1:n_confs_selected]
conf_and_timeseries_1 = [[confs_selected_1[i],timeseries_selected_1[i]] for i in 1:n_confs_selected]
result_congested_comparison_1 = JSON.parsefile(joinpath(results_folder,"comparison_results_$(name_file_1).json"))

confs_selected_2 = [parse(Int, sorted_confs_2[i][1]) for i in 1:n_confs_selected]
timeseries_selected_2 = [first(dict_confs_2["$(confs_selected[i])"]["timesteps"]) for i in 1:n_confs_selected]
conf_and_timeseries_2 = [[confs_selected_2[i],timeseries_selected_2[i]] for i in 1:n_confs_selected]
result_congested_comparison_2 = JSON.parsefile(joinpath(results_folder,"comparison_results_$(name_file_2).json"))


###############
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

###############
# Tailored part
##########

print_switch_configuration(conf_and_timeseries_1[1][1],dict_confs_1,test_case_updated_split_1,test_case_original)
print_switch_configuration(conf_and_timeseries_1[2][1],dict_confs_1,test_case_updated_split_1,test_case_original)
print_switch_configuration(conf_and_timeseries_1[3][1],dict_confs_1,test_case_updated_split_1,test_case_original)
print_switch_configuration(conf_and_timeseries_1[4][1],dict_confs_1,test_case_updated_split_1,test_case_original)

print_switch_configuration(conf_and_timeseries_2[1][1],dict_confs_2,test_case_updated_split_2,test_case_original)
print_switch_configuration(conf_and_timeseries_2[2][1],dict_confs_2,test_case_updated_split_2,test_case_original)
print_switch_configuration(conf_and_timeseries_2[3][1],dict_confs_2,test_case_updated_split_2,test_case_original)
print_switch_configuration(conf_and_timeseries_2[4][1],dict_confs_2,test_case_updated_split_2,test_case_original)






###########
function push_switch_configuration(test_case_updated_split,conf,confs_original)
    vect_sw = []
    for sw_id in 1:length(test_case_updated_split["switch"])
        push!(vect_sw,confs_original["$conf"]["configuration"][sw_id])
    end
    return vect_sw
end   

sw_1_1 = push_switch_configuration(test_case_updated_split_1,conf_and_timeseries_1[1][1],dict_confs_1)
sw_2_1 = push_switch_configuration(test_case_updated_split_1,conf_and_timeseries_1[2][1],dict_confs_1)
sw_3_1 = push_switch_configuration(test_case_updated_split_1,conf_and_timeseries_1[3][1],dict_confs_1)
sw_4_1 = push_switch_configuration(test_case_updated_split_1,conf_and_timeseries_1[4][1],dict_confs_1)

sw_1_2 = push_switch_configuration(test_case_updated_split_2,conf_and_timeseries_2[1][1],dict_confs_2)
sw_2_2 = push_switch_configuration(test_case_updated_split_2,conf_and_timeseries_2[2][1],dict_confs_2)
sw_3_2 = push_switch_configuration(test_case_updated_split_2,conf_and_timeseries_2[3][1],dict_confs_2)
sw_4_2 = push_switch_configuration(test_case_updated_split_2,conf_and_timeseries_2[4][1],dict_confs_2)

function compare_switch_configurations(sw_1, sw_2)
    differences = []
    for i in 1:length(sw_1)
        if sw_1[i] != sw_2[i]
            push!(differences, (i, sw_1[i], sw_2[i]))
        end
    end
    return differences
end
differences_sw_1_1_sw_1_2 = compare_switch_configurations(sw_1_1, sw_1_2)
differences_sw_1_2_sw_1_2 = compare_switch_configurations(sw_2_1, sw_2_2)
differences_sw_1_3_sw_1_3 = compare_switch_configurations(sw_3_1, sw_3_2)
differences_sw_1_4_sw_1_4 = compare_switch_configurations(sw_4_1, sw_4_2)

println("Differences between sw_1_1 and sw_1_2:")
for (idx, val1, val2) in differences_sw_1_1_sw_1_2
    println("Switch $idx: sw_1_1 = $val1, sw_1_2 = $val2")
end

#print_switch_configuration(conf_and_timeseries[5][1],dict_confs,test_case_updated_split,test_case_original)

# Top left is 4
# Top right is 2
# Bottom left is 4 
# Bottom right is 1

timesteps_OPF = number_configurations["OPF"]["timesteps"]
timesteps_1 =  number_configurations["$(conf_and_timeseries[1][1])"]["timesteps"]
timesteps_2 =  number_configurations["$(conf_and_timeseries[2][1])"]["timesteps"]
timesteps_3 =  number_configurations["$(conf_and_timeseries[3][1])"]["timesteps"]
timesteps_4 =  number_configurations["$(conf_and_timeseries[4][1])"]["timesteps"]

wind_OPF = [wind_cf[t] for t in timesteps_OPF]
wind_1 = [wind_cf[t] for t in timesteps_1]
wind_2 = [wind_cf[t] for t in timesteps_2]
wind_3 = [wind_cf[t] for t in timesteps_3]
wind_4 = [wind_cf[t] for t in timesteps_4]

load_OPF = [load[t] for t in timesteps_OPF]
load_1 = [load[t] for t in timesteps_1]
load_2 = [load[t] for t in timesteps_2]
load_3 = [load[t] for t in timesteps_3]
load_4 = [load[t] for t in timesteps_4]

opf_1 = sum(result_congested_comparison["$h"]["AC_OPF"]["objective"] for h in timesteps_1)
bus_1 = sum(result_congested_comparison["$h"]["$(conf_and_timeseries[1][1])"]["AC_OPF"]["objective"] for h in timesteps_1)
(1 - bus_1/opf_1)*100

opf_2 = sum(result_congested_comparison["$h"]["AC_OPF"]["objective"] for h in timesteps_2)
bus_2 = sum(result_congested_comparison["$h"]["$(conf_and_timeseries[2][1])"]["AC_OPF"]["objective"] for h in timesteps_2)
(1 - bus_2/opf_2)*100

opf_3 = sum(result_congested_comparison["$h"]["AC_OPF"]["objective"] for h in timesteps_3)
bus_3 = sum(result_congested_comparison["$h"]["$(conf_and_timeseries[3][1])"]["AC_OPF"]["objective"] for h in timesteps_3)
(1 - bus_3/opf_3)*100

opf_4 = sum(result_congested_comparison["$h"]["AC_OPF"]["objective"] for h in timesteps_4)
bus_4 = sum(result_congested_comparison["$h"]["$(conf_and_timeseries[4][1])"]["AC_OPF"]["objective"] for h in timesteps_4)
(1 - bus_4/opf_4)*100

load_combined = vcat(load_1,load_OPF,load_3)
wind_combined = vcat(wind_1,wind_OPF,wind_3)

# Add colors here
scatter(load_1, wind_1, label="Configuration 1", color=:blue)
scatter!(load_2, wind_2, label="Configuration 2", color=:green)
scatter!(load_3, wind_3, label="Configuration 3", color=:orange)
scatter!(load_4, wind_4, label="Configuration 4",color=:red)
scatter!(load_OPF,wind_OPF,label="OPF",color=:gray,xlabel= "Demand (p.u.)",ylabel="Wind capacity factor [-]",legend=:topright,grid =:none)

#scatter!(load_153 ,wind_153,label="Configuration 4",xlabel= "Demand [-]",ylabel="Wind capacity factor [-]",legend=:topright,grid =:none)
#scatter!(load_OPF,wind_OPF,label="OPF")

results_figures_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results/Figures"
savefig(joinpath(results_figures_folder,"$(name_file)_distribution.png"))
savefig(joinpath(results_figures_folder,"$(name_file)_distribution.pdf"))
savefig(joinpath(results_figures_folder,"$(name_file)_distribution.svg"))


###############

function opf_different_topologies(result_dict,grid_original,grid_bs,switches_couples,extremes_ZILs,timesteps_conf,conf,opf_ac,result_bs,wind,load,load_multiplier)
    for t in timesteps_conf
        println("Let's go with timestep $t")
        result_dict["$t"] = Dict{String,Any}()
        result_dict["$t"]["AC_OPF"] = deepcopy(opf_ac["$t"])
        for i in 1:length(conf_and_timeseries)
            test_case_split_conf = deepcopy(grid_bs)
            test_case_bs_check = deepcopy(test_case_split_conf)
            test_case_bs_check_auxiliary = deepcopy(test_case_split_conf)
            result_dict["$t"]["$(conf)"] = Dict{String,Any}()
            _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs["$(conf_and_timeseries[conf][2])"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,grid_original)
            test_case_bs_check["gen"]["30"]["pmax"] = grid_original["gen"]["30"]["pmax"] * wind[t]
            for (l_id,l) in test_case_bs_check["load"]
                if parse(Int64,l_id) < 100 
                    test_case_bs_check["load"][l_id]["pd"] = grid_original["load"][l_id]["pd"] * load_multiplier * load[t]
                end
            end
            results_ac_check_single = deepcopy(_PM.solve_opf(test_case_bs_check,ACPPowerModel,ipopt))
            result_dict["$t"]["$(conf_and_timeseries[conf][1])"] = deepcopy(results_ac_check_single)
        end
    end
end
result_different_topologies = Dict{String,Any}()
opf_different_topologies(result_different_topologies,test_case_original,test_case_updated_split,switches_couples,extremes_ZILs,timesteps_1,1,opf_ac_congested,bs_congested,wind_cf,load,load_multiplier)
opf_different_topologies(result_different_topologies,test_case_original,test_case_updated_split,switches_couples,extremes_ZILs,timesteps_2,2,opf_ac_congested,bs_congested,wind_cf,load,load_multiplier)
opf_different_topologies(result_different_topologies,test_case_original,test_case_updated_split,switches_couples,extremes_ZILs,timesteps_3,3,opf_ac_congested,bs_congested,wind_cf,load,load_multiplier)
opf_different_topologies(result_different_topologies,test_case_original,test_case_updated_split,switches_couples,extremes_ZILs,timesteps_4,4,opf_ac_congested,bs_congested,wind_cf,load,load_multiplier)

