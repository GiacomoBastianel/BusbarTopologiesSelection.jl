### Script analyzing results
# Upload packages
using PowerModelsTopologicalActionsII; const _PMTA = PowerModelsTopologicalActionsII
using PowerModels; const _PM = PowerModels
using PowerModelsACDC; const _PMACDC = PowerModelsACDC
using JuMP, Ipopt, JSON, HiGHS, Gurobi
using PowerPlots, CSV, DataFrames, StatsBase, Plots

####################
# Test case
test_case_folder = joinpath(dirname(@__DIR__),"test_cases")
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-6, "print_level" => 0)

rep = CSV.read(joinpath(test_case_folder,"Daily_K_Means_clustering.csv"),DataFrame)
load = rep.demand
wind_cf = rep.wind_cf
n_timesteps = length(load)


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
test_case["gen"]["138"]["cost"][1] = 9000

test_case_opf = deepcopy(test_case)
#splitted_bus_ac = [49,46]
splitted_bus_ac = [69,24]

n_timesteps = 365
test_case_updated_split = deepcopy(test_case)
test_case_result = deepcopy(test_case)
test_case_updated_split,  switches_couples_result,  extremes_ZILs_result  = _PMTA.AC_busbars_split(test_case_result,splitted_bus_ac)

####################
# Upload results
names = ["49_46_standard","49_46_congested","49_46_standard_data_center","49_46_congested_data_center","69_24_standard","69_24_congested","69_24_standard_data_center","69_24_congested_data_center"]
results_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results"
simulation_name = names[7]

bs_congested = JSON.parsefile(joinpath(results_folder,        "result_bs_$(simulation_name).json"))
opf_congested = JSON.parsefile(joinpath(results_folder,      "result_opf_$(simulation_name).json"))
opf_ac_congested = JSON.parsefile(joinpath(results_folder,"result_opf_ac_$(simulation_name).json"))

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

dict_confs, sorted_confs = sort_configurations(test_case_updated_split,bs_congested,n_timesteps)

#################
# Selecting configurations and timesteps
n_confs_selected = 4
confs_selected = [parse(Int, sorted_confs[i][1]) for i in 1:n_confs_selected]
timeseries_selected = [first(dict_confs["$(confs_selected[i])"]["timesteps"]) for i in 1:n_confs_selected]
conf_and_timeseries = [[confs_selected[i],timeseries_selected[i]] for i in 1:n_confs_selected]
result_congested_comparison = JSON.parsefile(joinpath(results_folder,"comparison_results_$(simulation_name).json"))

result_congested_comparison["1"]

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
number_configurations = count_configurations(n_timesteps,conf_and_timeseries,result_congested_comparison,opf_ac_congested)

###############

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

load_OPF =[load[t] for t in timesteps_OPF]
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

load_combined = vcat(load_4,load_OPF)
wind_combined = vcat(wind_4,load_OPF)


# Add colors here
scatter(load_1, wind_1,  label="Topology 1", color=:blue)
scatter!(load_2, wind_2, label="Topology 2", color=:green)
scatter!(load_3, wind_3, label="Topology 3", color=:orange)
scatter!(load_4, wind_4, label="Topology 4",color=:red)
scatter!(load_OPF,wind_OPF,label="Original",color=:gray,xlabel= "Demand (p.u.)",ylabel="Wind capacity factor [-]",legend=:topright,grid =:none)

#scatter!(load_153 ,wind_153,label="Configuration 4",xlabel= "Demand [-]",ylabel="Wind capacity factor [-]",legend=:topright,grid =:none)
#scatter!(load_OPF,wind_OPF,label="OPF")

results_figures_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results/Figures"
savefig(joinpath(results_figures_folder,"$(simulation_name)_distribution.png"))
savefig(joinpath(results_figures_folder,"$(simulation_name)_distribution.pdf"))
savefig(joinpath(results_figures_folder,"$(simulation_name)_distribution.svg"))


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

