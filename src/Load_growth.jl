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
load    = rep.demand
wind_cf = rep.wind_cf
weights = rep.weight
load_wind = [[load[i], wind_cf[i]] for i in 1:length(load)]
findmax(weights)
sorted_weights = sort(weights, rev = true)
sorted_load = sort(load, rev = true)
sorted_wind_cf = sort(wind_cf, rev = true)
sorted_load_wind_cf = sort(load_wind, by = x -> x[1], rev = true)

position = findfirst(x -> x == sorted_weights[3], weights)
position = findfirst(x -> x == 0.8632142294736842, load)
position = findfirst(x -> x ==  0.9923359757125989, wind_cf)

selected_timesteps = [235, 65, 234, 247]
load[247]
wind_cf[247]

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
name_file_1 = "49_46"
name_file_2 = "49_46_congested"

test_case_updated_split_1 = deepcopy(test_case_1)
test_case_updated_split_2 = deepcopy(test_case_2)

test_case_updated_split_1,  switches_couples_1,  extremes_ZILs_1  = _PMTA.AC_busbars_split(test_case_1,splitted_bus_ac)
test_case_updated_split_2,  switches_couples_2,  extremes_ZILs_2  = _PMTA.AC_busbars_split(test_case_2,splitted_bus_ac)

###
# Upload results
results_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results"
load_multiplier = 2

bs_congested_1 = JSON.parsefile(joinpath(results_folder,        "result_bs_49_46.json"))
opf_congested_1 = JSON.parsefile(joinpath(results_folder,      "result_opf_49_46.json"))
opf_ac_congested_1 = JSON.parsefile(joinpath(results_folder,"result_opf_ac_49_46.json"))

bs_congested_2 = JSON.parsefile(joinpath(results_folder,        "result_bs_congested_49_46.json"))
opf_congested_2 = JSON.parsefile(joinpath(results_folder,      "result_opf_congested_49_46.json"))
opf_ac_congested_2 = JSON.parsefile(joinpath(results_folder,"result_opf_ac_congested_49_46.json"))

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
# Og
confs_selected_1 = [parse(Int, sorted_confs_1[i][1]) for i in 1:n_confs_selected]
timeseries_selected_1 = [first(dict_confs_1["$(confs_selected[i])"]["timesteps"]) for i in 1:n_confs_selected]
conf_and_timeseries_1 = [[confs_selected_1[i],timeseries_selected_1[i]] for i in 1:n_confs_selected]
result_congested_comparison_1 = JSON.parsefile(joinpath(results_folder,"comparison_results_$(name_file_1).json"))

# Congested
confs_selected_2 = [parse(Int, sorted_confs_2[i][1]) for i in 1:n_confs_selected]
timeseries_selected_2 = [first(dict_confs_2["$(confs_selected[i])"]["timesteps"]) for i in 1:n_confs_selected]
conf_and_timeseries_2 = [[confs_selected_2[i],timeseries_selected_2[i]] for i in 1:n_confs_selected]
result_congested_comparison_2 = JSON.parsefile(joinpath(results_folder,"comparison_results_$(name_file_2).json"))

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
selected_timesteps = [235, 65, 234, 247]  
# 235 -> Low load, low wind  
# 65  -> High load, low wind 
# 234 -> Low load, high wind 
# 247 -> High load, high wind

load[247]
wind_cf[247]
weights[247]

###############

function print_switch_configuration(conf,confs_original,grid_bs,grid_original)
    for sw_id in 1:length(test_case_updated_split["switch"])
        if haskey(test_case_updated_split["switch"]["$sw_id"],"auxiliary")
            if test_case_updated_split["switch"]["$sw_id"]["auxiliary"] == "branch"
                println("Switch $sw_id between bus $(grid_bs["switch"]["$sw_id"]["f_bus"]) and bus $(grid_bs["switch"]["$sw_id"]["t_bus"]), auxiliary $(grid_bs["switch"]["$sw_id"]["auxiliary"]), from $(grid_original["branch"]["$(grid_bs["switch"]["$sw_id"]["original"])"]["f_bus"]),to $(grid_original["branch"]["$(grid_bs["switch"]["$sw_id"]["original"])"]["t_bus"]),status $(confs_original["$conf"]["configuration"][sw_id])")
            else
                println("Switch $sw_id between bus $(grid_bs["switch"]["$sw_id"]["f_bus"]) and bus $(grid_bs["switch"]["$sw_id"]["t_bus"]), auxiliary $(grid_bs["switch"]["$sw_id"]["auxiliary"]), $(grid_bs["switch"]["$sw_id"]["original"]),status $(confs_original["$conf"]["configuration"][sw_id])")        
            end
        else
            println("Switch $sw_id, f_bus $(grid_bs["switch"]["$sw_id"]["f_bus"]), t_bus $(grid_bs["switch"]["$sw_id"]["t_bus"]) status $(confs_original["$conf"]["configuration"][sw_id])") 
        end
    end 
end   

print_switch_configuration(conf_and_timeseries_1[1][1],dict_confs_1,test_case_updated_split_1,test_case_original_1)
print_switch_configuration(conf_and_timeseries_1[2][1],dict_confs_1,test_case_updated_split_1,test_case_original_1)
print_switch_configuration(conf_and_timeseries_1[3][1],dict_confs_1,test_case_updated_split_1,test_case_original_1)
print_switch_configuration(conf_and_timeseries_1[4][1],dict_confs_1,test_case_updated_split_1,test_case_original_1)


#################
# Adding load
test_case_opf_1["load"]["100"] = deepcopy(test_case["load"]["99"])
test_case_opf_1["load"]["100"]["source_id"][2] = 69
test_case_opf_1["load"]["100"]["load_bus"] = 69
test_case_opf_1["load"]["100"]["pd"] = deepcopy(test_case["load"]["97"]["pd"])

test_case_opf_2["load"]["100"] = deepcopy(test_case["load"]["99"])
test_case_opf_2["load"]["100"]["source_id"][2] = 69
test_case_opf_2["load"]["100"]["load_bus"] = 69
test_case_opf_2["load"]["100"]["pd"] = deepcopy(test_case["load"]["97"]["pd"])

test_case_updated_split_1["load"]["100"] = deepcopy(test_case["load"]["99"])
test_case_updated_split_1["load"]["100"]["source_id"][2] = 69
test_case_updated_split_1["load"]["100"]["load_bus"] = 69
test_case_updated_split_1["load"]["100"]["pd"] = deepcopy(test_case["load"]["97"]["pd"])

test_case_updated_split_2["load"]["100"] = deepcopy(test_case["load"]["99"])
test_case_updated_split_2["load"]["100"]["source_id"][2] = 69
test_case_updated_split_2["load"]["100"]["load_bus"] = 69
test_case_updated_split_2["load"]["100"]["pd"] = deepcopy(test_case["load"]["97"]["pd"])



growing_load = collect(0.0:0.2:18.4) # five times of the load (50% max load)
function opf_load_growth_selected_topology(result_dict,grid_original,grid_bs,switches_couples,extremes_ZILs,timestep,conf,result_bs,wind,load,load_multiplier,growing_load)
    result_dict["$timestep"] = Dict{String,Any}()
    for gl in growing_load
        result_dict["$timestep"]["$gl"] = Dict{String,Any}()
        test_case_split_conf = deepcopy(grid_bs)
        test_case_bs_check = deepcopy(test_case_split_conf)
        test_case_bs_check_auxiliary = deepcopy(test_case_split_conf)
        _PMTA.prepare_AC_feasibility_check_AC_busbars(result_bs["$(conf_and_timeseries[conf][2])"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples,extremes_ZILs,grid_original)
        test_case_bs_check["gen"]["30"]["pmax"] = grid_original["gen"]["30"]["pmax"] * wind[timestep]
        for (l_id,l) in test_case_bs_check["load"]
            if parse(Int64,l_id) < 100 
                test_case_bs_check["load"][l_id]["pd"] = grid_original["load"][l_id]["pd"] * load_multiplier * load[timestep]
            else
                test_case_bs_check["load"][l_id]["pd"] = deepcopy(gl)
                println("Load at bus $(l["load_bus"]) set to $(gl) MW")
                println("-----------")
            end
        end
        result_dict["$timestep"]["$gl"] = deepcopy(_PM.solve_opf(test_case_bs_check,ACPPowerModel,ipopt))
    end
end
function opf_load_growth_opf(result_dict,grid_original,timestep,wind,load,load_multiplier,growing_load)
    result_dict["$timestep"] = Dict{String,Any}()
    for gl in growing_load
        println("-----------")
        println("gl $gl")
        result_dict["$timestep"]["$gl"] = Dict{String,Any}()
        test_case = deepcopy(grid_original)
        test_case["gen"]["30"]["pmax"] = grid_original["gen"]["30"]["pmax"] * wind[timestep]
        for (l_id,l) in test_case["load"]
            if parse(Int64,l_id) < 100 
                test_case["load"][l_id]["pd"] = grid_original["load"][l_id]["pd"] * load_multiplier * load[timestep]
            elseif parse(Int64,l_id) == 100
                test_case["load"][l_id]["pd"] = deepcopy(gl)
                println("Load at bus $(l["load_bus"]) set to $(gl) MW")
                println("-----------")
            end
        end
        result_dict["$timestep"]["$gl"] = deepcopy(_PM.solve_opf(test_case,ACPPowerModel,ipopt))
    end
end

#-> OPF
#growing_load = collect(100:0.2:104) # five times of the load (50% max load)
load_multiplier_1 = 1
result_og = Dict{String,Any}()
#opf_load_growth_opf(result_og,test_case_opf_1,235,wind_cf,load,load_multiplier_1,growing_load)
opf_load_growth_selected_topology(result_og,test_case_opf_1,test_case_updated_split_1,switches_couples_1,extremes_ZILs_1,235, 1,bs_congested_1,wind_cf,load,load_multiplier_1,growing_load)
opf_load_growth_selected_topology(result_og,test_case_opf_1,test_case_updated_split_1,switches_couples_1,extremes_ZILs_1,65, 1,bs_congested_1,wind_cf,load,load_multiplier_1,growing_load)
#opf_load_growth_opf(result_og,test_case_opf_1,234,wind_cf,load,load_multiplier_1,growing_load)
opf_load_growth_selected_topology(result_og,test_case_opf_1,test_case_updated_split_1,switches_couples_1,extremes_ZILs_1,234,1,bs_congested_1,wind_cf,load,load_multiplier_1,growing_load)
opf_load_growth_selected_topology(result_og,test_case_opf_1,test_case_updated_split_1,switches_couples_1,extremes_ZILs_1,247,1,bs_congested_1,wind_cf,load,load_multiplier_1,growing_load)

load_multiplier_1 = 1
result_og_opf = Dict{String,Any}()
opf_load_growth_opf(result_og_opf,test_case_opf_1,235,wind_cf,load,load_multiplier_1,growing_load)
opf_load_growth_opf(result_og_opf,test_case_opf_1,65 ,wind_cf,load,load_multiplier_1,growing_load)
opf_load_growth_opf(result_og_opf,test_case_opf_1,234,wind_cf,load,load_multiplier_1,growing_load)
opf_load_growth_opf(result_og_opf,test_case_opf_1,247,wind_cf,load,load_multiplier_1,growing_load)


og_235      = [result_og["235"]["$gl"]["objective"] for gl in growing_load]
og_235_curt = [result_og["235"]["$gl"]["solution"]["gen"]["138"]["pg"] for gl in growing_load]
og_opf_235_curt = [result_og_opf["235"]["$gl"]["solution"]["gen"]["138"]["pg"] for gl in growing_load]

og_65      = [result_og["65"]["$gl"]["objective"] for gl in growing_load]
og_65_curt = [result_og["65"]["$gl"]["solution"]["gen"]["138"]["pg"] for gl in growing_load]
og_opf_65_curt = [result_og_opf["65"]["$gl"]["solution"]["gen"]["138"]["pg"] for gl in growing_load]

og_234      = [result_og["234"]["$gl"]["objective"] for gl in growing_load]
og_234_curt = [result_og["234"]["$gl"]["solution"]["gen"]["138"]["pg"] for gl in growing_load]
og_opf_234_curt = [result_og_opf["234"]["$gl"]["solution"]["gen"]["138"]["pg"] for gl in growing_load]

og_247      = [result_og["247"]["$gl"]["objective"] for gl in growing_load]
og_247_curt = [result_og["247"]["$gl"]["solution"]["gen"]["138"]["pg"] for gl in growing_load]
og_opf_247_curt = [result_og_opf["247"]["$gl"]["solution"]["gen"]["138"]["pg"] for gl in growing_load]


Plots.plot(growing_load/test_case_opf_2["load"]["100"]["pd"],og_235_curt./test_case_opf_2["load"]["100"]["pd"],label = "Low load, low wind",grid=:none,xlims = (-0.05,8.15),xticks = 0:0.5:8,yticks = 0:1:8,ylims=(-0.05,8.05),ylabel = "Curtailment / Biggest load in the test case [-]",xlabel = "Load growth bus 69 / Biggest load in the test case [-]")#,ylabelfontsize = 8,xlabelfontsize = 8,xtickfontsize = 6,ytickfontsize = 6)
Plots.plot!(growing_load/test_case_opf_2["load"]["100"]["pd"], og_65_curt./test_case_opf_2["load"]["100"]["pd"],label = "High load, low wind")
Plots.plot!(growing_load/test_case_opf_2["load"]["100"]["pd"],og_234_curt./test_case_opf_2["load"]["100"]["pd"],label = "Low load, high wind")
Plots.plot!(growing_load/test_case_opf_2["load"]["100"]["pd"],og_247_curt./test_case_opf_2["load"]["100"]["pd"],label = "High load, high wind")
Plots.plot!(growing_load/test_case_opf_2["load"]["100"]["pd"],og_opf_235_curt./test_case_opf_2["load"]["100"]["pd"],label = "Low load, low wind (OPF)")
Plots.plot!(growing_load/test_case_opf_2["load"]["100"]["pd"], og_opf_65_curt./test_case_opf_2["load"]["100"]["pd"],label = "High load, low wind (OPF)")
Plots.plot!(growing_load/test_case_opf_2["load"]["100"]["pd"],og_opf_234_curt./test_case_opf_2["load"]["100"]["pd"],label = "Low load, low wind (OPF)")
Plots.plot!(growing_load/test_case_opf_2["load"]["100"]["pd"], og_opf_247_curt./test_case_opf_2["load"]["100"]["pd"],label = "High load, low wind (OPF)")


results_figures_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results/Figures"
savefig(joinpath(results_figures_folder,"$(name_file_1)_curt_vs_load_growth_with_OPF.png"))
savefig(joinpath(results_figures_folder,"$(name_file_1)_curt_vs_load_growth_with_OPF.pdf"))
savefig(joinpath(results_figures_folder,"$(name_file_1)_curt_vs_load_growth_with_OPF.svg"))



##############################
load_multiplier_2 = 2
result_congested = Dict{String,Any}()
#opf_load_growth_opf(result_og,test_case_opf_2,235,wind_cf,load,load_multiplier_1,growing_load)
opf_load_growth_selected_topology(result_congested,test_case_opf_2,test_case_updated_split_2,switches_couples_2,extremes_ZILs_2,235,1,bs_congested_2,wind_cf,load,load_multiplier_2,growing_load)
opf_load_growth_selected_topology(result_congested,test_case_opf_2,test_case_updated_split_2,switches_couples_2,extremes_ZILs_2,65,2,bs_congested_2,wind_cf,load,load_multiplier_2,growing_load)
opf_load_growth_selected_topology(result_congested,test_case_opf_2,test_case_updated_split_2,switches_couples_2,extremes_ZILs_2,234,4,bs_congested_2,wind_cf,load,load_multiplier_2,growing_load)
opf_load_growth_selected_topology(result_congested,test_case_opf_2,test_case_updated_split_2,switches_couples_2,extremes_ZILs_2,247,4,bs_congested_2,wind_cf,load,load_multiplier_2,growing_load)
#opf_load_growth_opf(result_congested,test_case_opf_2,247,wind_cf,load,load_multiplier_2,growing_load)


result_congested_opf = Dict{String,Any}()
#opf_load_growth_opf(result_og,test_case_opf_2,235,wind_cf,load,load_multiplier_1,growing_load)
opf_load_growth_opf(result_congested_opf,test_case_opf_2,235,wind_cf,load,load_multiplier_2,growing_load)
opf_load_growth_opf(result_congested_opf,test_case_opf_2,65 ,wind_cf,load,load_multiplier_2,growing_load)
opf_load_growth_opf(result_congested_opf,test_case_opf_2,234,wind_cf,load,load_multiplier_2,growing_load)
opf_load_growth_opf(result_congested_opf,test_case_opf_2,247,wind_cf,load,load_multiplier_2,growing_load)




cong_235      = [result_congested["235"]["$gl"]["objective"] for gl in growing_load]
cong_235_curt = [result_congested["235"]["$gl"]["solution"]["gen"]["138"]["pg"] for gl in growing_load]
cong_235_curt_opf = [result_congested_opf["235"]["$gl"]["solution"]["gen"]["138"]["pg"] for gl in growing_load]

cong_65      = [result_congested["65"]["$gl"]["objective"] for gl in growing_load]
cong_65_curt = [result_congested["65"]["$gl"]["solution"]["gen"]["138"]["pg"] for gl in growing_load]
cong_65_ts   = [result_congested["65"]["$gl"]["termination_status"] for gl in growing_load]
cong_65_curt_opf   = [result_congested_opf["65"]["$gl"]["solution"]["gen"]["138"]["pg"] for gl in growing_load]

cong_234      = [result_congested["234"]["$gl"]["objective"] for gl in growing_load]
cong_234_curt = [result_congested["234"]["$gl"]["solution"]["gen"]["138"]["pg"] for gl in growing_load]
cong_234_curt_opf = [result_congested_opf["234"]["$gl"]["solution"]["gen"]["138"]["pg"] for gl in growing_load]

cong_247      = [result_congested["247"]["$gl"]["objective"] for gl in growing_load]
cong_247_curt = [result_congested["247"]["$gl"]["solution"]["gen"]["138"]["pg"] for gl in growing_load]
cong_247_curt_opf = [result_congested_opf["247"]["$gl"]["solution"]["gen"]["138"]["pg"] for gl in growing_load]

Plots.plot(growing_load/test_case_opf_2["load"]["100"]["pd"],cong_235_curt./test_case_opf_2["load"]["100"]["pd"], label = "Low load, low wind  ",grid=:none,xlims = (-0.05,4.05),xticks = 0:0.5:4,yticks = 0:1:8,ylims=(-0.05,8.05),ylabel = "Curtailment / Biggest load in the test case [-]",xlabel = "Load growth bus 69 / Biggest load in the test case [-]")#,xlabelfontsize = 8,ylabelfontsize = 8,xtickfontsize = 6,ytickfontsize = 6)
Plots.plot!(growing_load/test_case_opf_2["load"]["100"]["pd"], cong_65_curt./test_case_opf_2["load"]["100"]["pd"],label = "High load, low wind ")
Plots.plot!(growing_load/test_case_opf_2["load"]["100"]["pd"],cong_234_curt./test_case_opf_2["load"]["100"]["pd"],label = "Low load, high wind ")
Plots.plot!(growing_load/test_case_opf_2["load"]["100"]["pd"],cong_247_curt./test_case_opf_2["load"]["100"]["pd"],label = "High load, high wind")

Plots.plot!(growing_load/test_case_opf_2["load"]["100"]["pd"],cong_235_curt_opf./test_case_opf_2["load"]["100"]["pd"],label = "Low load, low wind (OPF)")
Plots.plot!(growing_load/test_case_opf_2["load"]["100"]["pd"],cong_65_curt_opf./test_case_opf_2["load"]["100"]["pd"],label  = "High load, low wind (OPF)")
Plots.plot!(growing_load/test_case_opf_2["load"]["100"]["pd"],cong_234_curt_opf./test_case_opf_2["load"]["100"]["pd"],label = "Low load, high wind (OPF)")
Plots.plot!(growing_load/test_case_opf_2["load"]["100"]["pd"],cong_247_curt_opf./test_case_opf_2["load"]["100"]["pd"],label = "High load, high wind (OPF)")


savefig(joinpath(results_figures_folder,"$(name_file_2)_curt_vs_load_growth_plus_OPF.png"))
savefig(joinpath(results_figures_folder,"$(name_file_2)_curt_vs_load_growth_plus_OPF.pdf"))
savefig(joinpath(results_figures_folder,"$(name_file_2)_curt_vs_load_growth_plus_OPF.svg"))
