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

n_hours = 8784
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
test_case["gen"]["138"]["cost"][1] = 9000

test_case_opf = deepcopy(test_case)
splitted_bus_ac = [49,46]

names = ["49_46_standard","49_46_congested","49_46_standard_data_center","49_46_congested_data_center","69_24_standard","69_24_congested","69_24_standard_data_center","69_24_congested_data_center"]
name_file = names[2]
busbars = "49_46"

test_case_updated_split_result = deepcopy(test_case)
test_case_result = deepcopy(test_case)
test_case_updated_split_result,  switches_couples_result,  extremes_ZILs_result  = _PMTA.AC_busbars_split(test_case_result,splitted_bus_ac)

###
# Upload results
results_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results/Validation/Full_results"
###############
hours = [[24,2160],[2160,4320],[4320,6480],[6480,8784]]

#=
for n in names
    results_dict_choose_topology = Dict{String,Any}()
    for h in 1:length(hours)
        start_hour = hours[h][1]
        end_hour = hours[h][2]
        result_dict = JSON.parsefile(joinpath(results_folder,"Choose_topology_$(n)_$(start_hour)_$(end_hour)_4_configurations.json"))

        if start_hour == 24        
            real_start = 1
            for t in real_start:end_hour
                results_dict_choose_topology["$t"] = result_dict[string(t)]
            end
        else
            for t in start_hour:end_hour
                results_dict_choose_topology["$t"] = result_dict[string(t)]
            end
        end
    end
    json_dict = JSON.json(results_dict_choose_topology)        
    open(joinpath(results_folder,"Choose_topology_$(n)_4_configurations.json"),"w") do f 
        write(f, json_dict) 
    end
end
=#
results_dict_opf = JSON.parsefile(joinpath(results_folder,"OPF_$(name_file)_24_8784.json"))
results_dict_bus = JSON.parsefile(joinpath(results_folder,"Bus_$(name_file)_24_8784.json"))
results_dict_choose_topology = JSON.parsefile(joinpath(results_folder,"Choose_topology_$(name_file)_4_configurations.json"))
results_dict_1_conf = JSON.parsefile(joinpath(dirname(results_folder),"1_configurations_$(name_file)_24_8784.json"))
results_dict_4_conf = JSON.parsefile(joinpath(results_folder,"4_configurations_$(name_file)_24_8784.json"))

confs = keys(results_dict_4_conf["1"])
#results_dict_4_conf["1"][conf]["objective"] for conf in confs


tot_opf = sum(results_dict_opf[string(t)]["objective"] for t in 1:8784)
objs_opf = [results_dict_opf[string(t)]["objective"] for t in 1:8784]
term_status_opf = [results_dict_opf[string(t)]["termination_status"] for t in 1:8784]
countmap(term_status_opf)
sum(objs_opf)

obj_choose_topology = sum(results_dict_choose_topology[string(t)]["objective"] for t in 1:8784)
term_status_choose_topology = [results_dict_choose_topology[string(t)]["termination_status"] for t in 1:8784]
countmap(term_status_choose_topology)

objs_1_conf = [results_dict_1_conf[string(t)]["114"]["objective"] for t in 1:8784]
sum(objs_1_conf)


objs_1_conf_total = []
count_bus = 0
count_opf = 0
for i in 1:n_hours
    diff_opf_1_conf = objs_opf[i] .- objs_1_conf[i]
    if diff_opf_1_conf > 0
        push!(objs_1_conf_total,objs_1_conf[i])
        count_bus += 1
    else
        push!(objs_1_conf_total,objs_opf[i])
        count_opf += 1
    end
end
tot_objs_1_conf = sum(objs_1_conf_total)
(1 - tot_objs_1_conf/sum(objs_opf))*100


term_status_1_conf = [results_dict_1_conf[string(t)]["114"]["termination_status"] for t in 1:8784]


objs_bus = [results_dict_bus[string(t)]["objective"] for t in 1:8784]
term_status_bus = [results_dict_bus[string(t)]["termination_status"] for t in 1:8784]
tot_bus = sum(objs_bus)

diff_opf_conf = objs_opf .- objs_bus
scatter(diff_opf_1_conf)
hline!([0.0])

(1 - tot_bus/tot_opf)*100

wind_plot = scatter(cap_factor_69,grid =:none,label = :none,xlabel="Hour", ylabel="Wind capacity factor")
load_plot = scatter(cap_factor_load,grid =:none, label = :none,xlabel="Hour", ylabel="Load capacity factor")

results_figures_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results/Figures"
savefig(wind_plot,joinpath(results_figures_folder,"Wind_standard.png"))
savefig(load_plot,joinpath(results_figures_folder,"Load_standard.png"))
savefig(wind_plot,joinpath(results_figures_folder,"Wind_standard.pdf"))
savefig(load_plot,joinpath(results_figures_folder,"Load_standard.pdf"))


keys_conf = eachindex(results_dict_4_conf["1"])
objs_4_conf_total = []
count_bus_4_confs = 0
count_opf_4_confs = 0

for t in 1:8784
    objs = []
    for conf in keys_conf
        push!(objs, results_dict_4_conf[string(t)][conf]["objective"])
    end
    obj_4_bus = minimum(objs)
    if obj_4_bus < objs_opf[t]
        push!(objs_4_conf_total,obj_4_bus)
        count_bus_4_confs += 1
    else
        push!(objs_4_conf_total,objs_opf[t])
        count_opf_4_confs += 1
    end
end
(1 - sum(objs_4_conf_total)/sum(objs_opf))*100



diff_opf_4_conf = (objs_opf .- obj_4_confs)./objs_opf*100
scatter(diff_opf_4_conf,xlabel="Hour", ylabel="Relative difference with OPF (%)", grid = :none)
tot_4_confs = sum(obj_4_confs)



sum(objs_opf) - tot_4_confs

names_1_4 = [names[1]]


function save_simulation_comparison_results(results_folder, names)
    results_dict = Dict{String,Any}()
    for n in names
        results_dict[n] = Dict{String,Any}()
        results_dict_opf = JSON.parsefile(joinpath(results_folder,"OPF_$(n)_24_8784.json"))
        if n == "49_46_congested"
            results_dict_bus = JSON.parsefile(joinpath(dirname(results_folder),"BuS_$(n)_24_8784.json"))
        else
            results_dict_bus = JSON.parsefile(joinpath(results_folder,"Bus_$(n)_24_8784.json"))
        end
        results_dict_choose_topology = JSON.parsefile(joinpath(results_folder,"Choose_topology_$(n)_4_configurations.json"))
        results_dict_1_conf = JSON.parsefile(joinpath(dirname(results_folder),"1_configurations_$(n)_24_8784_no_137.json"))
        results_dict_4_conf = JSON.parsefile(joinpath(results_folder,"4_configurations_$(n)_24_8784_no_137.json"))
        keys_conf = eachindex(results_dict_4_conf["1"])

        objs_opf = [results_dict_opf[string(t)]["objective"] for t in 1:8784]
        tot_opf = sum(results_dict_opf[string(t)]["objective"] for t in 1:8784)
        keys_1_conf = keys(results_dict_1_conf["1"])

        obj_1_confs = []
        for t in 1:8784
            objs = []
            for conf in keys_1_conf
                push!(objs, results_dict_1_conf[string(t)][conf]["objective"])
            end
            push!(objs,objs_opf[t])
            push!(obj_1_confs, minimum(objs))
        end
        tot_1_conf = sum(obj_1_confs)

        tot_bus = 0
        for t in 1:8784
            if results_dict_bus[string(t)]["termination_status"] == "OPTIMAL"
                println("Processing hour $t...")
                tot_bus += results_dict_bus[string(t)]["objective"]
            else 
                tot_bus += results_dict_opf[string(t)]["objective"]
            end
        end
        obj_4_confs = []
        for t in 1:8784
            objs = []
            for conf in keys_conf
                push!(objs, results_dict_4_conf[string(t)][conf]["objective"])
            end
            push!(objs,objs_opf[t])
            push!(obj_4_confs, minimum(objs))
        end
        diff_opf_4_conf = objs_opf .- obj_4_confs
        #scatter(diff_opf_4_conf)
        tot_4_confs = sum(obj_4_confs)
        results_dict[n] = Dict(
            "tot_opf"     => deepcopy(tot_opf),
            "tot_1_conf"  => deepcopy(tot_1_conf),
            "tot_bus"     => deepcopy(tot_bus),
            "tot_4_confs" => deepcopy(tot_4_confs)
        )
        json_dict = JSON.json(results_dict)
        open(joinpath(results_folder,"Comparison_results_$(n).json"),"w") do f 
            write(f, json_dict) 
        end
    end
    return results_dict
end

results_dict = save_simulation_comparison_results(results_folder, names_1_4)

results_dict["49_46_standard"]
results_dict["49_46_congested"]

(1 - results_dict["49_46_standard"]["tot_bus"]/results_dict["49_46_standard"]["tot_opf"])*100
(1 - results_dict["49_46_standard"]["tot_1_conf"]/results_dict["49_46_standard"]["tot_opf"])*100
(1 - results_dict["49_46_standard"]["tot_4_confs"]/results_dict["49_46_standard"]["tot_opf"])*100


(1 - results_dict["49_46_congested"]["tot_bus"]    /results_dict["49_46_congested"]["tot_opf"])*100
(1 - results_dict["49_46_congested"]["tot_1_conf"] /results_dict["49_46_congested"]["tot_opf"])*100
(1 - results_dict["49_46_congested"]["tot_4_confs"]/results_dict["49_46_congested"]["tot_opf"])*100


results_dict_bus = JSON.parsefile(joinpath(dirname(results_folder),"BuS_$(name_file)_24_8784.json"))

keys_conf_4 = eachindex(results_dict_4_conf["1"])
objs_opf_cong = [results_dict_opf[string(t)]["objective"] for t in 1:8784]

obj_4_confs = []
for t in 1:8784
    objs = []
    for conf in keys_conf_4
        push!(objs, results_dict_4_conf[string(t)][conf]["objective"])
    end
    push!(objs,objs_opf_cong[t])
    push!(obj_4_confs, minimum(objs))
end
diff_opf_4_conf = objs_opf_cong .- obj_4_confs

diff_opf_4_conf = (objs_opf_cong .- obj_4_confs)./objs_opf_cong*100
cong_49_46 = scatter(diff_opf_4_conf,xlabel="Hour", ylabel="Relative difference with OPF (%)", grid = :none,label = :none)

savefig(cong_49_46,joinpath(results_figures_folder,"4_confs_vs_opf_$(name_file).pdf"))
savefig(cong_49_46,joinpath(results_figures_folder,"4_confs_vs_opf_$(name_file).png"))

objs = Dict{String,Any}()
obj_1_confs = []
obj_topology = []
objs["count"] = Dict{String,Any}()
objs["count"]["OPF"] = 0
for key in keys_conf
    objs["count"][key] = 0
end
for t in 1:8784
    objs["$t"] = []
    for key in keys_conf
        push!(objs["$t"], [results_dict_4_conf[string(t)][key]["objective"],key])
    end
    push!(objs["$t"],[objs_opf[t],"OPF"])
    push!(obj_1_confs, minimum(objs["$t"]))
    push!(obj_topology, minimum(objs["$t"])[2])
    objs["count"][minimum(objs["$t"])[2]] += 1
end

timesteps = Dict{String,Any}() 
for t in 1:5
    timesteps["$t"] = Dict{String,Any}()
    timesteps["$t"]["wind"] = []
    timesteps["$t"]["load"] = []
    timesteps["$t"]["t"] = []
end


for t in 1:8784
    if obj_topology[t] == "227" 
        push!(timesteps["1"]["wind"], cap_factor_69[t])
        push!(timesteps["1"]["load"], cap_factor_load[t])
        push!(timesteps["1"]["t"],t)
    elseif obj_topology[t] == "20" 
        push!(timesteps["2"]["wind"], cap_factor_69[t])
        push!(timesteps["2"]["load"], cap_factor_load[t])
        push!(timesteps["2"]["t"],t)
    elseif obj_topology[t] == "168" 
        push!(timesteps["3"]["wind"], cap_factor_69[t])
        push!(timesteps["3"]["load"], cap_factor_load[t])
        push!(timesteps["3"]["t"],t)
    elseif obj_topology[t] == "3"
        push!(timesteps["4"]["wind"], cap_factor_69[t])
        push!(timesteps["4"]["load"], cap_factor_load[t])
        push!(timesteps["4"]["t"],t)
    elseif obj_topology[t] == "OPF"
        push!(timesteps["5"]["wind"], cap_factor_69[t])
        push!(timesteps["5"]["load"], cap_factor_load[t])
        push!(timesteps["5"]["t"],t)
    end
end
#=
for t in 1:8784
    if obj_topology[t] == "4" 
        push!(timesteps["1"]["wind"], cap_factor_69[t])
        push!(timesteps["1"]["load"], cap_factor_load[t])
        push!(timesteps["1"]["t"],t)
    elseif obj_topology[t] == "1" 
        push!(timesteps["2"]["wind"], cap_factor_69[t])
        push!(timesteps["2"]["load"], cap_factor_load[t])
        push!(timesteps["2"]["t"],t)
    elseif obj_topology[t] == "114" 
        push!(timesteps["3"]["wind"], cap_factor_69[t])
        push!(timesteps["3"]["load"], cap_factor_load[t])
        push!(timesteps["3"]["t"],t)
    elseif obj_topology[t] == "137"
        push!(timesteps["4"]["wind"], cap_factor_69[t])
        push!(timesteps["4"]["load"], cap_factor_load[t])
        push!(timesteps["4"]["t"],t)
    elseif obj_topology[t] == "OPF"
        push!(timesteps["5"]["wind"], cap_factor_69[t])
        push!(timesteps["5"]["load"], cap_factor_load[t])
        push!(timesteps["5"]["t"],t)
    end
end
=#
tl_4 = timesteps["4"]["load"]
tl_5 = timesteps["5"]["load"]

tw_4 = timesteps["4"]["wind"]
tw_5 = timesteps["5"]["wind"]

tt_4 = timesteps["4"]["t"]
tt_5 = timesteps["5"]["t"]


load_opf = vcat(tl_4,tl_5)
wind_opf = vcat(tw_4,tw_5)
t_opf    = vcat(tt_4,tt_5)

p_distribution = scatter([], [], label="Topology 1", color=:blue)
scatter!(p_distribution,timesteps["4"]["load"], timesteps["4"]["wind"] , label="Topology 2", color=:green)
scatter!(p_distribution,timesteps["3"]["load"], timesteps["3"]["wind"] , label="Topology 3", color=:orange)
scatter!(p_distribution,timesteps["1"]["load"], timesteps["1"]["wind"], label="Topology 4",color=:red)
scatter!(p_distribution,timesteps["5"]["load"], timesteps["5"]["wind"],label="Original",color=:gray,xlabel= "Demand (p.u.)",ylabel="Wind capacity factor [-]",legend=:topright,grid =:none,xticks = 0:0.1:1,xlims = (0.25,1))
#scatter!(p_distribution,load_opf,wind_opf,label="Original",color=:gray,xlabel= "Demand (p.u.)",ylabel="Wind capacity factor [-]",legend=:topright,grid =:none,xticks = 0:0.1:1,xlims = (0.25,1))

savefig(p_distribution,joinpath(results_figures_folder,"4_confs_vs_opf_$(name_file)_distribution.pdf"))
savefig(p_distribution,joinpath(results_figures_folder,"4_confs_vs_opf_$(name_file)_distribution.png"))
