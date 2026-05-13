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
#gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"MIPGap" => 1e-4, "time_limit" => 300)
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"MIPGap" => 5e-4, "time_limit" => 300,"BarQCPConvTol"=>1e-6,"QCPDual" => 1,"ScaleFlag"=>2, "NumericFocus"=>2)
#gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"BarQCPConvTol"=>1e-6,"QCPDual" => 1)#r, "ScaleFlag"=>2, "NumericFocus"=>2) 


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

for (b_id,b) in test_case_1["bus"]
    b["vmax"] = 1.1
    b["vmin"] = 0.9
end
for (b_id,b) in test_case_original_1["bus"]
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
test_case_1["gen"]["138"]["cost"][1] = 9000

test_case_opf_1 = deepcopy(test_case_1)
splitted_bus_ac = [49,46]
name_file_1 = "49_46_congested"
#name_file_2 = "49_46_standard"

load_multiplier = 2

test_case_updated_split_1_result = deepcopy(test_case_1)
test_case_1_result = deepcopy(test_case_1)
test_case_updated_split_1_result,  switches_couples_1_result,  extremes_ZILs_1_result  = _PMTA.AC_busbars_split(test_case_1_result,splitted_bus_ac)

###############
results_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_topologies_selection_results/Results_batch"

function batch_split_every_hour(first_hour,last_hour,size_batch,load_time_series,wind_time_series,optimizer,formulation,results_folder,grid_original,load_multiplier,name_file,splitted_bus_ac)
    for sample in Int64(first_hour):Int64(last_hour/size_batch)
        dict = Dict{String,Any}()
        start_idx = size_batch*(sample-1) + 1
        end_idx = size_batch*(sample-1) + 24
        result = Dict{String,Any}()
        s = Dict("output" => Dict("duals" => false))
        for t in start_idx:end_idx
            println("===================================")
            println("Let's go with timestep $t")
            println("===================================")
            dict["$t"] = Dict{String,Any}()
            test_case = deepcopy(grid_original)       
            test_case_updated_split = deepcopy(grid_original)
            test_case_updated_split,  switches_couples,  extremes_ZILs  = _PMTA.AC_busbars_split(test_case,splitted_bus_ac)    
            test_case_updated_split["gen"]["30"]["pmax"] = grid_original["gen"]["30"]["pmax"] * wind_time_series[t]
            for (l_id,l) in test_case_updated_split["load"]
                if parse(Int64,l_id) < 100 
                    test_case_updated_split["load"][l_id]["pd"] = grid_original["load"][l_id]["pd"] * load_multiplier * load_time_series[t]
                end
            end
            result_switches_lpac  = _PMTA.run_acdcsw_AC_grid(test_case_updated_split,LPACCPowerModel,gurobi)
            dict["$t"] = deepcopy(result_switches_lpac)
        end
        json_opf = JSON.json(dict)        
        open(joinpath(results_folder,"BuS_yearly_$(name_file)_$(start_idx)_$(end_idx).json"),"w") do f 
            write(f, json_opf) 
        end
    end
end

batch_split_every_hour(1,8784,24,time_series.demand,time_series.wind_cf,gurobi,LPACCPowerModel,results_folder,test_case_opf_1,1,name_file_1,splitted_bus_ac)
#batch_split_every_hour(720,8784,24,time_series.demand,time_series.wind_cf,gurobi,LPACCPowerModel,results_folder,test_case_opf_2,2,name_file_2,splitted_bus_ac)




