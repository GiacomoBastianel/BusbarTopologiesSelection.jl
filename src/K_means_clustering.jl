#=
using PowerModelsTopologicalActionsII; const _PMTA = PowerModelsTopologicalActionsII
using PowerModels; const _PM = PowerModels
using PowerModelsACDC; const _PMACDC = PowerModelsACDC
using JuMP, Ipopt, JSON, HiGHS
using Gurobi
using PowerPlots
using CSV
using DataFrames
using StatsBase
using Plots
=#
using DataFrames
using CSV
using Statistics
using Clustering
using LinearAlgebra


wind_series = CSV.read(joinpath(test_case_folder,"RTS_GMLC_data","RTS_Data","timeseries_data_files","WIND","DAY_AHEAD_wind.csv"),DataFrame)
wind_69 = wind_series[:,7]
maximum(wind_69)
cap_factor_69 = wind_69 ./ maximum(wind_69)

load_series = CSV.read(joinpath(test_case_folder,"RTS_GMLC_data","RTS_Data","timeseries_data_files","Load","DAY_AHEAD_regional_Load.csv"),DataFrame)
load_series_2 = load_series[:,6]
maximum(load_series_2)
cap_factor_load = load_series_2 ./ maximum(load_series_2)
maximum(cap_factor_load)
Plots.scatter(cap_factor_load)

time_series = DataFrame(
    demand = cap_factor_load,      # MW
    wind_cf = cap_factor_69     # ∈ [0,1]
)

k = 365  # number of representative time steps
X = hcat(
    time_series.demand,
    time_series.wind_cf
)
X_correct = transpose(X)

R = kmeans(X_correct, k; maxiter=300, display=:none)

N = nrow(time_series)
weights = [count(==(i), R.assignments) / N for i in 1:k]
demand_rep = R.centers[1, :]  #.* std(time_series.demand) .+ mean(time_series.demand)
wind_cf_rep = R.centers[2, :] #.* std(time_series.wind_cf) .+ mean(time_series.wind_cf)
rep = DataFrame(
    demand = demand_rep,
    wind_cf = wind_cf_rep,
    weight = weights
)
Plots.scatter(time_series.demand, time_series.wind_cf,label = "Original",xlabel = "Demand [-]",ylabel = "Wind capacity factor [-]",grid=:none)
Plots.scatter!(rep.demand, rep.wind_cf,color = :red,label = "K-Means clustering")
#dx = 0.01
#dy = 0.01
#for w in 1:length(rep.weight)
#    annotate!(
#        rep.demand[w] .+ dx,
#        rep.wind_cf[w] .+ dy,
#        text(rep.weight[w]), 8, :red)
#end
CSV.write(joinpath(test_case_folder,"Daily_K_Means_clustering.csv"),rep)
