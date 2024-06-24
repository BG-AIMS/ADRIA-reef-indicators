using CSV

using GLMakie, GeoMakie, GraphMakie

using YAXArrays, DimensionalData

using DataFrames, Statistics

using ADRIA, CoralBlox

include("common.jl")

context_layers = find_latest_file("../canonical-reefs/output/")
context_layers = GDF.read(context_layers)

rs = ADRIA.load_results("outputs/ADRIA-out/ReefMod Engine__RCPs_45__2024-06-19_10_58_55_647")

fig_opts = Dict(:size => (1600, 800))



# s_tac = ADRIA.metrics.scenario_total_cover(rs)
tac = ADRIA.metrics.total_absolute_cover(rs)

# reduce all the scenarios down to one series for each reef
tac_sites = mapslices_toFloat64(median, tac, :scenarios)

# Have to remove the first year as there seems to be an issue with that year's data
tac_sites_reduced = tac_sites[timesteps=2:79]

# calculate the relative site cover from the initial cover across timesteps for each reef
rel_cover = mapslices_toFloat64(relative_site_cover, tac_sites_reduced, :timesteps)

port_subregions = unique(context_layers.closest_port)

lagged_analysis_sub = DataFrame()
for subregion in port_subregions
    subregion_reefs = context_layers[(context_layers.closest_port .== subregion), :UNIQUE_ID]
    subregion_cover = rel_cover[:, rel_cover.sites .∈ [subregion_reefs]]

    subregion_lagged_analysis = lagged_region_analysis(subregion_cover, subregion, 1:10)
    lagged_analysis_sub = vcat(lagged_analysis_sub, subregion_lagged_analysis)
end

reefs_p8_l5 = lagged_analysis_sub[(lagged_analysis_sub.lag5 .>= 0.85), :UNIQUE_ID]
f, ga = plot_map(context_layers, :closest_port)
plot_map!(ga, context_layers[(context_layers.UNIQUE_ID .∈ [reefs_p8]),:], color=:black)
series!(Matrix(lagged_analysis_sub[:,3:12]))
