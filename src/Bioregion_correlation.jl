using CSV

using GLMakie, GeoMakie, GraphMakie

using YAXArrays, DimensionalData

using DataFrames, Statistics

using ADRIA, CoralBlox

include("common.jl")

context_layers = find_latest_file("../canonical-reefs/output/")
context_layers = GDF.read(context_layers)

bioregions = GDF.read("data/GBRMPA_reefal_bioregions.gpkg")
context_bioregion = find_intersections(context_layers, bioregions, :GBRMPA_ID, :DESCRIP, :SHAPE)
context_layers = leftjoin(context_layers, context_bioregion; on=:GBRMPA_ID, matchmissing=:notequal, order=:left)
rename!(context_layers, :area_ID => :bioregion)
context_layers.bioregion .= ifelse.(ismissing.(context_layers.bioregion), "NA", context_layers.bioregion)
context_layers.bioregion = convert.(String, context_layers.bioregion)

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

bioregions = unique(context_layers.bioregion)

lagged_analysis_bior = DataFrame()
for bioregion in bioregions
    bioregion_reefs = context_layers[(context_layers.bioregion .== bioregion), :UNIQUE_ID]
    bioregion_cover = rel_cover[:, rel_cover.sites .∈ [bioregion_reefs]]

    bioregion_lagged_analysis = lagged_region_analysis(bioregion_cover, bioregion, 1:10)
    lagged_analysis_bior = vcat(lagged_analysis_bior, bioregion_lagged_analysis)
end

reefs_p9_l5 = lagged_analysis_bior[(lagged_analysis_bior.lag5 .>= 0.9), :UNIQUE_ID]
f, ga = plot_map(context_layers, :bioregion)
plot_map!(ga, context_layers[(context_layers.UNIQUE_ID .∈ [reefs_p9_l5]),:], color=:orange)
plot_map(context_layers[(context_layers.UNIQUE_ID .∈ [reefs_p9_l5]), :], :closest_port)
