using CSV

using GLMakie, GeoMakie, GraphMakie

using YAXArrays, DimensionalData

using DataFrames, Statistics

using ADRIA, CoralBlox

include("common.jl")

context_layers = find_latest_file("../canonical-reefs/output/")
context_layers = GDF.read(context_layers)

rs = ADRIA.load_results("outputs/ADRIA-out/ReefMod Engine__RCPs_45__2024-06-19_10_58_55_647")
#s_tac = ADRIA.metrics.scenario_total_cover(rs)
tac = ADRIA.metrics.total_absolute_cover(rs)

tac_sites_init = mapslices(median, tac, dims=[:scenarios])
data = convert.(Float64, tac_sites_init.data)
tac_sites::YAXArray{Float64, 2} = YAXArray(dims(tac_sites_init), data)

tac_sites_reduced = tac_sites[timesteps=2:79]
rel_cover_init = mapslices(relative_site_cover, tac_sites_reduced, dims=[:timesteps])
data = convert.(Float64, rel_cover_init.data)
rel_cover::YAXArray{Float64, 2} = YAXArray(dims(rel_cover_init), data)

n_clusters = 5
clusters = ADRIA.analysis.cluster_scenarios(rel_cover, n_clusters)

tsc_fig = ADRIA.viz.clustered_scenarios(
    rel_cover, clusters
)

tac_clusters = ADRIA.metrics.per_loc(median, tac)
tsc_map_fig = ADRIA.viz.map(rs, tac_clusters, clusters)
