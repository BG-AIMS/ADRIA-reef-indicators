using CSV

using GLMakie, GeoMakie, GraphMakie

using YAXArrays, DimensionalData

using DataFrames, Statistics

using ADRIA, CoralBlox

include("common.jl")

context_layers = find_latest_file("../canonical-reefs/output/")
context_layers = GDF.read(context_layers)

rs = ADRIA.load_results("outputs/ADRIA-out/ReefMod Engine__RCPs_45__2024-06-19_10_58_55_647")
tac = ADRIA.metrics.total_absolute_cover(rs)

tac_sites = mapslices(median, tac, dims=[:scenarios])


n_clusters = 5
clusters = ADRIA.analysis.cluster_series(tac_sites, n_clusters)
tac_clusters = ADRIA.metrics.per_loc(median, tac)
tsc_map_fig = ADRIA.viz.map(rs, tac_clusters, clusters)
