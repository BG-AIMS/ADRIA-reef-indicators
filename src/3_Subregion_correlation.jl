"""
Perfom lagged correlation analysis across reefs within closest_port and bioregion subregions.
There are 15 closes_port subregions and there are ~30 bioregions across the reefs in the GBR.
"""

using GLMakie, GeoMakie, GraphMakie

using YAXArrays, DimensionalData

using DataFrames, Statistics, YAXArrays

using ADRIA, CoralBlox

import GeoDataFrames as GDF

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

# # s_tac = ADRIA.metrics.scenario_total_cover(rs)
tac = ADRIA.metrics.total_absolute_cover(rs)

# reduce all the scenarios down to one series for each reef
tac_sites = Float64.(mapslices(median, tac, dims=[:scenarios]))

# Have to remove the first year as there seems to be an issue with that year's data
tac_sites_reduced = tac_sites[timesteps=2:79]

# calculate the relative site cover from the initial cover across timesteps for each reef
rel_cover = Float64.(mapslices(relative_site_cover, tac_sites_reduced, dims=[:timesteps]))

# Apply analysis to closest_port subregions - 15 subregions
port_subregions = unique(context_layers.closest_port)
lagged_analysis_subregion = subregion_analysis(port_subregions, rel_cover, context_layers, :closest_port, 1:10)
target_reefs_subr = lagged_analysis_subregion[(lagged_analysis_subregion.lag6 .>= 0.8), :RME_UNIQUE_ID]

# Apply analysis to bioregion subregions - 31 subregions
bioregions = unique(context_layers.bioregion)
lagged_analysis_bior = subregion_analysis(bioregions, rel_cover, context_layers, :bioregion, 1:10)
target_reefs_bior = lagged_analysis_bior[(lagged_analysis_bior.lag6 .>= 0.8), :RME_UNIQUE_ID]

# Find common reefs betwen two subregion analyses
target_reefs = target_reefs_bior[(target_reefs_bior .∈ [target_reefs_subr])]

# Add identified reefs to context_layers and write to gpkg with bioregion names
context_layers.target_reefs_subr = context_layers.RME_UNIQUE_ID .∈ [target_reefs_subr]
context_layers.target_reefs_bior = context_layers.RME_UNIQUE_ID .∈ [target_reefs_bior]
context_layers.target_reefs = context_layers.RME_UNIQUE_ID .∈ [target_reefs]
GDF.write("data/context_layers_targetted.gpkg", context_layers; crs=crs=GFT.EPSG(7844))
