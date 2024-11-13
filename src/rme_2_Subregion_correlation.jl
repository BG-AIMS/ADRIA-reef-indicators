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

context_layers = find_latest_file("../../canonical-reefs/output/")
context_layers = GDF.read(context_layers)
canonical_reefs = GDF.read(find_latest_file("../../canonical-reefs/output/"))
# Add Bioregions to context layers
bioregions = GDF.read("../data/GBRMPA_reefal_bioregions.gpkg")
context_bioregion = find_intersections(context_layers, bioregions, :GBRMPA_ID, :DESCRIP, :SHAPE)
context_layers = leftjoin(context_layers, context_bioregion; on=:GBRMPA_ID, matchmissing=:notequal, order=:left)
rename!(context_layers, :area_ID => :bioregion)
context_layers.bioregion .= ifelse.(ismissing.(context_layers.bioregion), "NA", context_layers.bioregion)
context_layers.bioregion = convert.(String, context_layers.bioregion)

rs = open_dataset("../outputs/RME_result_stores/RME_SSP245_200reps/median_GCM_cover_and_evenness_2024_11_11.nc")

# Align the order of context_layers with the order of sites in RME
context_layers = context_layers[indexin(rs.scaled_taxa_evenness.sites, context_layers.RME_UNIQUE_ID), :]

cover_ts = rs.total_relative_cover
cover_ts = cover_ts[1:50, :]
taxa_evenness = rs.scaled_taxa_evenness
taxa_evenness = taxa_evenness[1:50, :]

# Apply analysis to management regions - 4 regions
# management_regions = unique(context_layers.management_area)
# lagged_analysis_mgmt = subregion_analysis(management_regions, cover_ts, context_layers, :management_area, 1:10)
# target_reefs_mgmt = lagged_analysis_mgmt[(lagged_analysis_mgmt[:,"lag4"] .>= 0.7), :RME_UNIQUE_ID]

# lagged_analysis_mgmt_evenness = subregion_analysis(management_regions, taxa_evenness, context_layers, :management_area, 1:10)
# target_reefs_mgmt_evenness = lagged_analysis_mgmt_evenness[(lagged_analysis_mgmt_evenness[:,"lag5"] .>= 0.7), :RME_UNIQUE_ID]

# # Apply analysis to closest_port subregions - 15 subregions
# port_subregions = unique(context_layers.closest_port)
# lagged_analysis_subregion = subregion_analysis(port_subregions, cover_ts, context_layers, :closest_port, 1:10)
# target_reefs_subr = lagged_analysis_subregion[(lagged_analysis_subregion[:,"lag4"] .>= 0.7), :RME_UNIQUE_ID]

# lagged_analysis_subregion_evenness = subregion_analysis(port_subregions, taxa_evenness, context_layers, :closest_port, 1:10)
# target_reefs_subr_evenness = lagged_analysis_subregion_evenness[(lagged_analysis_subregion_evenness[:,"lag5"] .>= 0.7), :RME_UNIQUE_ID]

# Apply analysis to bioregion subregions - 31 subregions
bioregions = unique(context_layers.bioregion)

lagged_analysis_bior = subregion_analysis(bioregions, cover_ts, context_layers, :bioregion, 1:10)
context_layers = leftjoin(context_layers, lagged_analysis_bior[:, ["lag4", "RME_UNIQUE_ID"]], on=:RME_UNIQUE_ID, order=:left)
rename!(context_layers, :lag4 => :lag4_bior)
context_layers.bior_cover_qc_flag = ifelse.(isnan.(context_layers.lag4_bior), 1, 0)
context_layers.lag4_bior = ifelse.(isnan.(context_layers.lag4_bior), 0.0, context_layers.lag4_bior)

lagged_analysis_bior_evenness = subregion_analysis(bioregions, taxa_evenness, context_layers, :bioregion, 1:10)
context_layers = leftjoin(context_layers, lagged_analysis_bior_evenness[:, ["lag4", "RME_UNIQUE_ID"]], on=:RME_UNIQUE_ID, order=:left)
rename!(context_layers, :lag4 => :lag4_bior_evenness)
context_layers.bior_evenness_qc_flag = ifelse.(isnan.(context_layers.lag4_bior_evenness), 1, 0)
context_layers.lag4_bior_evenness = ifelse.(isnan.(context_layers.lag4_bior_evenness), 0.0, context_layers.lag4_bior_evenness)

# Find common reefs betwen two subregion analyses
# target_reefs = target_reefs_bior[(target_reefs_bior .∈ [target_reefs_subr])]
# target_reefs_evenness = target_reefs_bior_evenness[(target_reefs_bior_evenness .∈ [target_reefs_subr_evenness])]

# Add identified reefs to context_layers and write to gpkg with bioregion names
# context_layers.target_reefs_mgmt = context_layers.RME_UNIQUE_ID .∈ [target_reefs_mgmt]
# context_layers.target_reefs_subr = context_layers.RME_UNIQUE_ID .∈ [target_reefs_subr]
context_layers.target_reefs_bior = context_layers.lag4_bior .> 0.7
# context_layers.target_reefs = context_layers.RME_UNIQUE_ID .∈ [target_reefs]

# context_layers.target_reefs_mgmt_evenness = context_layers.RME_UNIQUE_ID .∈ [target_reefs_mgmt_evenness]
# context_layers.target_reefs_subr_evenness = context_layers.RME_UNIQUE_ID .∈ [target_reefs_subr_evenness]
context_layers.target_reefs_bior_evenness = context_layers.lag4_bior_evenness .> 0.7
# context_layers.target_reefs_evenness = context_layers.RME_UNIQUE_ID .∈ [target_reefs_evenness]
GDF.write("../data/context_layers_targetted_rme.gpkg", context_layers; crs=GFT.EPSG(7844), overwrite=true)
