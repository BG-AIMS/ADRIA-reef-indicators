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

for GCM in gcms
    rs = open_dataset("../outputs/RME_result_stores/RME_SSP245_300reps/RME_SSP245_300reps_$(GCM)/cover_and_evenness_2024_12_18.nc"; driver=:netcdf)


    # Align the order of context_layers with the order of sites in RME
    context_layers = context_layers[indexin(rs.scaled_taxa_evenness.sites, context_layers.RME_UNIQUE_ID), :]

    cover_ts = rs.total_relative_cover_median
    cover_ts = cover_ts[1:50, :]
    taxa_evenness = rs.scaled_taxa_evenness
    taxa_evenness = taxa_evenness[1:50, :]

    context_layers[:, "$(GCM)_lag0_bior"] .= 0.0
    context_layers[:, "$(GCM)_lag3_bior"] .= 0.0
    context_layers[:, "$(GCM)_lag4_bior"] .= 0.0
    context_layers[:, "$(GCM)_lag5_bior"] .= 0.0

    # Apply analysis to bioregion subregions - 31 subregions
    context_layers[:, "$(GCM)_lag0_bior"] = cluster_correlation(context_layers.bioregion, cover_ts, 0, CF)
    context_layers[:, "$(GCM)_lag3_bior"] = cluster_correlation(context_layers.bioregion, cover_ts, 3, CF)
    context_layers[:, "$(GCM)_lag4_bior"] = cluster_correlation(context_layers.bioregion, cover_ts, 4, CF)
    context_layers[:, "$(GCM)_lag5_bior"] = cluster_correlation(context_layers.bioregion, cover_ts, 5, CF)
    # lagged_analysis_bior = subregion_analysis(bioregions, cover_ts, context_layers, :bioregion, 1:10)
    # context_layers = leftjoin(context_layers, lagged_analysis_bior[:, ["lag4", "RME_UNIQUE_ID"]], on=:RME_UNIQUE_ID, order=:left)
    # rename!(context_layers, "lag4" => "$(GCM)_lag4_bior")
    context_layers[:, "$(GCM)_bior_cover_qc_flag"] = ifelse.(isnan.(context_layers[:, "$(GCM)_lag4_bior"]), 1, 0)
    context_layers[:, "$(GCM)_lag4_bior"] = Float64.(
        ifelse.(
            isnan.(context_layers[:, "$(GCM)_lag4_bior"]),
            0.0,
            context_layers[:, "$(GCM)_lag4_bior"]
        )
    )

    context_layers[:, "$(GCM)_lag0_bior_evenness"] .= 0.0
    context_layers[:, "$(GCM)_lag3_bior_evenness"] .= 0.0
    context_layers[:, "$(GCM)_lag4_bior_evenness"] .= 0.0
    context_layers[:, "$(GCM)_lag5_bior_evenness"] .= 0.0

    # Apply analysis to bioregion subregions - 31 subregions
    context_layers[:, "$(GCM)_lag0_bior_evenness"]= cluster_correlation(context_layers.bioregion, taxa_evenness, 0, CF)
    context_layers[:, "$(GCM)_lag3_bior_evenness"]= cluster_correlation(context_layers.bioregion, taxa_evenness, 3, CF)
    context_layers[:, "$(GCM)_lag4_bior_evenness"]= cluster_correlation(context_layers.bioregion, taxa_evenness, 4, CF)
    context_layers[:, "$(GCM)_lag5_bior_evenness"] = cluster_correlation(context_layers.bioregion, taxa_evenness, 5, CF)

    #gged_analysis_bior_evenness = subregion_analysis(bioregions, taxa_evenness, context_layers, :bioregion, 1:10)
    # context_layers = leftjoin(context_layers, lagged_analysis_bior_evenness[:, ["lag4", "RME_UNIQUE_ID"]], on=:RME_UNIQUE_ID, order=:left)
    # rename!(context_layers, "lag4" => "$(GCM)_lag4_bior_evenness")
    context_layers[:, "$(GCM)_bior_evenness_qc_flag"] = ifelse.(isnan.(context_layers[:, "$(GCM)_lag4_bior_evenness"]), 1, 0)
    context_layers[:, "$(GCM)_lag4_bior_evenness"] = Float64.(
        ifelse.(
            isnan.(context_layers[:, "$(GCM)_lag4_bior_evenness"]),
            0.0,
            context_layers[:, "$(GCM)_lag4_bior_evenness"]
        )
    )

    context_layers[:, "$(GCM)_target_reefs_bior"] = (
        #(context_layers[:, "$(GCM)_lag4_bior"] .<= quantile(context_layers[:, "$(GCM)_lag4_bior"], 0.20)) .&
        ((context_layers[:, "$(GCM)_lag4_bior"] .- context_layers[:, "$(GCM)_lag0_bior"]) .< -0.2)
    )
    context_layers[:, "$(GCM)_target_reefs_bior_cat"] = ifelse.(context_layers[:, "$(GCM)_target_reefs_bior"], "bellwether", "non-bellwether")

    context_layers[:, "$(GCM)_target_reefs_bior_evenness"] = (
        (context_layers[:, "$(GCM)_lag4_bior_evenness"] .<= quantile(context_layers[:, "$(GCM)_lag4_bior_evenness"], 0.20)) .&
        ((context_layers[:, "$(GCM)_lag4_bior_evenness"] .- context_layers[:, "$(GCM)_lag0_bior_evenness"]) .> 0)
    )
    context_layers[:, "$(GCM)_target_reefs_bior_evenness_cat"] = ifelse.(context_layers[:, "$(GCM)_target_reefs_bior_evenness"], "bellwether", "non-bellwether")

    context_layers[!, "$(GCM)_lag4_bior"] = convert.(Float64, context_layers[:, "$(GCM)_lag4_bior"])
    context_layers[!, "$(GCM)_lag4_bior_evenness"] = convert.(Float64, context_layers[:, "$(GCM)_lag4_bior_evenness"])
end

GDF.write("../data/context_layers_targetted_rme_CF.gpkg", context_layers; crs=GFT.EPSG(7844), overwrite=true)
