"""
Access the domain connectivity and calculate the strength and number of outgoing/incoming
connections for each reef. Run after running 3_.jl to include the bioregions and identified
reef layers.
"""

using DataFrames, Statistics, YAXArrays

using GLMakie, GeoMakie, GraphMakie

using ADRIA

import GeoDataFrames as GDF

import GeoFormatTypes as GFT

include("../common.jl")
include("custom_RMEDomain_dataloading.jl")
gcms = ["EC_Earth3_Veg", "UKESM1_0_LL", "GFDL-ESM4", "CNRM_ESM2_1"]

# Load context layers with bioregions and target reefs
context_layers = GDF.read("../data/context_layers_targetted_rme_CF.gpkg")
context_layers.RME_GBRMPA_ID = ifelse.((context_layers.RME_GBRMPA_ID .== "20198"), "20-198", context_layers.RME_GBRMPA_ID)

# GBR wide domain
gbr_domain_path = "c:/Users/bgrier/Documents/Projects/ADRIA_Domains/CUSTOM_rme_v1033_w_1028conn/"
gbr_dom = ADRIA.load_domain(RMEDomain, gbr_domain_path, "45")

# 1. Attach connectivity data to context_layers
connectivity_matrix = gbr_dom.conn

# Calculate connectivity metrics for context analysis
reefs = collect(getAxis("Source", connectivity_matrix).val)

col_types = [
    String,
    Float64,
    Int64,
    Float64,
    Float64,
    Int64,
    Float64,
    Float64,
    Int64,
    Float64,
    Float64
]

source_to_sink = DataFrame(
    [T[] for T in col_types],
    [
    :RME_GBRMPA_ID,
    :income_strength,
    :income_count,
    :income_comb,
    :out_strength,
    :out_count,
    :out_comb,
    :total_strength,
    :total_count,
    :total_comb,
    :so_to_si
    ]
)

for reef in eachindex(reefs)
    reef_id = reefs[reef]
    outgoing = connectivity_matrix[connectivity_matrix.Source .== reef_id, :]
    incoming = connectivity_matrix[:, connectivity_matrix.Sink .== reef_id]

    income_strength = sum(incoming)
    income_count = count(incoming .> 0)
    income_comb = income_strength / income_count

    out_strength = sum(outgoing)
    out_count = count(outgoing .> 0)
    out_comb = out_strength / out_count

    total_strength = income_strength + out_strength
    total_count = income_count + out_count
    total_comb = total_strength / total_count

    so_to_si = out_comb / income_comb

    push!(
        source_to_sink,
        [
        reef_id,
        income_strength,
        income_count,
        income_comb,
        out_strength,
        out_count,
        out_comb,
        total_strength,
        total_count,
        total_comb,
        so_to_si
        ]
    )
end
context_layers = leftjoin(context_layers, source_to_sink; on=:RME_GBRMPA_ID, order=:left)

# Attaching GCM-dependent context layers
for GCM in gcms

    rs = open_dataset("../outputs/RME_result_stores/RME_SSP245_300reps/RME_SSP245_300reps_$(GCM)/cover_and_evenness_2024_12_18.nc"; driver=:netcdf)
    # 2. Attach initial coral cover data to context_layers
    initial_cc = rs.total_relative_cover_median[1, :]
    context_layers[:, "$(GCM)_initial_coral_cover"] = initial_cc.data

    initial_evenness = rs.scaled_taxa_evenness[1, :]
    context_layers[:, "$(GCM)_initial_evenness"] = initial_evenness

    # 3. Attach DHW data to context_layers
    # mean DHW data for sites
    #dhw_time_scens = load_DHW("../../RME/rme_ml_2024_06_13/data_files/", "245", GCM)
    #dhw_time = Float64.(mapslices(median, dhw_time_scens, dims=[:scenarios]))
    #dhw_locs = Float64.(mapslices(mean, dhw_time[1:50, :], dims=[:timesteps])) # Select only the years in relevant time
    dhw_time = rs.dhw
    dhw_locs = Float64.(mapslices(mean, dhw_time[1:50, :], dims=[:timesteps]))

    year_when_8DHW = zeros(Int64, length(dhw_time.sites))
    for loc in 1:size(dhw_time, 2)
        year = findfirst(x -> x > 15, dhw_time[:, loc])
        if isnothing(year)
            year_when_8DHW[loc] = 2100

            continue
        end
        year_when_8DHW[loc] =  dhw_time.timesteps[year]
    end

    if all(dhw_locs.sites .== context_layers.RME_UNIQUE_ID)
        context_layers[:, "$(GCM)_mean_dhw"] = dhw_locs
        context_layers[:, "$(GCM)_year_8DHW"] = Int64.(year_when_8DHW)
    else
        reordered_indices = indexin(context_layers.RME_UNIQUE_ID, String.(dhw_locs.sites))
        if all(dhw_locs.sites[reordered_indices] .== context_layers.RME_UNIQUE_ID) .& (any(.!isnothing.(reordered_indices)))
            context_layers[:, "$(GCM)_mean_dhw"] = dhw_locs[Int64.(reordered_indices)]
            context_layers[:, "$(GCM)_year_8DHW"] = Int64.(year_when_8DHW[reordered_indices])
        else
            @warn "Warning, mean DHW data not joined to context layers due to incorrect ordering"
        end
    end

    rel_cover = readcubedata(rs.total_relative_cover_median)
    evenness_data = readcubedata(rs.scaled_taxa_evenness)
    taxa_cover = readcubedata(rs.taxa_cover)
    cots_mortality = readcubedata(rs.cots_mortality)
    cyc_mortality = readcubedata(rs.cyc_mortality)
    dhw_mortality = readcubedata(rs.dhw_mortality)

    dhw_time = readcubedata(rs.negative_normal_dhw)

    dhw_cover_cor = zeros(Float64, length(context_layers.RME_UNIQUE_ID))
    dhw_evenness_cor = zeros(Float64, length(context_layers.RME_UNIQUE_ID))
    dhw_species1 = zeros(Float64, length(context_layers.RME_UNIQUE_ID))
    dhw_species2 = zeros(Float64, length(context_layers.RME_UNIQUE_ID))
    dhw_species3 = zeros(Float64, length(context_layers.RME_UNIQUE_ID))
    dhw_species4 = zeros(Float64, length(context_layers.RME_UNIQUE_ID))
    dhw_species5 = zeros(Float64, length(context_layers.RME_UNIQUE_ID))
    dhw_species6 = zeros(Float64, length(context_layers.RME_UNIQUE_ID))

    if all(context_layers.RME_UNIQUE_ID .== dhw_time.sites)
        for ind in 1:length(context_layers.RME_UNIQUE_ID)
            reef_total_cover = rel_cover[sites = ind]
            reef_evenness = evenness_data[sites = ind]
            reef_dhw = dhw_time[sites = ind]

            dhw_cover_cor[ind] = cross_correlation(reef_total_cover, reef_dhw, 0, CF, true)
            dhw_evenness_cor[ind] = cross_correlation(reef_evenness, reef_dhw, 0, CF, true)
            dhw_species1[ind] = cross_correlation(taxa_cover[sites = ind, species=1], reef_dhw, 0, CF, true)
            dhw_species2[ind] = cross_correlation(taxa_cover[sites = ind, species=2], reef_dhw, 0, CF, true)
            dhw_species3[ind] = cross_correlation(taxa_cover[sites = ind, species=3], reef_dhw, 0, CF, true)
            dhw_species4[ind] = cross_correlation(taxa_cover[sites = ind, species=4], reef_dhw, 0, CF, true)
            dhw_species5[ind] = cross_correlation(taxa_cover[sites = ind, species=5], reef_dhw, 0, CF, true)
            dhw_species6[ind] = cross_correlation(taxa_cover[sites = ind, species=6], reef_dhw, 0, CF, true)
            # for taxa in 1:6
            #     reef_cover = taxa_cover[sites = At(reef.RME_UNIQUE_ID), species= At(taxa)]

            #     reef["$(GCM)_dhw_species$(taxa)_cover_cor"] = cross_correlation(reef_cover, reef_dhw, 0, CF, true)
            # end

        end
    else
        throw("context_layers and dhw_time orders do not match")
    end

    context_layers[!, "$(GCM)_dhw_cover_cor"] = dhw_cover_cor
    context_layers[!, "$(GCM)_dhw_evenness_cor"] = dhw_evenness_cor
    context_layers[!, "$(GCM)_dhw_species1_cover_cor"] = dhw_species1
    context_layers[!, "$(GCM)_dhw_species2_cover_cor"] = dhw_species2
    context_layers[!, "$(GCM)_dhw_species3_cover_cor"] = dhw_species3
    context_layers[!, "$(GCM)_dhw_species4_cover_cor"] = dhw_species4
    context_layers[!, "$(GCM)_dhw_species5_cover_cor"] = dhw_species5
    context_layers[!, "$(GCM)_dhw_species6_cover_cor"] = dhw_species6

    if all(context_layers.RME_UNIQUE_ID .== rs.cots_mortality.sites)
        context_layers[:, "$(GCM)_mean_cots_mortality"]  = vec(mapslices(mean, cots_mortality.data; dims=1)')
        context_layers[:, "$(GCM)_mean_cyc_mortality"] = vec(mapslices(mean, cyc_mortality.data; dims=1)')
        context_layers[:, "$(GCM)_mean_dhw_mortality"] = vec(mapslices(mean, dhw_mortality.data; dims=1)')
    end

    context_layers[!, "$(GCM)_year_8DHW"] .= convert.(Float64, context_layers[:, "$(GCM)_year_8DHW"])
end

# 5. Calculate the number of reefs in each bioregion
context_layers.bioregion_count .= 1
for reef in eachrow(context_layers)
    reef.bioregion_count = count(context_layers.bioregion .== [reef.bioregion])
end

# 6. Calculate eigenvector_centrality scores for each reef - GBR-wide connections weighted by subregions
mean_conn, std_conn = load_connectivity(conn_path, context_layers.RME_UNIQUE_ID)

mean_conn[diagind(mean_conn)] .= 0.0

conn_whole_gbr = connectivity_scoring(mean_conn)
conn_mgmt_area = connectivity_scoring(mean_conn; gdf=context_layers, context_layer=:management_area, conn_col_name=:conn_score_mgmt, by_layer=true)
conn_subr_area = connectivity_scoring(mean_conn; gdf=context_layers, context_layer=:closest_port, conn_col_name=:conn_score_subr, by_layer=true)
conn_bior_area = connectivity_scoring(mean_conn; gdf=context_layers, context_layer=:bioregion, conn_col_name=:conn_score_bior, by_layer=true)

if all(conn_whole_gbr.RME_UNIQUE_ID .== context_layers.RME_UNIQUE_ID)
    context_layers.conn_score = conn_whole_gbr.conn_score
    context_layers.conn_score_mgmt = conn_mgmt_area.conn_score_mgmt
    context_layers.conn_score_subr = conn_subr_area.conn_score_subr
    context_layers.conn_score_bior = conn_bior_area.conn_score_bior
else
    @warn "Warning, connectivity data not joined to context layers due to incorrect ordering"
end

context_layers = weight_by_context(context_layers, :conn_score_mgmt, :management_area, :conn_ranked_mgmt)
context_layers = weight_by_context(context_layers, :conn_score_subr, :closest_port, :conn_ranked_subr)
context_layers = weight_by_context(context_layers, :conn_score_bior, :bioregion, :conn_ranked_bior)

# Format columns for writing to geopackage
context_layers.income_strength = Float64.(context_layers.income_strength)
context_layers.income_count .= convert.(Int64, context_layers.income_count)
context_layers.income_comb .= convert.(Float64, context_layers.income_comb)
context_layers.out_strength .= convert.(Float64, context_layers.out_strength)
context_layers.out_count .= convert.(Int64, context_layers.out_count)
context_layers.out_comb .= convert.(Float64, context_layers.out_comb)
context_layers.total_strength .= convert.(Float64, context_layers.total_strength)
context_layers.total_count .= convert.(Int64, context_layers.total_count)
context_layers.total_comb .= convert.(Float64, context_layers.total_comb)
context_layers.so_to_si .= convert.(Float64, context_layers.so_to_si)
# context_layers.dhw_cor = ifelse.(isnan.(context_layers.dhw_cor), 0.0, context_layers.dhw_cor)
# context_layers.area .= convert.(Float64, context_layers.area)
context_layers.conn_score .= convert.(Float64, context_layers.conn_score)
context_layers.conn_score_mgmt .= convert.(Float64, context_layers.conn_score_mgmt)
context_layers.conn_score_subr .= convert.(Float64, context_layers.conn_score_subr)
context_layers.conn_score_bior .= convert.(Float64, context_layers.conn_score_bior)


GDF.write("../data/analysis_context_layers_rme_CF.gpkg", context_layers; crs=GFT.EPSG(7844), overwrite=true)
