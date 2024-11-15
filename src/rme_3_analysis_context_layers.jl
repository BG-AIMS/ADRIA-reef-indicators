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

include("common.jl")
include("custom_RMEDomain_dataloading.jl")

# Load context layers with bioregions and target reefs
context_layers = GDF.read("../data/context_layers_targetted_rme.gpkg")
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

rs = open_dataset("../outputs/RME_result_stores/RME_SSP245_300reps/cover_and_evenness_2024_11_14.nc"; driver=:netcdf)
# 2. Attach initial coral cover data to context_layers
initial_cc = rs.total_relative_cover_median[1, :]
context_layers.initial_coral_cover = initial_cc.data

initial_evenness = rs.scaled_taxa_evenness[1, :]
context_layers.initial_evenness = initial_evenness

# 3. Attach DHW data to context_layers
# mean DHW data for sites
dhw_time_scens = load_DHW("../../RME/rme_ml_2024_06_13/data_files/", "SSP245")
dhw_time = Float64.(mapslices(median, dhw_time_scens, dims=[:scenarios]))
dhw_locs = Float64.(mapslices(mean, dhw_time[1:50, :], dims=[:timesteps])) # Select only the years in relevant time

year_when_8DHW = []
for loc in 1:size(dhw_time, 2)
    year = findfirst(x -> x > 7, dhw_time[:, loc])
    push!(year_when_8DHW, dhw_time.timesteps[year])
end

if all(dhw_locs.locs .== context_layers.RME_GBRMPA_ID)
    context_layers.mean_dhw = dhw_locs
    context_layers.year_8DHW = year_when_8DHW
else
    @warn "Warning, mean DHW data not joined to context layers due to incorrect ordering"
end

rel_cover = rs.total_relative_cover_median
evenness_data = rs.scaled_taxa_evenness

context_layers.dhw_cover_cor .= 0.0
context_layers.dhw_evenness_cor .= 0.0

for reef in eachrow(context_layers)
    reef_cover = rel_cover[sites = At(reef.RME_UNIQUE_ID)]
    reef_evenness = evenness_data[sites = At(reef.RME_UNIQUE_ID)]

    reef_dhw = dhw_time[locs = At(reef.RME_GBRMPA_ID)]
    reef_dhw.data .= -1 * reef_dhw.data

    cover_correlation = cross_correlation(reef_cover, reef_dhw, [0], true)
    evenness_correlation = cross_correlation(reef_evenness, reef_dhw, [0], true)
    reef.dhw_cover_cor = first(cover_correlation)
    reef.dhw_evenness_cor = first(evenness_correlation)
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

if all(context_layers.RME_UNIQUE_ID .== rs.cots_mortality.sites)
    context_layers.mean_cots_mortality = Float64.(mapslices(mean, rs.cots_mortality; dims=[:timesteps]))
    context_layers.mean_cyc_mortality = Float64.(mapslices(mean, rs.cyc_mortality; dims=[:timesteps]))
    context_layers.mean_dhw_mortality = Float64.(mapslices(mean, rs.dhw_mortality; dims=[:timesteps]))
end

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
context_layers.mean_dhw .= convert.(Float64, context_layers.mean_dhw)
# context_layers.dhw_cor = ifelse.(isnan.(context_layers.dhw_cor), 0.0, context_layers.dhw_cor)
# context_layers.area .= convert.(Float64, context_layers.area)
context_layers.conn_score .= convert.(Float64, context_layers.conn_score)
context_layers.conn_score_mgmt .= convert.(Float64, context_layers.conn_score_mgmt)
context_layers.conn_score_subr .= convert.(Float64, context_layers.conn_score_subr)
context_layers.conn_score_bior .= convert.(Float64, context_layers.conn_score_bior)
context_layers.year_8DHW = Float64.(context_layers.year_8DHW)

GDF.write("../data/analysis_context_layers_rme.gpkg", context_layers; crs=GFT.EPSG(7844), overwrite=true)
