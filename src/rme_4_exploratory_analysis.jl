"""
Access the domain connectivity and calculate the strength and number of outgoing/incoming
connections for each reef. Run immediately after 4_.jl to use the updated context_layers with
bioregions.
"""

include("common.jl")
include("plotting_functions.jl")
include("custom_RMEDomain_dataloading.jl")

using DataFrames, Statistics, YAXArrays

using
    #GLMakie,
    GeoMakie,
    GraphMakie,
    CairoMakie

using GLM, MixedModels

using ADRIA

import GeoDataFrames as GDF

CairoMakie.activate!(type = "png")


# Load updated context layer data
context_layers = GDF.read("../data/analysis_context_layers_rme.gpkg",)
context_layers = context_layers[context_layers.bioregion .!= "NA", :]
context_layers.target_reefs_bior_cat = ifelse.(context_layers.target_reefs_bior, "bellwether", "non-bellwether")
context_layers.target_reefs_bior_evenness_cat = ifelse.(context_layers.target_reefs_bior_evenness, "bellwether", "non-bellwether")

# analysis_levels = ["target_reefs_mgmt", "target_reefs_subr", "target_reefs_bior", "target_reefs_mgmt_evenness", "target_reefs_subr_evenness","target_reefs_bior_evenness"]

# context_target_reefs = DataFrame(
#     analysis_level = [],
#     n_reefs = [],
#     n_mgmt_regions = [],
#     n_bioregions = [],
#     average_so_to_si = [],
#     average_total_strength = [],
#     average_dhw = []
# )

# for analysis_level in analysis_levels
#     context = context_layers[context_layers[:, analysis_level], :]
#     n_reefs = size(context, 1)
#     n_mgmt_regions = length(unique(context.management_area))
#     n_bioregions = length(unique(context.bioregion))
#     average_so_to_si = mean(context.so_to_si)
#     average_total_strength = mean(context.total_strength)
#     average_dhw = mean(context.mean_dhw)

#     push!(context_target_reefs, [analysis_level,n_reefs, n_mgmt_regions, n_bioregions, average_so_to_si, average_total_strength, average_dhw])
# end

# Looking at the 21 target_reefs_bior and 776 target_reefs_bior_evenness
# 11 of the 21 target reefs are also in the 776 target_reefs_bior_evenness

# Looking at 21 target_reefs_bior
bior_reefs = context_layers[Bool.(context_layers.target_reefs_bior), :]
# 2/21 bioregion target reefs are COTS_priority T - target
# No reefs are EcoRRAP photogrammetry sites
# Cairns-Cooktown, Townsville-Whitsunday and Mackay-Capricorn regions
# 7 Green zones, 9 Dark Blue zones, 4 Yellow zone, 1 Orange zone
# 5 reefs in Eastern Kuku Yalanji IPA
# 8 bioregions (
# "Sheltered Mid Shelf Reefs", "Outer Barrier Reefs", "Outer Shelf Reefs", "High Continental Island Reefs",
# "Exposed Mid Shelf Reefs", "NA", "High Tidal Fringing Reefs", "Coastal Southern Reefs")

# Looking at the 776 target_reefs_bior_evenness
bior_evenness_reefs = context_layers[context_layers.target_reefs_bior_evenness, :]
# 8/776 bioregion target reefs are COTS priority T targets, 26 are P targets
# No EcoRRAP photogrammetry reefs
# 32/776 FN, 65/776 CC, 260/776 TSV-Wht, 419/776 Mac-Cap
# 239/776 Green zones, 474/776 Dark Blue zones, 20/776 Yellow zones, 29/776 Light Blue Zone, 10/776 Pink zone, 3/776 Orange zone, 1/776 Olive Green Zone
# 6/776 in Wuthathi IPA, 17/776 in Eastern Kuku Yalanji IPA, 4/776 Girringun
# Bioregions: (
# "Outer Barrier Reefs", "Far Northern Protected Mid Shelf Reefs and Shoals", "Far Northern Outer Mid Shelf Reefs",
# "Far Northern Open Lagoon Reefs", "Sheltered Mid Shelf Reefs", "Outer Shelf Reefs", "Exposed Mid Shelf Reefs",
# "Northern Open Lagoon Reefs", "High Continental Island Reefs", "Strong Tidal Outer Shelf Reefs",
# "Strong Tidal Mid Shelf Reefs (West)", "Coral Sea Swains-Northern Reefs", "High Tidal Fringing Reefs",
# "Strong Tidal Inner Mid Shelf Reefs", "Strong Tidal Mid Shelf Reefs (East)", "Swains Mid Reefs",
# "Coastal Southern Fringing Reefs", "Capricorn Bunker Mid Shelf Reefs"
# )


# Coral Cover Analysis
# Load rel_cover data and identify reefs that crash in the first timeseries year

rs = open_dataset("../outputs/RME_result_stores/RME_SSP245_300reps/cover_and_evenness_2024_11_14.nc")
rel_cover = rs.total_relative_cover_median
#crash_reefs = collect(getAxis("sites", rel_cover).val)[rel_cover.data[2,:] .< 0.5]

# Remove crashing reefs
#no_crash_gbr = context_layers[context_layers.UNIQUE_ID .∉ [crash_reefs],:]

# Filter to only include reefs that are within the same bioregion/closest_port subregion as target reefs
filtered_bior = context_layers[context_layers.bioregion .∈ [unique(context_layers[context_layers.target_reefs_bior_cat .== "bellwether", :bioregion])], :]
filtered_bior = filtered_bior[filtered_bior.so_to_si .< quantile(filtered_bior.so_to_si, 0.95), :]
filtered_bior = filtered_bior[filtered_bior.conn_score .< quantile(filtered_bior.conn_score, 0.95), :]
filtered_bior = filtered_bior[filtered_bior.bioregion .∈ [unique(filtered_bior[filtered_bior.target_reefs_bior_cat .== "bellwether", :bioregion])], :]

# Basic exploratory models of factors
#glm_allreefs = glm(@formula(target_reefs_bior ~ out_comb + dhw_cor ), no_crash_gbr, Binomial())
glm_bioregions = glm(@formula(target_reefs_bior ~ mean_dhw + so_to_si + conn_score + total_strength + initial_coral_cover + dhw_cover_cor + mean_cots_mortality), filtered_bior, Binomial())
#glm_subregions = glm(@formula(target_reefs_subr ~ total_comb + so_to_si + initial_coral_cover + bioregion_count), filtered_subr, Binomial())

# aic(glm_allreefs), aic(glm_bioregions), aic(glm_subregions)

glmm_form = @formula(
    target_reefs_bior ~
    mean_dhw +
    so_to_si +
    conn_score +
    total_strength +
    initial_coral_cover +
    initial_evenness +
    dhw_cover_cor +
    mean_cots_mortality +
    (1|bioregion)
)
glmm_fit = fit(MixedModel, glmm_form, filtered_bior, Bernoulli())

# aic(glmm_fit) # seems to be an improvement when (1|bioregion) is used

# Plotting bioregion raincloud plots for each variable
mean_dhw_raincloud = bioregion_grouped_boxplots(
    filtered_bior,
    :target_reefs_bior_cat,
    :bioregion, :mean_dhw, 4,
    ylabel="mean DHW"
);
save("../figs/mean_dhw_raincloud_cover.png", mean_dhw_raincloud)

so_to_si_raincloud = bioregion_grouped_boxplots(
    filtered_bior,
    :target_reefs_bior_cat,
    :bioregion, :so_to_si, 4,
    ylabel="Source to sink ratio"
);
save("../figs/so_to_si_raincloud_cover.png", so_to_si_raincloud)

conn_score_raincloud = bioregion_grouped_boxplots(
    filtered_bior,
    :target_reefs_bior_cat,
    :bioregion, :conn_score, 4,
    ylabel="Connectivity eigenvector centrality"
);
save("../figs/conn_score_raincloud_cover.png", conn_score_raincloud)

total_strength_raincloud = bioregion_grouped_boxplots(
    filtered_bior,
    :target_reefs_bior_cat,
    :bioregion, :total_strength, 4,
    ylabel="Total connectivity strength"
);
save("../figs/total_strength_raincloud_cover.png", total_strength_raincloud)

initial_cover_raincloud = bioregion_grouped_boxplots(
    filtered_bior,
    :target_reefs_bior_cat,
    :bioregion, :initial_coral_cover, 4,
    ylabel="Initial coral cover (%)"
);
save("../figs/initial_cover_raincloud_cover.png", initial_cover_raincloud)

dhw_cover_cor_raincloud = bioregion_grouped_boxplots(
    filtered_bior,
    :target_reefs_bior_cat,
    :bioregion, :dhw_cover_cor, 4,
    ylabel="DHW and cover correlation"
);
save("../figs/dhw_cover_cor_raincloud_cover.png", dhw_cover_cor_raincloud)

mean_COTS_mortality_raincloud = bioregion_grouped_boxplots(
    filtered_bior,
    :target_reefs_bior_cat,
    :bioregion, :mean_cots_mortality, 4,
    ylabel="CoTS mortality"
)
save("../figs/mean_COTS_mortality_raincloud_cover.png", mean_COTS_mortality_raincloud)

rel_cover = rs.total_relative_cover_median
rel_cover_reefs = extract_timeseries(rel_cover, context_layers, [:bioregion, :lag4_bior, :target_reefs_bior_cat, :management_area])
rel_cover_reefs = rel_cover_reefs[.!(Bool.(context_layers.bior_cover_qc_flag)), :]
rel_cover_reefs = rel_cover_reefs[rel_cover_reefs.RME_UNIQUE_ID .∈ [filtered_bior.RME_UNIQUE_ID], :]

timeseries_bior_reefs = timeseries_plot(rel_cover_reefs, :target_reefs_bior_cat, (1,50), "Total relative coral cover (%)", 0.05)
save("../figs/bior_reefs_timeseries_plot.png", timeseries_bior_reefs)

lagged_ts_bior_reefs = lagged_timeseries_plot(rel_cover_reefs, :target_reefs_bior_cat, (1,50), "Year", "Total relative coral cover (%)", 0.1, 4)
save("../figs/bior_reefs_lagged_timeseries.png", lagged_ts_bior_reefs)

lagged_ts_combined_bior = lagged_timeseries_plot_combined(rel_cover_reefs, :target_reefs_bior_cat, (1,50), "Year", "Total relative coral cover (%)", 0.1, 4)
save("../figs/bior_reefs_combined_timeseries.png", lagged_ts_combined_bior)

bioregion_grouped_timeseries, bior_reefs_bioregion_map = bioregion_grouped_lagged_timeseries(rel_cover_reefs, :target_reefs_bior_cat, (1,50), :bioregion, 0, 4; ylabel="Total relative coral cover (%)")
save("../figs/bioregion_cover_timeseries.png", bioregion_grouped_timeseries)
bioregion_grouped_timeseries_lagged = bioregion_grouped_lagged_timeseries(rel_cover_reefs, :target_reefs_bior_cat, (1,50), :bioregion, 4, 4; ylabel="Total relative coral cover (%)")[1]
save("../figs/bioregion_cover_timeseries_lagged.png", bioregion_grouped_timeseries_lagged)

dhw_data = load_DHW("../../RME/rme_ml_2024_06_13/data_files/", "SSP245")
dhw_ts = Float64.(mapslices(median, dhw_data, dims=[:scenarios]))
dhw_ts = dhw_ts[locs = At(context_layers.RME_GBRMPA_ID)]
axlist = (
    Dim{:timesteps}(2022:2099),
    Dim{:sites}(context_layers.RME_UNIQUE_ID)
)
dhw_ts = rebuild(dhw_ts, dims=axlist)

bior_reefs_dhw = extract_timeseries(dhw_ts, filtered_bior, [:bioregion, :lag4_bior, :target_reefs_bior_cat])
bior_reefs_dhw_plot = timeseries_plot(bior_reefs_dhw, :target_reefs_bior_cat, (1,50), "median scenario DHW", 0.05)
save("../figs/evenness_bellwether_dhw_timeseries.png", bior_evenness_dhw)

CC_filtered = filtered_bior[filtered_bior.management_area_short .== "Cairns-Cooktown", :]
cover_CC_map = plot_map_bellwether_reefs(CC_filtered, :bioregion, :target_reefs_bior_cat)[1]
save("../figs/cover_CairnsCooktown_map.png", cover_CC_map)

TSV_filtered = filtered_bior[filtered_bior.management_area_short .== "Townsville-Whitsunday", :]
cover_TSV_map = plot_map_bellwether_reefs(TSV_filtered, :bioregion, :target_reefs_bior_cat)[1]
save("../figs/cover_TownsvilleWhitsunday_map.png", cover_TSV_map)

MCap_filtered = filtered_bior[filtered_bior.management_area_short .== "Mackay-Capricorn", :]
cover_MCap_map = plot_map_bellwether_reefs(MCap_filtered, :bioregion, :target_reefs_bior_cat)[1]
save("../figs/cover_MackayCapricorn_map.png", cover_MCap_map)


# Evenness
filtered_bior_evenness = context_layers[context_layers.bioregion .∈ [unique(context_layers[context_layers.target_reefs_bior_evenness_cat .== "bellwether", :bioregion])], :]
filtered_bior_evenness = filtered_bior_evenness[filtered_bior_evenness.so_to_si .< quantile(filtered_bior_evenness.so_to_si, 0.95), :]
filtered_bior_evenness = filtered_bior_evenness[filtered_bior_evenness.conn_score .< quantile(filtered_bior_evenness.conn_score, 0.95), :]
filtered_bior_evenness = filtered_bior_evenness[filtered_bior_evenness.bioregion .∈ [unique(filtered_bior_evenness[filtered_bior_evenness.target_reefs_bior_evenness_cat .== "bellwether", :bioregion])], :]

evenness = rs.scaled_taxa_evenness
evenness_ts = extract_timeseries(evenness, context_layers, [:bioregion, :lag4_bior_evenness, :target_reefs_bior_evenness_cat, :management_area])

# Remove reefs that have a qc flag indicating they have NaN/0 evenness values for every year
filtered_bior_evenness_no_zeros = filtered_bior_evenness[.!(Bool.(filtered_bior_evenness.bior_evenness_qc_flag)) , :]
evenness_no_zeros = evenness_ts[indexin(filtered_bior_evenness_no_zeros.RME_UNIQUE_ID, evenness_ts.RME_UNIQUE_ID), :]


# Basic exploratory models of factors

glm_bioregions = glm(@formula(target_reefs_bior_evenness ~ mean_dhw + so_to_si + conn_score + total_strength + initial_coral_cover + dhw_evenness_cor + mean_cots_mortality), filtered_bior_evenness_no_zeros, Binomial(), LogitLink())


glmm_form = @formula(target_reefs_bior_evenness ~ mean_dhw + so_to_si + conn_score + total_strength + initial_coral_cover + dhw_evenness_cor + mean_cots_mortality + (1|bioregion))
glmm_fit = fit(MixedModel, glmm_form, filtered_bior_evenness_no_zeros, Bernoulli())
aic(glmm_fit), aic(glm_bioregions) # seems to be an improvement when (1|bioregion) is used

# Plotting bioregion raincloud plots for each variable
mean_dhw_raincloud_evenness = bioregion_grouped_boxplots(
    filtered_bior_evenness,
    :target_reefs_bior_evenness_cat,
    :bioregion, :mean_dhw, 5,
    ylabel="mean DHW"
);
save("../figs/mean_dhw_raincloud_evenness.png", mean_dhw_raincloud_evenness)

so_to_si_raincloud_evenness = bioregion_grouped_boxplots(
    filtered_bior_evenness,
    :target_reefs_bior_evenness_cat,
    :bioregion, :so_to_si, 5,
    ylabel="Source to sink ratio"
);
save("../figs/so_to_si_raincloud_evenness.png", so_to_si_raincloud_evenness)

conn_score_raincloud_evenness = bioregion_grouped_boxplots(
    filtered_bior_evenness,
    :target_reefs_bior_evenness_cat,
    :bioregion, :conn_score, 5,
    ylabel="Connectivity eigenvector centrality"
);
save("../figs/conn_score_raincloud_evenness.png", conn_score_raincloud_evenness)

total_strength_raincloud_evenness = bioregion_grouped_boxplots(
    filtered_bior_evenness,
    :target_reefs_bior_evenness_cat,
    :bioregion, :total_strength, 5,
    ylabel="Total connectivity strength"
);
save("../figs/total_strength_raincloud_evenness.png", total_strength_raincloud_evenness)

initial_cover_raincloud_evenness = bioregion_grouped_boxplots(
    filtered_bior_evenness,
    :target_reefs_bior_evenness_cat,
    :bioregion, :initial_coral_cover, 5,
    ylabel="Initial coral cover (%)"
);
save("../figs/initial_cover_raincloud_evenness.png", initial_cover_raincloud_evenness)

dhw_evenness_cor_raincloud = bioregion_grouped_boxplots(
    filtered_bior_evenness_no_zeros,
    :target_reefs_bior_evenness_cat,
    :bioregion, :dhw_evenness_cor, 5,
    ylabel="Evenness and DHW correlation",
);
save("../figs/dhw_evenness_cor_raincloud_evenness.png", dhw_evenness_cor_raincloud)

mean_COTS_mortality_raincloud_evenness = bioregion_grouped_boxplots(
    filtered_bior_evenness_no_zeros,
    :target_reefs_bior_evenness_cat,
    :bioregion, :mean_cots_mortality, 5,
    ylabel="CoTS mortality",
);
save("../figs/mean_COTS_mortality_raincloud_evenness.png", mean_COTS_mortality_raincloud_evenness)


timeseries_bior_evenness = timeseries_plot(evenness_no_zeros, :target_reefs_bior_evenness_cat, (1,50), "Scaled Taxa Evenness", 0.05)
save("../figs/bior_evenness_timeseries_plot.png", timeseries_bior_evenness)

lagged_ts_bior_evenness = lagged_timeseries_plot(evenness_no_zeros, :target_reefs_bior_evenness_cat, (1,50), "Year", "Scaled Taxa Evenness", 0.05, 4)
save("../figs/bior_evenness_lagged_timeseries.png", lagged_ts_bior_evenness)

combined_ts_bior_evenness = lagged_timeseries_plot_combined(evenness_no_zeros, :target_reefs_bior_evenness_cat, (1,50), "Year", "Scaled taxa evenness", 0.1, 4)
save("../figs/bior_evenness_combined_timeseries.png", combined_ts_bior_evenness)

bioregion_grouped_timeseries_evenness, bior_evenness_group_map = bioregion_grouped_lagged_timeseries(evenness_no_zeros, :target_reefs_bior_evenness_cat, (1,50), :bioregion, 0, 5; ylabel="Cover-scaled taxa evenness")
save("../figs/bioregion_evenness_timeseries.png", bioregion_grouped_timeseries_evenness)
bioregion_grouped_timeseries_lagged_evenness = bioregion_grouped_lagged_timeseries(evenness_no_zeros, :target_reefs_bior_evenness_cat, (1,50), :bioregion, 4, 5; ylabel="Cover-scaled taxa evenness")[1]
save("../figs/bioregion_evenness_timeseries_lagged.png", bioregion_grouped_timeseries_lagged_evenness)
# # Mean_DHW plot
# mean_dhw_evenness_plot = explore_regression_scatter_plots(
#     filtered_bior_evenness_no_zeros.mean_dhw,
#     filtered_bior_evenness_no_zeros.lag4_bior_evenness .^2;
#     color = categorical(filtered_bior_evenness_no_zeros.target_reefs_bior_evenness).refs,
#     xlab="mean DHW",
#     ylab="(correlation at lag 5) ^2",
# )
# save("../figs/mean_dhw_evenness_correlation_scatter.png", mean_dhw_evenness_plot)

# initial_cover_plot = explore_regression_scatter_plots(
#     filtered_bior_evenness_no_zeros.initial_coral_cover,
#     filtered_bior_evenness_no_zeros.lag4_bior_evenness
# )

# summarise_scores_evenness = combine(DataFrames.groupby(all_reefs_no_zeros, :bioregion), :lag4_bior_evenness => mean, :total_strength => mean, :ReefMod_area_m2 => mean, :conn_score => mean)
# highest_bioregion = summarise_scores_evenness.bioregion[argmax(summarise_scores_evenness.lag4_bior_evenness_mean)]
# lowest_bioregion = summarise_scores_evenness.bioregion[argmin(summarise_scores_evenness[summarise_scores_evenness.lag4_bior_evenness_mean .> 0.0, :lag4_bior_evenness_mean])]


# looking at the timing of dhw stress
bior_evenness_dhw = extract_timeseries(dhw_ts, filtered_bior_evenness, [:bioregion, :lag4_bior_evenness, :target_reefs_bior_evenness_cat])
bior_evenness_dhw = timeseries_plot(bior_reefs_dhw, :target_reefs_bior_evenness_cat, (1,50), "median scenario DHW", 0.05)
save("../figs/evenness_bellwether_dhw_timeseries.png", bior_evenness_dhw)

FN_filtered_evenness = filtered_bior_evenness_no_zeros[filtered_bior_evenness_no_zeros.management_area_short .== "FarNorthern", :]
evenness_FN_map = plot_map_bellwether_reefs(FN_filtered_evenness, :bioregion, :target_reefs_bior_evenness_cat)[1]
save("../figs/evenness_FarNorthern_map.png", evenness_FN_map)

CC_filtered_evenness = filtered_bior_evenness_no_zeros[filtered_bior_evenness_no_zeros.management_area_short .== "Cairns-Cooktown", :]
evenness_CC_map = plot_map_bellwether_reefs(CC_filtered_evenness, :bioregion, :target_reefs_bior_evenness_cat)[1]
save("../figs/evenness_CairnsCooktown_map.png", evenness_CC_map)

TSV_filtered_evenness = filtered_bior_evenness_no_zeros[filtered_bior_evenness_no_zeros.management_area_short .== "Townsville-Whitsunday", :]
evenness_TSV_map = plot_map_bellwether_reefs(TSV_filtered_evenness, :bioregion, :target_reefs_bior_evenness_cat)[1]
save("../figs/evenness_TownsvilleWhitsunday_map.png", evenness_TSV_map)

MCap_filtered_evenness = filtered_bior_evenness_no_zeros[filtered_bior_evenness_no_zeros.management_area_short .== "Mackay-Capricorn", :]
evenness_MCap_map = plot_map_bellwether_reefs(MCap_filtered_evenness, :bioregion, :target_reefs_bior_evenness_cat)[1]
save("../figs/evenness_MackayCapricorn_map.png", evenness_MCap_map)
