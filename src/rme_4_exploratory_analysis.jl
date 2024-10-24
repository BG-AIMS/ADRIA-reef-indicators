"""
Access the domain connectivity and calculate the strength and number of outgoing/incoming
connections for each reef. Run immediately after 4_.jl to use the updated context_layers with
bioregions.
"""

include("common.jl")
include("plotting_functions.jl")

using DataFrames, Statistics, YAXArrays

using GLMakie, GeoMakie, GraphMakie

using GLM

using ADRIA

import GeoDataFrames as GDF

# Load updated context layer data
context_layers = GDF.read("../data/analysis_context_layers_rme.gpkg",)
analysis_levels = ["target_reefs_mgmt", "target_reefs_subr", "target_reefs_bior", "target_reefs_mgmt_evenness", "target_reefs_subr_evenness","target_reefs_bior_evenness"]

context_target_reefs = DataFrame(
    analysis_level = [],
    n_reefs = [],
    n_mgmt_regions = [],
    n_bioregions = [],
    average_so_to_si = [],
    average_total_strength = [],
    average_dhw = []
)

for analysis_level in analysis_levels
    context = context_layers[context_layers[:, analysis_level], :]
    n_reefs = size(context, 1)
    n_mgmt_regions = length(unique(context.management_area))
    n_bioregions = length(unique(context.bioregion))
    average_so_to_si = mean(context.so_to_si)
    average_total_strength = mean(context.total_strength)
    average_dhw = mean(context.mean_dhw)

    push!(context_target_reefs, [analysis_level,n_reefs, n_mgmt_regions, n_bioregions, average_so_to_si, average_total_strength, average_dhw])
end

# Looking at the 6 target_reefs_bior and 98 target_reefs_bior_evenness
# 4 of the 6 target reefs are also in the 98 target_reefs_bior_evenness

# Looking at 6 target_reefs_bior
bior_reefs = context_layers[context_layers.target_reefs_bior, :]
# 2/6 bioregion target reefs are COTS_priority N non-priority
# No reefs are EcoRRAP photogrammetry sites
# Cairns-Cooktown and Townsville-Whitsunday regions
# 2 Green zones, 3 Dark Blue zones, 1 Yellow zone
# 3 reefs in Eastern Kuku Yalanji IPA
# 4 bioregions (Outer barrier reefs n=260, Sheltered mid shelf reefs n=104, high continental island reefs n=126, NA n=273)
# so_to_si 3.28, 0.06, 0.91, 0.53, 0.59, 0.78
# mean_dhw 6.84, 5.98, 6, 10, 10, 12.9

# Looking at the 98 target_reefs_bior_evenness
bior_evenness_reefs = context_layers[context_layers.target_reefs_bior_evenness, :]
# 4/98 bioregion target reefs are COTS priority P targets, 1/98 are T targets
# No EcoRRAP photogrammetry reefs
# 11/98 FN, 17/98 CC, 45/98 TSV-Wht, 25/98 Mac-Cap
# 22/98 Green zones, 45/98 Dark Blue zones, 28 Yellow zones, 1 Light Blue Zone, 1 Pink zone, 1 Orange zone
# 4/98 in Wuthathi IPA, 6/98 in Eastern Kuku Yalanji IPA
# Bioregions: (
# Outer Barrier Reefs, Sheltered Mid Shelf Reefs, Outer Shelf Reefs, Exposed Mid Shelf Reefs,
# Far Northern Open Lagoon Reefs, High Continental Island Reefs, Coastal Southern Reefs,
# Strong Tidal Outer Shelf Reefs, Strong Tidal Mid Shelf Reefs (West),
# High Tidal Fringing Reefs, Strong Tidal Mid Shelf Reefs (East), Capricorn Bunker Mid Shelf Reefs
# )
# median so_to_si 0.88
# median mean_dhw 9.39

# Looking at all 3806 reefs
# median so_to_si 0.998
# median mean_dhw 11.38


# Coral Cover Analysis
# Load rel_cover data and identify reefs that crash in the first timeseries year

rs = open_dataset("../outputs/RME_result_stores/RME_SSP245_20reps/cover_and_evenness_2024_10_18.nc")
rel_cover = rs.total_relative_cover_median
#crash_reefs = collect(getAxis("sites", rel_cover).val)[rel_cover.data[2,:] .< 0.5]

# Remove crashing reefs
no_crash_gbr = context_layers[context_layers.UNIQUE_ID .∉ [crash_reefs],:]

# Filter to only include reefs that are within the same bioregion/closest_port subregion as target reefs
filtered_bior = no_crash_gbr[no_crash_gbr.bioregion .∈ [unique(no_crash_gbr[no_crash_gbr.target_reefs_bior, :bioregion])], :]
#filtered_subr = no_crash_gbr[no_crash_gbr.closest_port .∈ [unique(no_crash_gbr[no_crash_gbr.target_reefs_subr, :closest_port])], :]
filtered_bior = filtered_bior[filtered_bior.so_to_si .< quantile(filtered_bior.so_to_si, 0.95), :]

# Basic exploratory models of factors
#glm_allreefs = glm(@formula(target_reefs_bior ~ out_comb + dhw_cor ), no_crash_gbr, Binomial())
glm_bioregions = glm(@formula(target_reefs_bior ~ mean_dhw + so_to_si), filtered_bior, Binomial())
#glm_subregions = glm(@formula(target_reefs_subr ~ total_comb + so_to_si + initial_coral_cover + bioregion_count), filtered_subr, Binomial())

aic(glm_allreefs), aic(glm_bioregions), aic(glm_subregions)

glmm_form = @formula(target_reefs_bior ~ initial_coral_cover + so_to_si + mean_dhw + (1|bioregion))
glmm_fit = fit(MixedModel, glmm_form, filtered_bior, Binomial(), ProbitLink())
aic(glmm_fit) # seems to be an improvement when (1|bioregion) is used

bior_reefs_dhw = basic_target_reef_boxplot(
    categorical(filtered_bior.target_reefs_bior).refs,
    filtered_bior.mean_dhw;
    ylabel="mean DHW",
    title="Reefs mean DHW"
)
save("../figs/bior_reefs_mean_dhw.png", bior_reefs_dhw)

bior_reefs_so_to_si = basic_target_reef_boxplot(
    categorical(filtered_bior.target_reefs_bior).refs,
    filtered_bior.so_to_si;
    ylabel="Source to Sink Ratio (< 70)",
    title="Reefs Source to Sink Ratio"
)
save("../figs/bior_reefs_so_to_si.png", bior_reefs_so_to_si)

bior_reefs_conn_centrality = basic_target_reef_boxplot(
    categorical(filtered_bior.target_reefs_bior).refs,
    filtered_bior.conn_score;
    ylabel="Connectivity Centrality Score",
    title="Reefs Connectivity Centrality Score"
)
save("../figs/bior_reefs_conn_centrality.png", bior_reefs_conn_centrality)

bior_reefs_initial_cover = basic_target_reef_boxplot(
    categorical(filtered_bior.target_reefs_bior).refs,
    filtered_bior.initial_coral_cover;
    ylabel="Initial Coral Cover"
)
save("../figs/bior_reefs_initial_coral_cover.png", bior_reefs_initial_cover)

rel_cover = rs.total_relative_cover_median
df = DataFrame(rel_cover.data, collect(getAxis("sites", rel_cover).val))
df.year = [string(i) for i in 1:size(df,1)]
select!(df, :year, Not(:year))

data = permutedims(df, 1, "UNIQUE_ID")
data = data[data.UNIQUE_ID .∈ [filtered_bior.UNIQUE_ID],:]
data = leftjoin(data, filtered_bior[:, [:UNIQUE_ID, :target_reefs_bior, :bioregion, :management_area]], on=:UNIQUE_ID)
data.target_reefs_bior = ifelse.(data.target_reefs_bior, "bellwether", "non-bellwether")

timeseries_bior_reefs = timeseries_plot(data, :target_reefs_bior, (1,78), "Total relative coral cover (%)", 0.05)
save("../figs/bior_reefs_timeseries_plot.png", timeseries_bior_reefs)

lagged_ts_bior_reefs = lagged_timeseries_plot(data, :target_reefs_bior, (1,78), "Year", "Total relative coral cover", 0.05, 4)
save("../figs/bior_reefs_lagged_timeseries.png", lagged_ts_bior_reefs)

# Evenness

evenness = rs.scaled_taxa_evenness
#crash_reefs = collect(getAxis("sites", rel_cover).val)[rel_cover.data[2,:] .< 0.5]

# Remove crashing reefs
# no_crash_gbr = context_layers[context_layers.UNIQUE_ID .∉ [crash_reefs],:]

# Filter to only include reefs that are within the same bioregion/closest_port subregion as target reefs
filtered_bior_evenness = context_layers[context_layers.bioregion .∈ [unique(context_layers[context_layers.target_reefs_bior_evenness, :bioregion])], :]
filtered_bior_evenness = filtered_bior_evenness[filtered_bior_evenness.so_to_si .< quantile(filtered_bior_evenness.so_to_si, 0.95), :]


# Basic exploratory models of factors
#glm_allreefs = glm(@formula(target_reefs_bior ~ out_comb + dhw_cor ), no_crash_gbr, Binomial())
glm_bioregions = glm(@formula(target_reefs_bior_evenness ~ so_to_si + mean_dhw + bioregion + 0), filtered_bior_evenness, Binomial())
#glm_subregions = glm(@formula(target_reefs_subr ~ total_comb + so_to_si + initial_coral_cover + bioregion_count), filtered_subr, Binomial())

#aic(glm_allreefs), aic(glm_bioregions), aic(glm_subregions)

glmm_form = @formula(target_reefs_bior_evenness ~ initial_coral_cover + so_to_si + mean_dhw + (1|bioregion))
glmm_fit = fit(MixedModel, glmm_form, filtered_bior_evenness, Binomial(), ProbitLink())
aic(glmm_fit) # seems to be an improvement when (1|bioregion) is used



bior_evenness_dhw = basic_target_reef_boxplot(
    categorical(filtered_bior_evenness.target_reefs_bior_evenness).refs,
    filtered_bior_evenness.mean_dhw;
    ylabel="mean DHW",
    title="Reef mean DHW"
)
save("../figs/bior_evenness_dhw.png", bior_evenness_dhw)

bior_evenness_so_to_si = basic_target_reef_boxplot(
    categorical(filtered_bior_evenness.target_reefs_bior_evenness).refs,
    filtered_bior_evenness.so_to_si, [0,1,2,3,4,5,6,7,8,9,10];
    ylabel="Source to Sink Ratio (< 10)"
)
save("../figs/bior_evenness_so_to_si.png", bior_evenness_so_to_si)

bior_evenness_conn_centrality = basic_target_reef_boxplot(
    categorical(filtered_bior_evenness.target_reefs_bior_evenness).refs,
    filtered_bior_evenness.conn_score;
    ylabel="Connectivity Centrality Score"
)
save("../figs/bior_evenness_conn_centrality.png", bior_evenness_conn_centrality)

bior_evenness_total_conn_strength = basic_target_reef_boxplot(
    categorical(filtered_bior_evenness.target_reefs_bior_evenness).refs,
    filtered_bior_evenness.total_strength;
    ylabel="Total Connectivity Strength"
)
save("../figs/bior_evenness_total_conn_strength.png", bior_evenness_total_conn_strength)

bior_evenness_initial_cover = basic_target_reef_boxplot(
    categorical(filtered_bior_evenness.target_reefs_bior_evenness).refs,
    filtered_bior_evenness.initial_coral_cover;
    ylabel="Initial Coral Cover"
)
save("../figs/bior_evenness_initial_coral_cover.png", bior_evenness_initial_cover)


evenness = rs.scaled_taxa_evenness
df = DataFrame(evenness.data, collect(getAxis("sites", evenness).val))
df.year = [string(i) for i in 1:size(df,1)]
select!(df, :year, Not(:year))

data = permutedims(df, 1, "RME_UNIQUE_ID")
data = data[data.RME_UNIQUE_ID .∈ [filtered_bior_evenness.RME_UNIQUE_ID],:]
data = leftjoin(data, filtered_bior_evenness[:, [:RME_UNIQUE_ID, :target_reefs_bior_evenness, :bioregion, :management_area]], on=:RME_UNIQUE_ID)
data = dropmissing(data)
data.target_reefs_bior_evenness = ifelse.(data.target_reefs_bior_evenness, "bellwether", "non-bellwether")

timeseries_bior_evenness = timeseries_plot(data, :target_reefs_bior_evenness, (1,78), "Scaled Taxa Evenness", 0.05)
save("../figs/bior_evenness_timeseries_plot.png", timeseries_bior_evenness)

lagged_ts_bior_evenness = lagged_timeseries_plot(data, :target_reefs_bior_evenness, (1,78), "Year", "Scaled_Taxa Evenness", 0.05, 5)
save("../figs/bior_evenness_lagged_timeseries.png", lagged_ts_bior_evenness)

#data = data[indexin(context_layers.RME_UNIQUE_ID, data.RME_UNIQUE_ID), :]
all_zeros = data[[all(collect(dat[2:79]) .== 0.0) for dat in eachrow(data)], :RME_UNIQUE_ID]
context_layers_no_zeros = context_layers[context_layers.RME_UNIQUE_ID .∉ [all_zeros], :]

all_zeros = [all(collect(dat[2:79]) .== 0.0) for dat in eachrow(data)]
bior_reefs_no_zeros = filtered_bior_evenness[.!all_zeros, :]
bior_reefs_no_zeros = reefs_no_zeros[reefs_no_zeros.lag5_bior_evenness .!= 0.0, :]
#context_layers_all_cor = context_layers[context_layers.bioregion .∉ [unique((context_layers[context_layers.lag5_bior_evenness .== 0.0, :bioregion]))], :]
summarise_scores_evenness = combine(DataFrames.groupby(reefs_no_zeros, :bioregion), :lag5_bior_evenness => mean, :total_strength => mean, :ReefMod_area_m2 => mean, :conn_score => mean)
highest_bioregion = summarise_scores_evenness.bioregion[argmax(summarise_scores_evenness.lag5_bior_evenness_mean)]
lowest_bioregion = summarise_scores_evenness.bioregion[argmin(summarise_scores_evenness[summarise_scores_evenness.lag5_bior_evenness_mean .> 0.0, :lag5_bior_evenness_mean])]

data = permutedims(df, 1, "RME_UNIQUE_ID")
data = data[data.UNIQUE_ID .∈ [context_layers_no_zeros[context_layers_no_zeros.bioregion .== "Coral Sea Swains-Northern Reefs", :].UNIQUE_ID],:]
data = leftjoin(data, context_layers_no_zeros[:, [:UNIQUE_ID, :target_reefs_bior_evenness, :bioregion, :management_area, :lag5_bior_evenness]], on=:UNIQUE_ID)
data = dropmissing(data)
data.target_reefs_bior_evenness = ifelse.(data.target_reefs_bior_evenness, "bellwether", "non-bellwether")

lagged_ts = plot_lines(data, :UNIQUE_ID, 1, 78, "Year", "Scaled Taxa Evenness", 0.5)
