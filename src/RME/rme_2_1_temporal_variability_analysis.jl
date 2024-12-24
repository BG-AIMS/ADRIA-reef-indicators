include("../common.jl")
GCM = gcms[1]
rs = open_dataset("../outputs/RME_result_stores/RME_SSP245_300reps/RME_SSP245_300reps_$(GCM)/cover_and_evenness_2024_12_18.nc"; driver=:netcdf)

# Calculate coral cover relative to initial starting cover
total_cover = rs.total_relative_cover_median
rel_cover = Float64.(mapslices(relative_site_cover, total_cover, dims=[:timesteps]))

# Filter out reefs with less than 5% coral cover in initial timesteps
icc = total_cover[1,:].data
initial_cover_indices = icc .> 5

# Calculate the split timeseries and temporal variability metrics for the relative coral cover timeseries
temporal_scores = [temporal_variability(naive_split_metric(@view(rel_cover[:, i]), 3, mean)) for i in 1:length(rel_cover.sites)]
hist(temporal_scores)

# Plot lines and the distribution of temporal variability scores
rel_cover_vec = [vec(rel_cover.data[:, i]) for i in 1:length(rel_cover.sites)]
sub_2_5 = initial_cover_indices .& (temporal_scores .< 2.5)

line_scores_zip = zip(rel_cover_vec[sub_2_5], temporal_scores[sub_2_5])
extremes = extrema(temporal_scores[sub_2_5])
lines_color((a, c)) = lines!(a; color=c, colorrange=extremes, alpha = 0.3)

fig = Figure()
Axis(fig[1,1])
map(lines_color, line_scores_zip)
Colorbar(fig[1,2], limits = extremes)

"""

"""
function cluster_temporal_variability(reef_clusters, temporal_variability_results, threshold)
    bellwether_reefs = zeros(Bool, length(reef_clusters))

    for (ind, cluster) in enumerate(reef_clusters)
        indices_not_target = findall(reef_clusters .== cluster)
        indices_not_target = indices_not_target[indices_not_target .!= ind]

        target_temporal_variability = temporal_variability_results[ind]
        non_target_temp_variability = median(temporal_variability_results[indices_not_target])

        bellwether_reefs[ind] = (non_target_temp_variability - target_temporal_variability) > threshold
    end

    return bellwether_reefs
end

context_layers = context_layers[temporal_scores .< 5, :]
context_layers.temporal_scores = temporal_scores[temporal_scores .< 5]

context_layers[:, "$(GCM)_target_reefs_bior"] = cluster_temporal_variability(context_layers.management_area, context_layers.temporal_scores, 0.25)
context_layers[:, "$(GCM)_target_reefs_bior_cat"] = ifelse.(context_layers[:, "$(GCM)_target_reefs_bior"], "bellwether", "non-bellwether")

context_layers.gbr .= "gbr"
bellwethers_gbr = cluster_temporal_variability(context_layers.gbr, context_layers.temporal_scores, 0.1)
bellwethers_management = cluster_temporal_variability(context_layers.management_area, context_layers.temporal_scores, 0.1)
bellwethers_bioregion = cluster_temporal_variability(context_layers.bioregion, context_layers.temporal_scores, 0.1)

consistent_bellwether_reefs = bellwethers_gbr .& bellwethers_management .& bellwethers_bioregion
context_layers.consistent_bellwethers = consistent_bellwether_reefs
context_layers.consistent_bellwethers_cat = ifelse.(context_layers.consistent_bellwethers, "bellwether", "non-bellwether")


consistent_reefs_dhw = extract_timeseries(dhw_ts, context_layers, ["bioregion", "$(GCM)_lag4_bior", "consistent_bellwethers_cat", "management_area", "gbr"])
consistent_reefs_rel_cover = extract_timeseries(rel_cover, context_layers, ["bioregion", "$(GCM)_lag4_bior", "consistent_bellwethers_cat", "management_area", "gbr"])
consistent_reefs_total_cover = extract_timeseries(cover_ts, context_layers, ["bioregion", "$(GCM)_lag4_bior", "consistent_bellwethers_cat", "management_area", "gbr"])

consistent_reefs_dhw_plot = grouped_timeseries_plots(
    consistent_reefs_total_cover,
    "consistent_bellwethers_cat",
    "$(GCM)_lag4_bior",
    :gbr,
    (1,50), 0;
    ylabel="Degree Heating Weeks (°C)",
    xlabel="Year",
    x_fig_size = 2000,
    y_fig_size = 700
)
