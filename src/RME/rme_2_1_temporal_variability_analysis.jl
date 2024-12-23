include("../common.jl")

# Calculate coral cover relative to initial starting cover
total_cover = rs.total_relative_cover_median
rel_cover = Float64.(mapslices(relative_site_cover, total_cover, dims=[:timesteps]))

# Filter out reefs with less than 5% coral cover in initial timesteps
icc = total_cover[1,:].data
initial_cover_indices = icc .> 5

# Calculate the split timeseries and temporal variability metrics for the relative coral cover timeseries
temporal_scores = [temporal_variability(naive_split_metric(rel_cover[:, i].data, 10, mean)) for i in 1:length(rel_cover.sites)]
hist(temporal_scores[initial_cover_indices])

# Plot lines and the distribution of temporal variability scores
rel_cover_vec = [vec(rel_cover.data[:, i]) for i in 1:length(rel_cover.sites)]
sub_2_5 = initial_cover_indices .& (temporal_scores .< 2.5)

line_scores_zip = zip(rel_cover_vec[sub_2_5], temporal_scores[sub_2_5])
extremes = extrema(temporal_scores[sub_2_5])
lines_color((a, c)) = lines!(a; color=c, colorrange=extremes, alpha = 0.3)

fig = Figure()
Axis(fig[1,2])
map(lines_color, line_scores_zip)
Colorbar(fig[1,2], limits = extremes)
