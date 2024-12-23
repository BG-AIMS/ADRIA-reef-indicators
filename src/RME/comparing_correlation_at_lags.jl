EC_earth3_Veg_covr_matrix = zeros(Float64, (length(context_layers.bioregion), 10))
for i in 1:10
    EC_earth3_Veg_covr_matrix[:, i] = cluster_correlation(context_layers.bioregion, taxa_evenness, i, CF)
end

initial_lag = cluster_correlation(context_layers.bioregion, taxa_evenness, 0, CF)
diff_matrix = copy(EC_earth3_Veg_covr_matrix)
for i in 1:10
    diff_matrix[:, i] = cluster_correlation(context_layers.bioregion, taxa_evenness, i, CF) .- initial_lag
end

fig = Figure()
ax1, hm = heatmap(fig[1,1], EC_earth3_Veg_covr_matrix)
Colorbar(fig[1,2], hm)

diff_matrix = ifelse.(diff_matrix .> 0, missing, diff_matrix)

ax2, hm = heatmap(fig[2,1], (diff_matrix))
Colorbar(fig[2,2], hm)

linkaxes!(ax1, ax2)
