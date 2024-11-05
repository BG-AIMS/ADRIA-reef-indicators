"""
plots of the number of reefs that are designated as leading at the different lag points.
This is using rme_2_subregion_correlation.jl and RME v.1.0.33.
"""

include("plotting_functions.jl")

# Plotting the number of bellwether reefs in each lag and each analysis level.
lags = 1:10
lags = ["lag$(a)" for a in eachindex(lags)]
lagged_numbers = DataFrame(
    lags = lags,
    mgmt_leading_reefs = zeros(Int64,10),
    bior_leading_reefs = zeros(Int64,10),
    subr_leading_reefs = zeros(Int64,10),
    mgmt_leading_reefs_evenness = zeros(Int64,10),
    bior_leading_reefs_evenness = zeros(Int64,10),
    subr_leading_reefs_evenness = zeros(Int64,10)
)
for row in eachrow(lagged_numbers)
    lag = row.lags
    row.mgmt_leading_reefs = sum(lagged_analysis_mgmt[:,lag] .>= 0.6)
    row.bior_leading_reefs = sum(lagged_analysis_bior[:, lag] .>= 0.6)
    row.subr_leading_reefs = sum(lagged_analysis_subregion[:, lag] .>= 0.6)

    row.mgmt_leading_reefs_evenness = sum(lagged_analysis_mgmt_evenness[:,lag] .>= 0.6)
    row.bior_leading_reefs_evenness = sum(lagged_analysis_bior_evenness[:, lag] .>= 0.6)
    row.subr_leading_reefs_evenness = sum(lagged_analysis_subregion_evenness[:, lag] .>= 0.6)
end

lagged_numbers_long = stack(lagged_numbers, 2:7)
lagged_numbers_long.lags = vcat(fill(1:10, 6)...)
lagged_numbers_long = lagged_numbers_long[lagged_numbers_long.lags .> 1 , :]
lagged_numbers_long_7plus = lagged_numbers_long[lagged_numbers_long.lags .> 6, :]
lagged_numbers_long_7plus.value_log = log10.((lagged_numbers_long_7plus.value .+ 1))



# Plotting each reef at lags to ensure there are no missing lags that would be incorrect or hazardous.
lags = 1:10
lags = ["lag$(a)" for a in eachindex(lags)]

analysis_dfs = (
    mgmt_leading_reefs = lagged_analysis_mgmt,
    bior_leading_reefs = lagged_analysis_bior,
    subr_leading_reefs = lagged_analysis_subregion,
    mgmt_leading_reefs_evenness = lagged_analysis_mgmt_evenness,
    bior_leading_reefs_evenness = lagged_analysis_bior_evenness,
    subr_leading_reefs_evenness = lagged_analysis_subregion_evenness
)

for analysis_level in keys(analysis_dfs)
    analysis_df = analysis_dfs[analysis_level]
    target_reefs_names = unique(vcat([analysis_df[analysis_df[:, "lag$(x)"] .> 0.8 , :RME_UNIQUE_ID] for x in 4:10]...))

    lagged_reef_names = DataFrame(
        reef_ids = target_reefs_names,
        lag4 = zeros(Int64, length(target_reefs_names)),
        lag5 = zeros(Int64, length(target_reefs_names)),
        lag6 = zeros(Int64, length(target_reefs_names)),
        lag7 = zeros(Int64, length(target_reefs_names)),
        lag8 = zeros(Int64, length(target_reefs_names)),
        lag9 = zeros(Int64, length(target_reefs_names)),
        lag10 = zeros(Int64, length(target_reefs_names))
    )

    for reef in eachrow(lagged_reef_names)
        for lag in ["lag$(x)" for x in 4:10]
            reef[lag] = Int64(reef.reef_ids ∈   String.(analysis_df[(analysis_df[:, lag] .> 0.8), :RME_UNIQUE_ID]))
        end
    end

    for col in names(lagged_reef_names)[2:end]
        if !any(Bool.(lagged_reef_names[:, col]))
            lagged_reef_names = DataFrames.select(lagged_reef_names, Not(col))
            continue
        end

        lagged_reef_names[Bool.(lagged_reef_names[:, col]), col] .= parse(Int64, split(col, "g")[2])
    end

    lags = names(lagged_reef_names)
    lags = lags[contains.(lags, "lag")]
    all_lags = vcat(["lag2", "lag3", lags]...)

    figure = radarplot_df(lagged_reef_names, lags, all_lags, ["no_label"]; labelsize = 0.5)
    save("../figs/lag_radarplot_$(String(analysis_level)).png", figure)
end
