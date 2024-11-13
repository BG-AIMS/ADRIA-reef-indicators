function process_median_ReefMod_results(
    result_store_dir,
    location_ids,
    start_year,
    end_year,
    median_indices,
    evenness_method;
    evenness_weight=1,
    cover_weight=1
)
    results = open_dataset("$(result_store_dir)/results_no_duplicates.nc")
    total_cover = results.total_cover[:,:,median_indices]
    taxa_cover = results.total_taxa_cover[:,:,:,median_indices]

    bootstrap_median = x -> first(bootstrap(median, x, BasicSampling(20)).t0)

    #new
    cover_median = ADRIA.metrics.summarize(total_cover, [:scenarios], bootstrap_median)
    taxa_median = ADRIA.metrics.summarize(taxa_cover, [:scenarios], bootstrap_median)

    # old
    # cover_median = Float64.(mapslices(median, total_cover, dims=[:scenarios]))
    # taxa_median = Float64.(mapslices(median, taxa_cover, dims=[:scenarios]))
    #
    taxa_evenness = _coral_evenness(taxa_median; method=evenness_method, evenness_weight=evenness_weight, cover_weight=cover_weight)

    axlist = (
                Dim{:timesteps}(start_year:end_year),
                Dim{:sites}(location_ids)
    )
    cover_median = rebuild(cover_median, dims=axlist)
    taxa_evenness = rebuild(taxa_evenness, dims=axlist)
    arrays = Dict(
        :total_relative_cover_median => cover_median,
        :scaled_taxa_evenness => taxa_evenness
    )
    ds = Dataset(; arrays...)

    return ds
end
result_dir = "../outputs/RME_result_stores/RME_SSP245_200reps/RME_SSP245_200reps_UKESM1_0_LL/"

scenario_sample_ranges = [5, 10, 20, 40, 70, 80, 100, 150, 170, 180, 200]
scenario_numbers = DataFrame(
    scenarios = scenario_sample_ranges,
    bior_leading_reefs = zeros(Int64,11),
    bior_leading_reefs_evenness = zeros(Int64,11),
)

for index_range in scenario_sample_ranges
    if index_range <= 40
        median_indices = randsubseq(1:200, 0.3)[rand(1:index_range, index_range)]
    else
        median_indices = rand(1:200, index_range)
    end

    rs_data = process_median_ReefMod_results(result_dir, canonical_reefs.RME_UNIQUE_ID, start_year, end_year, median_indices, "scaled_evenness_additive")

    cover_ts = rs_data.total_relative_cover_median
    cover_ts = cover_ts[1:50, :]
    taxa_evenness = rs_data.scaled_taxa_evenness
    taxa_evenness = taxa_evenness[1:50, :]

    # Apply analysis to bioregion subregions - 31 subregions
    bioregions = unique(context_layers.bioregion)
    lagged_analysis_bior = subregion_analysis(bioregions, cover_ts, context_layers, :bioregion, 1:10)
    scenario_numbers[scenario_numbers.scenarios .== index_range, :bior_leading_reefs] .= sum(lagged_analysis_bior[:, :lag4] .> 0.7)

    lagged_analysis_bior_evenness = subregion_analysis(bioregions, taxa_evenness, context_layers, :bioregion, 1:10)
    scenario_numbers[scenario_numbers.scenarios .== index_range, :bior_leading_reefs_evenness] .= sum(lagged_analysis_bior_evenness[:, :lag4] .> 0.7)
end
scenario_numbers_long = stack(scenario_numbers, 2:3)

scatter(
    scenario_numbers_long[scenario_numbers_long.variable .== "bior_leading_reefs", :scenarios],
    scenario_numbers_long[scenario_numbers_long.variable .== "bior_leading_reefs", :value];
    color=(:green, 0.5)
)
scatter!(
    scenario_numbers_long[scenario_numbers_long.variable .== "bior_leading_reefs_evenness", :scenarios],
    scenario_numbers_long[scenario_numbers_long.variable .== "bior_leading_reefs_evenness", :value];
    color=(:blue, 0.5)
)
