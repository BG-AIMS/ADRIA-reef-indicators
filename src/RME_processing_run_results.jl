using ReefModEngine, ADRIA
using YAXArrays, NetCDF
using GLMakie

total_scens = 300
gcms = ["EC_Earth3_Veg", "UKESM1_0_LL", "GFDL-ESM4", "CNRM_ESM2_1"]
star_year = 2022
end_year = 2099

include("common.jl")

for gcm in gcms
    prepare_ReefMod_results(
        "../outputs/RME_result_stores/RME_SSP245_300reps/RME_SSP245_$(total_scens)reps_$(gcm)/",
        canonical_reefs.RME_UNIQUE_ID,
        start_year, end_year,
        "cover_and_evenness_$(Dates.format(today(), "Y_m_d")).nc",
        "scaled_evenness_additive",
        evenness_weight=1,
        cover_weight=1
    )
end

EC_Earth3_Veg_results = open_dataset("../outputs/RME_result_stores/RME_SSP245_300reps/RME_SSP245_$(total_scens)reps_EC_Earth3_Veg/cover_and_evenness_2024_11_13.nc"; driver=:netcdf)
UKESM1_0_LL_results = open_dataset("../outputs/RME_result_stores/RME_SSP245_300reps/RME_SSP245_$(total_scens)reps_UKESM1_0_LL/cover_and_evenness_2024_11_13.nc"; driver=:netcdf)
GFDL_ESM4_results = open_dataset("../outputs/RME_result_stores/RME_SSP245_300reps/RME_SSP245_$(total_scens)reps_GFDL-ESM4/cover_and_evenness_2024_11_13.nc"; driver=:netcdf)
CNRM_ESM2_1_results = open_dataset("../outputs/RME_result_stores/RME_SSP245_300reps/RME_SSP245_$(total_scens)reps_CNRM_ESM2_1/cover_and_evenness_2024_11_13.nc"; driver=:netcdf)

GCM_axis = Dim{:GCM}(gcms)
combined_cover = concatenatecubes(
    [
    EC_Earth3_Veg_results.total_relative_cover_median,
    UKESM1_0_LL_results.total_relative_cover_median,
    GFDL_ESM4_results.total_relative_cover_median,
    CNRM_ESM2_1_results.total_relative_cover_median
    ],
    GCM_axis
)
median_cover = Float64.(mapslices(median, combined_cover, dims=[:GCM]))

combined_evenness = concatenatecubes(
    [
    EC_Earth3_Veg_results.scaled_taxa_evenness,
    UKESM1_0_LL_results.scaled_taxa_evenness,
    GFDL_ESM4_results.scaled_taxa_evenness,
    CNRM_ESM2_1_results.scaled_taxa_evenness
    ],
    GCM_axis
)
median_evenness = Float64.(mapslices(median, combined_evenness, dims=[:GCM]))

ds = Dataset(; Dict(:total_relative_cover => median_cover, :scaled_taxa_evenness => median_evenness)...)

savedataset(ds; path="../outputs/RME_result_stores/RME_SSP245_300reps/median_GCM_cover_and_evenness_$(Dates.format(today(), "Y_m_d")).nc", driver=:netcdf, overwrite=true)
