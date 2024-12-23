"""
Processing of runs that were created on the remote desktop and have been organised by their
GCM.
"""

using ReefModEngine, ADRIA
using YAXArrays, NetCDF
using GLMakie

total_scens = 300
gcms = ["EC_Earth3_Veg", "UKESM1_0_LL", "GFDL-ESM4", "CNRM_ESM2_1"]
start_year = 2022
end_year = 2099

include("../common.jl")

EC_Earth3_Veg_results = open_dataset("../outputs/RME_result_stores/RME_SSP245_300reps/RME_SSP245_$(total_scens)reps_EC_Earth3_Veg/results_no_duplicates.nc"; driver=:netcdf)
UKESM1_0_LL_results = open_dataset("../outputs/RME_result_stores/RME_SSP245_300reps/RME_SSP245_$(total_scens)reps_UKESM1_0_LL/results_no_duplicates.nc"; driver=:netcdf)
GFDL_ESM4_results = open_dataset("../outputs/RME_result_stores/RME_SSP245_300reps/RME_SSP245_$(total_scens)reps_GFDL-ESM4/results_no_duplicates.nc"; driver=:netcdf)
CNRM_ESM2_1_results = open_dataset("../outputs/RME_result_stores/RME_SSP245_300reps/RME_SSP245_$(total_scens)reps_CNRM_ESM2_1/results_no_duplicates.nc"; driver=:netcdf)

# Create a combined dataset with the median over GCMs
# combined_results = concat_RME_netcdfs(EC_Earth3_Veg_results, UKESM1_0_LL_results, GFDL_ESM4_results, CNRM_ESM2_1_results)
# savedataset(combined_results, path="../outputs/RME_result_stores/RME_SSP245_300reps/results_no_duplicates.nc", driver=:netcdf)
combined_results = open_dataset("../outputs/RME_result_stores/RME_SSP245_300reps/results_no_duplicates.nc"; driver=:netcdf)
prepare_ReefMod_results(
    combined_results,
    "../outputs/RME_result_stores/RME_SSP245_300reps/",
    canonical_reefs.RME_UNIQUE_ID,
    start_year, end_year,
    "cover_and_evenness_$(Dates.format(today(), "Y_m_d")).nc",
    "scaled_evenness_additive",
    evenness_weight=1,
    cover_weight=1
)

# Create processed datasets for each GCM separately.
for GCM in gcms
    GCM_dir = "../outputs/RME_result_stores/RME_SSP245_300reps/RME_SSP245_300reps_$(GCM)"
    results = open_dataset(joinpath(GCM_dir, "results_no_duplicates.nc"); driver=:netcdf)
    prepare_ReefMod_results(
        results,
        GCM_dir,
        canonical_reefs.RME_UNIQUE_ID,
        start_year, end_year,
        "cover_and_evenness_$(Dates.format(today(), "Y_m_d")).nc",
        "scaled_evenness_additive",
        evenness_weight=1,
        cover_weight=1
    )
end
