using ReefModEngine, ADRIA
using YAXArrays, NetCDF
using GLMakie

include("common.jl")

prepare_ReefMod_results(
    "../outputs/RME_result_stores/RME_SSP245_100reps_$(gcm)/",
    canonical_reefs.RME_UNIQUE_ID,
    start_year, end_year,
    "cover_and_evenness_$(Dates.format(today(), "Y_m_d")).nc",
    "scaled_evenness_additive",
    evenness_weight=0.6,
    cover_weight=0.4
)
