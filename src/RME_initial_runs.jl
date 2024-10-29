using ReefModEngine, ADRIA
using YAXArrays, NetCDF
using GLMakie

include("common.jl")

# Initialize RME (may take a minute or two)
init_rme("C:/Users/bgrier/Documents/Projects/RME/rme_ml_2024_06_13")
# [ Info: Loaded RME 1.0.28

# Set to use four threads
set_option("thread_count", 6)
set_option("use_fixed_seed", 1)  # Turn on use of a fixed seed value
set_option("fixed_seed", 123.0)  # Set the fixed seed value

# Get list of reef ids as specified by ReefMod Engine
reef_id_list = reef_ids()

# Get reef areas from RME
reef_area_km² = reef_areas()

name = "reef-indicators SSP2-4.5 CNRM_ESM2_1"       # Name to associate with this set of runs
start_year = 2022
end_year = 2099
RCP_scen = "SSP 2.45"  # RCP/SSP scenario to use
gcm = "CNRM_ESM2_1"    # The base Global Climate Model (GCM)
n_reefs = length(reef_id_list)

reps = 10  # Used here to indicate total number of runs

# Initialize result store
result_store = result_store = ResultStore(start_year, end_year, n_reefs)

reset_rme()  # Reset RME to clear any previous runs

# Note: if the Julia runtime crashes, check that the specified data file location is correct
@RME runCreate(name::Cstring, start_year::Cint, end_year::Cint, RCP_scen::Cstring, gcm::Cstring, reps::Cint)::Cint

# Initialize RME runs as defined above
run_init()

# Run all years and all reps
@RME runProcess()::Cint

concat_results!(result_store, start_year, end_year, reps)

ReefModEngine.save_result_store(result_store, "../outputs/RME_result_stores/RME_SSP245_20reps")

prepare_ReefMod_results(
    "../outputs/RME_result_stores/RME_SSP245_20reps/",
    canonical_reefs.RME_UNIQUE_ID,
    start_year, end_year,
    "cover_and_evenness_$(Dates.format(today(), "Y_m_d")).nc",
    "scaled_evenness_additive",
    evenness_weight=0.6,
    cover_weight=0.4
)
