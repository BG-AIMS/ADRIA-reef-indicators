using ReefModEngine
using YAXArrays, NetCDF

cd("src")
include("common.jl")

# Initialize RME (may take a minute or two)
init_rme("D:/development/RME/rme_ml_2024_06_13_v1.0.33/rme_ml_2024_06_13")
output_dir = "../outputs/RME_Results/"

# Set to use four threads
set_option("thread_count", 8)
set_option("use_fixed_seed", 1)  # Turn on use of a fixed seed value

# Get list of reef ids as specified by ReefMod Engine
reef_id_list = reef_ids()

# Get reef areas from RME
reef_area_km² = reef_areas()

start_year = 2022
end_year = 2099
RCP_scen = "SSP 2.45"  # RCP/SSP scenario to use
n_reefs = length(reef_id_list)

gcms = ["EC_Earth3_Veg", "UKESM1_0_LL", "GFDL-ESM4", "CNRM_ESM2_1"]

reps = 50  # Used here to indicate total number of runs
total_reps = 300

for gcm in gcms
    name = "reef-indicators SSP2-4.5 $(gcm)" # Name to associate with this set of runs
    # Initialize result store
    result_store = result_store = ResultStore(start_year, end_year, n_reefs)

    reset_rme()  # Reset RME to clear any previous runs

    for j in 1:6 # Have to run reps in batches of 50, so in total run 4 batches of 50 = 200 reps
        @info "running batch $(j) of 4 of $(reps) reps using model $(gcm)"

        set_option("fixed_seed", j)  # Set the fixed seed value

        # Note: if the Julia runtime crashes, check that the specified data file location is correct
        @RME runCreate(name::Cstring, start_year::Cint, end_year::Cint, RCP_scen::Cstring, gcm::Cstring, reps::Cint)::Cint

        # Initialize RME runs as defined above
        run_init()

        # Run all years and all reps
        @RME runProcess()::Cint

        concat_results!(result_store, start_year, end_year, reps)
        @info "results of batch $(j) of model $(gcm) concattenated"
    end

    @info "Saving results for model $(gcm) to result store"
    ReefModEngine.save_result_store(joinpath(output_dir, "RME_SSP245_$(total_reps)reps_$(gcm)"), result_store)
end

for gcm in gcms
    result_set = open_dataset(joinpath(output_dir, "RME_SSP245_$(reps*2)reps_$(gcm)/results.nc"), driver=:netcdf)
    result_set = remove_duplicate_reps(result_set, start_year, end_year, canonical_reefs.RME_UNIQUE_ID, reps*2)
    savedataset(result_set, path=joinpath(output_dir, "RME_SSP245_$(reps*2)reps_$(gcm)/results_no_duplicates.nc"), driver=:netcdf, overwrite=true)
end

for gcm in gcms
    results_200 = open_dataset(joinpath(output_dir, "RME_SSP245_200reps_$(gcm)/results_no_duplicates.nc"), driver=:netcdf)
    results_100 = open_dataset(joinpath(output_dir, "RME_SSP245_100reps_$(gcm)/results_no_duplicates.nc"), driver=:netcdf)
    results_300 = concat_RME_netcdfs(results_200, results_100)
    savedataset(results_300, path=joinpath(output_dir, "RME_SSP245_300reps_$(gcm)/results_no_duplicates.nc"), driver=:netcdf, overwrite=true)
end
