using RCall

url = "s3://gbr-dms-data-public/aims-ltmp-mmp-coralreef-model/data.parquet"
df = R"arrow::s3_bucket($(url))"
df = R"arrow::open_dataset($(df))"
variables = ["HARD CORAL", "SOFT CORAL"]
R"library(dplyr)"

LTMP = R"collect($(df))"
test = rcopy(LTMP)
test = select(test, Not(:geometry))
Parquet.write_parquet("../data/LTMP_all_data.parquet", test)



test = Parquet.read_parquet("../data/LTMP_all_data.parquet")
test = DataFrame(test)
test = test[test.variable .∈ [["HARD CORAL", "SOFT CORAL"]], :]

summed_coral_cover = combine(groupby(test, [:domain_name, :date]), :mean => sum)

LTMP_data = GDF.read("../../LTMP_Formatting/output/ltmp_reef_level_mean_data.gpkg")
manta_tow = GDF.read("../../LTMP_Formatting/output/manta_tow_data_reef_lvl.gpkg")

reefs_with_sites = ((context_layers.RME_UNIQUE_ID .∈ [unique(LTMP_data.RME_UNIQUE_ID)]) .| (context_layers.RME_UNIQUE_ID .∈ [unique(manta_tow.RME_UNIQUE_ID)]))
reefs_with_sites = ifelse.(ismissing.(reefs_with_sites), 0, reefs_with_sites)
context_layers.LTMP_photo_or_manta = reefs_with_sites

ltmp_array = YAXArray(
    (
        Dim{:locations}(unique(LTMP_data.RME_UNIQUE_ID)),
        Dim{:years}(unique(LTMP_data.YEAR_CODE))
    ),
    Matrix{Union{Missing, Float64}}(missing, length(unique(LTMP_data.RME_UNIQUE_ID)), length(unique(LTMP_data.YEAR_CODE)))
)

for (i, reef) in enumerate(ltmp_array.locations)
    for (j, year) in enumerate(ltmp_array.years)
        dat = LTMP_data[(LTMP_data.RME_UNIQUE_ID .== reef) .& (LTMP_data.YEAR_CODE .== year) .& (LTMP_data.GROUP_CODE .== "Hard Coral"), :]
        if size(dat, 1) == 1
            ltmp_array[i, j] = first(dat.COVER_mean)
        end
    end
end

ltmp_bellwether = ltmp_array.locations .∈ [filtered_bior_evenness[filtered_bior_evenness.target_reefs_bior_evenness, :RME_UNIQUE_ID]]
sum(ltmp_bellwether)
