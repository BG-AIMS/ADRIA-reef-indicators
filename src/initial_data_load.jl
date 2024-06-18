rs = ADRIA.load_results("outputs/ReefMod Engine__RCPs_45__2024-06-17_15_06_38_869")
include("functions.jl")
context_layers = find_latest_file("../canonical-reefs/output/")
context_layers = GDF.read(context_layers)

using DataFrames, Statistics

using CSV

using GLMakie, GeoMakie, GraphMakie

using ADRIA, DynamicCoralCoverModel

moore_domain_path = "c:/Users/bgrier/Documents/Projects/ADRIA_Domains/Moore_2024-02-14_v060_rc1/"
# Local scale domain
moore_dom = ADRIA.load_domain(moore_domain_path)

# GBR wide domain
gbr_domain_path = "c:/Users/bgrier/Documents/Projects/ADRIA_Domains/rme_ml_2024_01_08/"
gbr_dom = ADRIA.load_domain(RMEDomain, gbr_domain_path, "45")

# generate 1024 sample scenarios from counterfactual scenarios
scens = ADRIA.sample_cf(gbr_dom, 16)

# Run sampled scenarios for a given RCP
rs = ADRIA.run_scenarios(gbr_dom, scens, "45")

# ... or repeat scenario runs across multiple RCPs
#rs = ADRIA.run_scenarios(dom, scens, ["45", "60", "85"])

# Analysis
s_tac = ADRIA.metrics.scenario_total_cover(rs, sites=1)
ADRIA.viz.scenarios(rs, s_tac)

tac = ADRIA.metrics.total_absolute_cover(rs)
tac_site_series = ADRIA.metrics.loc_trajectory(median, tac)

df = DataFrame(tac_site_series.data, collect(getAxis("sites", tac_site_series).val))
CSV.write("outputs/tac_site_series_19_06_2024.csv", df)


n_clusters = 5
clusters = ADRIA.analysis.cluster_scenarios(tac_site_series, n_clusters)
tac_sites = ADRIA.metrics.per_loc(median, tac)
tsc_map_fig = ADRIA.viz.map(rs, tac_sites, clusters)







df.year = [i for i in 1:size(df,1)]
select!(df, :year, Not(:year))
df.year = string.(df.year)

data = permutedims(df, 1, "UNIQUE_ID")

#remove the
test_vec = rand(15)
function relative_site_cover(x)
    init = x[1]
    for (index, step) in enumerate(x)
        x[index] = x[index]/init
    end

    return x
end


test = mapcols(relative_site_cover, df[:,2:end])
test.year = [string(i) for i in 1:size(df,1)]
select!(test, :year, Not(:year))

data = permutedims(test, 1, "UNIQUE_ID")
data = leftjoin(data, context_layers[:, [:UNIQUE_ID, :management_area]], on=:UNIQUE_ID)
filtered = data[(data[:,:UNIQUE_ID] .∉ [["15084100104", "14134100104", "15028100104"]]),:]

f, ax = plot_lines(data, :management_area, 1, 79, "Year", "Proportion of reef cover Year1")
hlines!(1, color=:green)
save("figs/all_reefs_initial_man_areas.png", f)


f, ax = plot_lines(filtered, :management_area, 1, 79, "Year", "Proportion of reef Year1")
hlines!(1, color=:green)
save("figs/filtered_reefs_initial_man_areas.png", f)
