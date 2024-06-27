"""
Load the initially created result set and calculate the metrics needed for later analyses.
(Currently these metrics are recalculated in later scripts.)
"""
using CSV

using GLMakie, GeoMakie, GraphMakie

using YAXArrays, DimensionalData

using DataFrames, Statistics

using ADRIA, CoralBlox

include("src/common.jl")

rs = ADRIA.load_results("outputs/ADRIA-out/ReefMod Engine__RCPs_45__2024-06-19_10_58_55_647")

context_layers = find_latest_file("../canonical-reefs/output/")
context_layers = GDF.read(context_layers)

# Analysis
s_tac = ADRIA.metrics.scenario_total_cover(rs)
ADRIA.viz.scenarios(rs, s_tac)

tac = ADRIA.metrics.total_absolute_cover(rs)
#tac_site_series = ADRIA.metrics.loc_trajectory(median, tac)

tac_sites = mapslices(median, tac, dims=[:scenarios])

#df = DataFrame(tac_sites.data, collect(getAxis("sites", tac_sites).val))
#CSV.write("outputs/tac_site_series_19_06_2024_depth_new.csv", df)
df = CSV.read("outputs/tac_site_series_19_06_2024_depth_new.csv", DataFrame)

rel_cover = mapcols(relative_site_cover, df[2:79,:])
rel_cover.year = [string(i+1) for i in 1:size(rel_cover,1)]
select!(rel_cover, :year, Not(:year))

data = permutedims(rel_cover, 1, "UNIQUE_ID")

data = leftjoin(data, context_layers[:, [:UNIQUE_ID, :management_area]], on=:UNIQUE_ID)
filtered = data[(data[:,:UNIQUE_ID] .∉ [["15084100104", "14134100104", "15028100104"]]),:]

data = leftjoin(data, context_layers[:, [:UNIQUE_ID, :closest_port]], on=:UNIQUE_ID)

f, ax = plot_lines(data, :management_area, 2, 79, "Year", "Proportion of reef cover Year2", 0.2)
hlines!(1, color=:green)
save("figs/all_reefs_initial_man_areas.png", f)


f, ax = plot_lines(data, :management_area, 2, 17, "Year", "Proportion of reef cover Year2", 0.1)
hlines!(1, color=:green)
save("figs/initial_man_areas_15_years.png", f)
