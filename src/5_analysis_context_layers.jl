"""
Access the domain connectivity and calculate the strength and number of outgoing/incoming
connections for each reef. Run immediately after 4_.jl to use the updated context_layers with
bioregions.
"""

using DataFrames, Statistics, YAXArrays

using GLMakie, GeoMakie, GraphMakie

using ADRIA

import GeoDataFrames as GDF

import GeoFormatTypes as GFT

include("common.jl")

# GBR wide domain
gbr_domain_path = "c:/Users/bgrier/Documents/Projects/ADRIA_Domains/rme_ml_2024_01_08/"
gbr_dom = ADRIA.load_domain(RMEDomain, gbr_domain_path, "45")

# 1. Attach connectivity data to context_layers
connectivity_matrix = gbr_dom.conn

# Calculate connectivity metrics for context analysis
reefs = collect(getAxis("Source", connectivity_matrix).val)
source_to_sink = DataFrame(
    RME_GBRMPA_ID = [],
    income_strength = [],
    income_count = [],
    income_comb = [],
    out_strength = [],
    out_count = [],
    out_comb = [],
    total_strength = [],
    total_count = [],
    total_comb = [],
    so_to_si = []
)

for reef in eachindex(reefs)
    reef_id = reefs[reef]
    outgoing = connectivity_matrix[connectivity_matrix.Source .== reef_id, :]
    incoming = connectivity_matrix[:, connectivity_matrix.Sink .== reef_id]

    income_strength = sum(incoming)
    income_count = count(incoming .> 0)
    income_comb = income_strength / income_count

    out_strength = sum(outgoing)
    out_count = count(outgoing .> 0)
    out_comb = out_strength / out_count

    total_strength = income_strength + out_strength
    total_count = income_count + out_count
    total_comb = total_strength / total_count

    so_to_si = out_comb / income_comb

    push!(source_to_sink,
    [
        reef_id,
        income_strength,
        income_count,
        income_comb,
        out_strength,
        out_count,
        out_comb,
        total_strength,
        total_count,
        total_comb,
        so_to_si
    ]
    )
end
context_layers = leftjoin(context_layers, source_to_sink; on=:RME_GBRMPA_ID, order=:left)

# 2. Attach initial coral cover data to context_layers
initial_cc = gbr_dom.init_coral_cover
total_icc = Float64.(mapslices(sum, initial_cc, dims=[:species]))
context_layers.initial_coral_cover = total_icc.data

# species =
# [
# "tabular Acropora_2_1", "tabular Acropora_2_2",
# "tabular Acropora_2_3", "tabular Acropora_2_4", "tabular Acropora_2_5", "tabular Acropora_2_6", "tabular Acropora_2_7",
# "corymbose Acropora_3_1", "corymbose Acropora_3_2", "corymbose Acropora_3_3", "corymbose Acropora_3_4",
# "corymbose Acropora_3_5", "corymbose Acropora_3_6", "corymbose Acropora_3_7", "corymbose non-Acropora_4_1", "corymbose non-Acropora_4_2",
# "corymbose non-Acropora_4_3", "corymbose non-Acropora_4_4", "corymbose non-Acropora_4_5", "corymbose non-Acropora_4_6", "corymbose non-Acropora_4_7",
# "Small massives_5_1", "Small massives_5_2", "Small massives_5_3", "Small massives_5_4", "Small massives_5_5",
# "Small massives_5_6", "Small massives_5_7", "Large massives_6_1", "Large massives_6_2", "Large massives_6_3",
# "Large massives_6_4", "Large massives_6_5", "Large massives_6_6", "Large massives_6_7"
# ]
# spec_groups = getindex.(species, [1:15])

# icc = DataFrame(Matrix(initial_cc.data), collect(getAxis("locs", initial_cc).val))
# icc.spec_groups = spec_groups
# gdf = DataFrames.groupby(icc, :spec_groups)
# mapcols(gdf, sum)


# 3. Attach DHW data to context_layers
# mean DHW data for sites
dhw_time_scens = gbr_dom.dhw_scens
dhw_time = Float64.(mapslices(median, dhw_time_scens, dims=[:scenarios]))
dhw_locs = Float64.(mapslices(mean, dhw_time, dims=[:timesteps]))
dhw_locs = DataFrame(mean_dhw = dhw_locs.data, RME_GBRMPA_ID = collect(getAxis("locs", dhw_locs).val))

context_layers = leftjoin(context_layers, dhw_locs, on=:RME_GBRMPA_ID, order=:left)

# DHW correlation for sites
#rel_dhw = Float64.(mapslices(relative_site_cover, dhw_time, dims=[:timesteps]))
dhw_cor = DataFrame(UNIQUE_ID = [], dhw_cor = [])
for (ind, reef) in enumerate(eachcol(rel_cover))

    reef_cover = rel_cover[:,ind]
    id = collect(getAxis("sites", rel_cover).val)[ind]
    reef_dhw = dhw_time[:,ind]
    reef_dhw.data .= -1 * reef_dhw.data

    correlation = cross_correlation(reef_cover, reef_dhw, [0], true)
    push!(dhw_cor, [id; correlation])

end
context_layers = leftjoin(context_layers, dhw_cor, on=:RME_UNIQUE_ID, order=:left)

# 4. Attach site data/area to context_layers
context_layers = leftjoin(context_layers, gbr_dom.site_data[:,[:RME_UNIQUE_ID, :area]]; on=:RME_UNIQUE_ID, order=:left)

# 5. Calculate the number of reefs in each bioregion
context_layers.bioregion_count .= 1
for reef in eachrow(context_layers)
    reef.bioregion_count = count(context_layers.bioregion .== [reef.bioregion])
end

GDF.write("../data/analysis_context_layers.gpkg", context_layers; crs=GFT.EPSG(7844))
