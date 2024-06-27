"""
Access the domain connectivity and calculate the strength and number of outgoing/incoming
connections for each reef.
"""

using DataFrames, Statistics

using CSV, YAXArrays

using GLMakie, GeoMakie, GraphMakie

using Infiltrator

using ADRIA, CoralBlox

#moore_domain_path = "c:/Users/bgrier/Documents/Projects/ADRIA_Domains/Moore_2024-02-14_v060_rc1/"
# Local scale domain
#moore_dom = ADRIA.load_domain(moore_domain_path)

# GBR wide domain
gbr_domain_path = "c:/Users/bgrier/Documents/Projects/ADRIA_Domains/rme_ml_2024_01_08/"
gbr_dom = ADRIA.load_domain(RMEDomain, gbr_domain_path, "45")


connectivity_matrix = gbr_dom.conn

reefs = collect(getAxis("Source", connectivity_matrix).val)
source_to_sink = DataFrame(GBRMPA_ID = [], income_strength = [], income_count = [], income_comb = [], out_strength = [], out_count = [], out_comb = [], total_connect = [], total_comb = [])
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

    total_connect = income_count + out_count
    total_comb = out_comb / income_comb

    push!(source_to_sink, [reef_id, income_strength, income_count, income_comb, out_strength, out_count, out_comb, total_connect, total_comb])
end

source_to_sink_context = leftjoin(context_layers, source_to_sink; on=:GBRMPA_ID, order=:left)

hist(source_to_sink.out_in_strength; bins=[0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.5,4,5,6])
vlines!(1, color=:red)

hist(source_to_sink.out_in_count, bins=[0,0.5,1,1.5,2,2.5,3,3.5,4,5,7.5,10,12.5,15,20,25,30,35,40,50,60,70,80,90,150])
vlines!(1, color=:red)

hist(source_to_sink, bins=[0,0.5,1,1.5,2,2.5,3,3.5,4,5,7.5,10,12.5,15,20,25,30,35,40,50,60,70,80,90,150])

target_reefs_context = context_layers[context_layers.UNIQUE_ID .∈ [target_reefs], :]
target_reefs_context = leftjoin(target_reefs_context, source_to_sink, on=:GBRMPA_ID, order=:left)


hist(source_to_sink.total_comb; bins=[0,0.5,1,1.5,2,2.5,3,3.5,4,5,7.5,10,12.5,15,20,25,30,35,40,50,60,70,80,90,150])
vlines!(1, color=:red)

show(target_reefs_context.total_comb)
hist(target_reefs_context.total_comb)
vlines!(1, color=:red)


hist(source_to_sink.total_connect)
show(target_reefs_context.total_connect)


hist(source_to_sink.income_strength)
hist(target_reefs_context.income_strength)

hist(source_to_sink.income_count)
hist(target_reefs_context.income_count)
