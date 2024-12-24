"""
Access the domain connectivity and calculate the strength and number of outgoing/incoming
connections for each reef. Run immediately after 4_.jl to use the updated context_layers with
bioregions.
"""

include("../common.jl")
include("../plotting_functions.jl")
include("custom_RMEDomain_dataloading.jl")

using DataFrames, Statistics, YAXArrays

using
    #GLMakie,
    GeoMakie,
    GraphMakie,
    CairoMakie

using GLM, MixedModels

# using ADRIA

import GeoDataFrames as GDF

CairoMakie.activate!(type = "png")


# Load updated context layer data
context_layers = GDF.read("../data/analysis_context_layers_rme_CF.gpkg",)
context_layers = context_layers[context_layers.bioregion .!= "NA", :]

gcms = ["EC_Earth3_Veg", "UKESM1_0_LL", "GFDL-ESM4", "CNRM_ESM2_1"]


for GCM in gcms

    # Coral Cover Analysis
    # Load rel_cover data and identify reefs that crash in the first timeseries year

    rs = open_dataset("../outputs/RME_result_stores/RME_SSP245_300reps/RME_SSP245_300reps_$(GCM)/cover_and_evenness_2024_12_18.nc"; driver=:netcdf)
    rel_cover = rs.total_relative_cover_median
    #crash_reefs = collect(getAxis("sites", rel_cover).val)[rel_cover.data[2,:] .< 0.5]

    dhw_data = load_DHW("../../RME/rme_ml_2024_06_13/data_files/", "245", GCM)
    dhw_ts = Float64.(mapslices(median, dhw_data, dims=[:scenarios]))
    dhw_ts = dhw_ts[locs = At(context_layers.RME_GBRMPA_ID)]
    axlist = (
        Dim{:timesteps}(2022:2099),
        Dim{:sites}(context_layers.RME_UNIQUE_ID)
    )
    dhw_ts = rebuild(dhw_ts, dims=axlist)

    # Remove crashing reefs
    #no_crash_gbr = context_layers[context_layers.UNIQUE_ID .∉ [crash_reefs],:]
    if sum(context_layers[:, "$(GCM)_target_reefs_bior"]) > 1
        # Filter to only include reefs that are within the same bioregion/closest_port subregion as target reefs
        filtered_bior = context_layers[context_layers.management_area .∈ [unique(context_layers[context_layers[:, "$(GCM)_target_reefs_bior_cat"] .== "bellwether", :management_area])], :]
        filtered_bior = filtered_bior[.!(ismissing.(filtered_bior[:, "$(GCM)_dhw_species2_cover_cor"])), :]
        filtered_bior[!, "$(GCM)_dhw_species2_cover_cor"] = Float64.(filtered_bior[:, "$(GCM)_dhw_species2_cover_cor"])
        filtered_bior = filtered_bior[filtered_bior.so_to_si .< quantile(filtered_bior.so_to_si, 0.95), :]
        filtered_bior = filtered_bior[filtered_bior.conn_score .< quantile(filtered_bior.conn_score, 0.95), :]
        filtered_bior = filtered_bior[filtered_bior.bioregion .∈ [unique(filtered_bior[filtered_bior[:, "$(GCM)_target_reefs_bior_cat"] .== "bellwether", :management_area])], :]
        filtered_bior = filtered_bior[filtered_bior.bioregion .∈ [unique(filtered_bior[filtered_bior[:, "$(GCM)_target_reefs_bior_cat"] .== "bellwether", :management_area])], :]

        removed_bioregions = bioregion_counts(
            filtered_bior,
            "$(GCM)_target_reefs_bior",
            "../outputs/RME_result_stores/RME_SSP245_300reps/RME_SSP245_300reps_$(GCM)/cover_reef_counts.csv",
            GCM
        )
        filtered_bior = filtered_bior[filtered_bior.bioregion .∉ [removed_bioregions], :]

        @info "$(removed_bioregions) bioregions have been removed from $(GCM) cover analysis as they are lacking bellwether reef numbers"
        # Basic exploratory models of factors
        #glm_allreefs = glm(@formula(target_reefs_bior ~ out_comb + dhw_cor ), no_crash_gbr, Binomial())
        #glm_bioregions = glm(@formula(target_reefs_bior ~ mean_dhw + year_8DHW + so_to_si + conn_score + total_strength + initial_coral_cover + dhw_cover_cor + mean_cots_mortality), filtered_bior, Binomial())
        #glm_subregions = glm(@formula(target_reefs_subr ~ total_comb + so_to_si + initial_coral_cover + bioregion_count), filtered_subr, Binomial())

        # aic(glm_allreefs), aic(glm_bioregions), aic(glm_subregions)

        # y_term_cover = term("$(GCM)_target_reefs_bior")
        # mean_dhw_term = term("$(GCM)_mean_dhw")
        # icc_term = term("$(GCM)_initial_coral_cover")
        # #dhw_cover_cor_term = term("$(GCM)_dhw_cover_cor")
        # cots_mortality_term = term("$(GCM)_mean_cots_mortality")

        # glmm_form =
        #     y_term_cover ~
        #     mean_dhw_term +
        #     term(:so_to_si) +
        #     term(:conn_score) +
        #     term(:total_strength) +
        #     icc_term +
        #     #dhw_cover_cor_term +
        #     cots_mortality_term +
        #     (term(1)|term(:bioregion))
        # glmm_fit = fit(MixedModel, glmm_form, filtered_bior, Bernoulli())
        # coef_table = coeftable(glmm_fit)
        # results = DataFrame(
        #     Names=coef_table.rownms,
        #     Coef=coef_table.cols[1],
        #     std_error=coef_table.cols[2],
        #     z=coef_table.cols[3],
        #     p=coef_table.cols[4]
        # )
        # CSV.write("../outputs/RME_result_stores/RME_SSP245_300reps/RME_SSP245_300reps_$(GCM)/cover_model_results.csv", results)

        # aic(glmm_fit) # seems to be an improvement when (1|bioregion) is used

        filtered_bior_bioregions = length(unique(filtered_bior.bioregion))
        # Plotting bioregion violin plots for each variable
        mean_dhw_violin = grouped_violin_plots(
            filtered_bior,
            Symbol("$(GCM)_target_reefs_bior_cat"),
            :bioregion, Symbol("$(GCM)_mean_dhw");
            ylabel="mean DHW", xlabel="Bellwether Reefs"
        );
        save("../figs/$(GCM)/mean_dhw_violin_cover.png", mean_dhw_violin)

        so_to_si_violin = grouped_violin_plots(
            filtered_bior,
            Symbol("$(GCM)_target_reefs_bior_cat"),
            :bioregion, :so_to_si;
            ylabel="Source to sink ratio", xlabel="Bellwether Reefs"
        );
        save("../figs/$(GCM)/so_to_si_violin_cover.png", so_to_si_violin)

        conn_score_violin = grouped_violin_plots(
            filtered_bior,
            Symbol("$(GCM)_target_reefs_bior_cat"),
            :bioregion, :conn_score;
            ylabel="Connectivity eigenvector centrality", xlabel="Bellwether Reefs"
        );
        save("../figs/$(GCM)/conn_score_violin_cover.png", conn_score_violin)

        total_strength_violin = grouped_violin_plots(
            filtered_bior,
            Symbol("$(GCM)_target_reefs_bior_cat"),
            :bioregion, :total_strength;
            ylabel="Total connectivity strength", xlabel="Bellwether Reefs"
        );
        save("../figs/$(GCM)/total_strength_violin_cover.png", total_strength_violin)

        initial_cover_violin = grouped_violin_plots(
            filtered_bior,
            Symbol("$(GCM)_target_reefs_bior_cat"),
            :bioregion, Symbol("$(GCM)_initial_coral_cover");
            ylabel="Initial coral cover", xlabel="Bellwether Reefs"
        );
        save("../figs/$(GCM)/initial_cover_violin_cover.png", initial_cover_violin)

        dhw_cover_cor_violin = grouped_violin_plots(
            filtered_bior,
            Symbol("$(GCM)_target_reefs_bior_cat"),
            :bioregion, Symbol("$(GCM)_dhw_cover_cor");
            ylabel="Total coral cover - DHW correlation", xlabel="Bellwether Reefs"
        );
        save("../figs/$(GCM)/dhw_cover_cor_violin.png", dhw_cover_cor_violin)

        for fg in 1:6

            fg_dhw_cover_cor_violin = grouped_violin_plots(
                filtered_bior,
                Symbol("$(GCM)_target_reefs_bior_cat"),
                :bioregion, Symbol("$(GCM)_dhw_species$(fg)_cover_cor");
                ylabel="Functional Group-$(fg) DHW - cover correlation", xlabel="Bellwether Reefs"
            );
            save("../figs/$(GCM)/dhw_fg$(fg)_cover_cor_violin_cover.png", fg_dhw_cover_cor_violin)
        end

        mean_COTS_mortality_violin = grouped_violin_plots(
            filtered_bior,
            Symbol("$(GCM)_target_reefs_bior_cat"),
            :bioregion, Symbol("$(GCM)_mean_cots_mortality");
            ylabel="mean CoTS mortality", xlabel="Bellwether Reefs"
        );
        save("../figs/$(GCM)/mean_COTS_mortality_violin_cover.png", mean_COTS_mortality_violin)

        rel_cover = rs.total_relative_cover_median
        rel_cover_reefs = extract_timeseries(rel_cover, context_layers, ["bioregion", "$(GCM)_lag4_bior", "$(GCM)_target_reefs_bior_cat", "management_area"])
        rel_cover_reefs = rel_cover_reefs[rel_cover_reefs.RME_UNIQUE_ID .∈ [filtered_bior.RME_UNIQUE_ID], :]

        # timeseries_bior_reefs = timeseries_plot(rel_cover_reefs, :target_reefs_bior_cat, (1,50), "Total relative coral cover (%)", 0.05)
        # save("../figs/bior_reefs_timeseries_plot.png", timeseries_bior_reefs)

        # lagged_ts_bior_reefs = lagged_timeseries_plot(rel_cover_reefs, :target_reefs_bior_cat, (1,50), "Year", "Total relative coral cover (%)", 0.1, 4)
        # save("../figs/bior_reefs_lagged_timeseries.png", lagged_ts_bior_reefs)

        # lagged_ts_combined_bior = lagged_timeseries_plot_combined(rel_cover_reefs, :target_reefs_bior_cat, (1,50), "Year", "Total relative coral cover (%)", 0.1, 4)
        # save("../figs/bior_reefs_combined_timeseries.png", lagged_ts_combined_bior)

        bioregion_grouped_timeseries = grouped_timeseries_plots(
            rel_cover_reefs,
            "$(GCM)_target_reefs_bior_cat",
            "$(GCM)_lag4_bior",
            :management_area,
            (1,50), 0;
            ylabel="Total coral cover",
            xlabel="Year"
        )
        save("../figs/$(GCM)/bioregion_cover_timeseries.png", bioregion_grouped_timeseries)
        bioregion_grouped_timeseries_lagged = grouped_timeseries_plots(
            rel_cover_reefs,
            "$(GCM)_target_reefs_bior_cat",
            "$(GCM)_lag4_bior",
            :bioregion,
            (1,50), 4;
            ylabel="Total coral cover",
            xlabel="Year"
        )
        save("../figs/$(GCM)/bioregion_cover_timeseries_lagged.png", bioregion_grouped_timeseries_lagged)

        bior_reefs_dhw = extract_timeseries(dhw_ts, filtered_bior, ["bioregion", "$(GCM)_lag4_bior", "$(GCM)_target_reefs_bior_cat", "management_area"])
        bior_reefs_dhw_plot = grouped_timeseries_plots(
            bior_reefs_dhw,
            "$(GCM)_target_reefs_bior_cat",
            "$(GCM)_lag4_bior",
            :management_area,
            (1,50), 0;
            ylabel="Degree Heating Weeks (°C)",
            xlabel="Year",
            x_fig_size = 2000,
            y_fig_size = 700
        )
        save("../figs/$(GCM)/cover_bellwether_dhw_timeseries.png", bior_reefs_dhw_plot)

        # CC_filtered = filtered_bior[filtered_bior.management_area_short .== "Cairns-Cooktown", :]
        # cover_CC_map = plot_map_bellwether_reefs(CC_filtered, :bioregion, :target_reefs_bior_cat)[1]
        # save("../figs/cover_CairnsCooktown_map.png", cover_CC_map)

        # TSV_filtered = filtered_bior[filtered_bior.management_area_short .== "Townsville-Whitsunday", :]
        # cover_TSV_map = plot_map_bellwether_reefs(TSV_filtered, :bioregion, :target_reefs_bior_cat)[1]
        # save("../figs/cover_TownsvilleWhitsunday_map.png", cover_TSV_map)

        # MCap_filtered = filtered_bior[filtered_bior.management_area_short .== "Mackay-Capricorn", :]
        # cover_MCap_map = plot_map_bellwether_reefs(MCap_filtered, :bioregion, :target_reefs_bior_cat)[1]
        # save("../figs/cover_MackayCapricorn_map.png", cover_MCap_map)
    end

    if sum(context_layers[:, "$(GCM)_target_reefs_bior_evenness"]) > 1
        # Evenness
        filtered_bior_evenness = context_layers[context_layers.bioregion .∈ [unique(context_layers[context_layers[:, "$(GCM)_target_reefs_bior_evenness_cat"] .== "bellwether", :bioregion])], :]
        filtered_bior_evenness = filtered_bior_evenness[.!(ismissing.(filtered_bior_evenness[:, "$(GCM)_dhw_species2_cover_cor"])), :]
        filtered_bior_evenness[!, "$(GCM)_dhw_species2_cover_cor"] = Float64.(filtered_bior_evenness[:, "$(GCM)_dhw_species2_cover_cor"])
        filtered_bior_evenness = filtered_bior_evenness[filtered_bior_evenness.so_to_si .< quantile(filtered_bior_evenness.so_to_si, 0.95), :]
        filtered_bior_evenness = filtered_bior_evenness[filtered_bior_evenness.conn_score .< quantile(filtered_bior_evenness.conn_score, 0.95), :]
        filtered_bior_evenness = filtered_bior_evenness[filtered_bior_evenness.bioregion .∈ [unique(filtered_bior_evenness[filtered_bior_evenness[:, "$(GCM)_target_reefs_bior_evenness_cat"] .== "bellwether", :bioregion])], :]

        evenness = rs.scaled_taxa_evenness

        # Remove reefs that have a qc flag indicating they have NaN/0 evenness values for every year
        filtered_bior_evenness_no_zeros = filtered_bior_evenness[.!(Bool.(filtered_bior_evenness[:, "$(GCM)_bior_cover_qc_flag"])) , :]
        filtered_bior_evenness_no_zeros[!, "$(GCM)_dhw_evenness_cor"] = Float64.(filtered_bior_evenness_no_zeros[:, "$(GCM)_dhw_evenness_cor"])
        filtered_bior_evenness_no_zeros = filtered_bior_evenness_no_zeros[filtered_bior_evenness_no_zeros.bioregion .∈ [unique(filtered_bior_evenness_no_zeros[filtered_bior_evenness_no_zeros[:, "$(GCM)_target_reefs_bior_evenness_cat"] .== "bellwether", :bioregion])], :]

        removed_bioregions_evenness = bioregion_counts(
            filtered_bior_evenness_no_zeros,
            "$(GCM)_target_reefs_bior_evenness",
            "../outputs/RME_result_stores/RME_SSP245_300reps/RME_SSP245_300reps_$(GCM)/evennss_reef_counts.csv",
            GCM
        )
        filtered_bior_evenness_no_zeros = filtered_bior_evenness_no_zeros[filtered_bior_evenness_no_zeros.bioregion .∉ [removed_bioregions_evenness], :]
        @info "$(removed_bioregions_evenness) bioregions have been removed from $(GCM) evenness analysis as they are lacking bellwether reef numbers"

        evenness_ts = extract_timeseries(evenness, filtered_bior_evenness_no_zeros, ["bioregion", "$(GCM)_lag4_bior_evenness", "$(GCM)_target_reefs_bior_evenness_cat", "management_area"])
        evenness_no_zeros = evenness_ts[indexin(filtered_bior_evenness_no_zeros.RME_UNIQUE_ID, evenness_ts.RME_UNIQUE_ID), :]

        # Basic exploratory models of factors

        #glm_bioregions = glm(@formula(target_reefs_bior_evenness ~ mean_dhw + so_to_si + conn_score + total_strength + initial_coral_cover + dhw_evenness_cor + mean_cots_mortality), filtered_bior_evenness_no_zeros, Binomial(), LogitLink())


        # y_term_evenness = term("$(GCM)_target_reefs_bior_evenness")
        # mean_dhw_term = term("$(GCM)_mean_dhw")
        # icc_term = term("$(GCM)_initial_coral_cover")
        # #dhw_cover_cor_term = term("$(GCM)_dhw_cover_cor")
        # cots_mortality_term = term("$(GCM)_mean_cots_mortality")

        # glmm_form =
        #     y_term_evenness ~
        #     mean_dhw_term +
        #     term(:so_to_si) +
        #     term(:conn_score) +
        #     term(:total_strength) +
        #     icc_term +
        #     #dhw_cover_cor_term +
        #     cots_mortality_term +
        #     (term(1)|term(:bioregion))
        # glmm_fit = fit(MixedModel, glmm_form, filtered_bior_evenness_no_zeros, Bernoulli())

        # #aic(glmm_fit), aic(glm_bioregions) # seems to be an improvement when (1|bioregion) is used
        # coef_table = coeftable(glmm_fit)
        # results = DataFrame(
        #     Names=coef_table.rownms,
        #     Coef=coef_table.cols[1],
        #     std_error=coef_table.cols[2],
        #     z=coef_table.cols[3],
        #     p=coef_table.cols[4]
        # )
        # CSV.write("../outputs/RME_result_stores/RME_SSP245_300reps/RME_SSP245_300reps_$(GCM)/evenness_model_results.csv", results)

        evenness_bioregions = length(unique(filtered_bior_evenness_no_zeros.bioregion))
        # Plotting bioregion violin plots for each variable
        mean_dhw_violin_evenness = grouped_violin_plots(
            filtered_bior_evenness_no_zeros,
            Symbol("$(GCM)_target_reefs_bior_evenness_cat"),
            :bioregion, Symbol("$(GCM)_mean_dhw");
            ylabel="mean DHW", xlabel="Bellwether Reefs"
        );
        save("../figs/$(GCM)/mean_dhw_violin_evenness.png", mean_dhw_violin_evenness)

        so_to_si_violin_evenness = grouped_violin_plots(
            filtered_bior_evenness_no_zeros,
            Symbol("$(GCM)_target_reefs_bior_evenness_cat"),
            :bioregion, :so_to_si;
            ylabel="Source to sink ratio", xlabel="Bellwether Reefs"
        );
        save("../figs/$(GCM)/so_to_si_violin_evenness.png", so_to_si_violin_evenness)

        conn_score_violin_evenness = grouped_violin_plots(
            filtered_bior_evenness_no_zeros,
            Symbol("$(GCM)_target_reefs_bior_evenness_cat"),
            :bioregion, :conn_score;
            ylabel="Connectivity eigenvector centrality", xlabel="Bellwether Reefs"
        );
        save("../figs/$(GCM)/conn_score_violin_evenness.png", conn_score_violin_evenness)

        total_strength_violin_evenness = grouped_violin_plots(
            filtered_bior_evenness_no_zeros,
            Symbol("$(GCM)_target_reefs_bior_evenness_cat"),
            :bioregion, :total_strength;
            ylabel="Total strength", xlabel="Bellwether Reefs"
        );
        save("../figs/$(GCM)/total_strength_violin_evenness.png", total_strength_violin_evenness)

        initial_cover_violin_evenness = grouped_violin_plots(
            filtered_bior_evenness_no_zeros,
            Symbol("$(GCM)_target_reefs_bior_evenness_cat"),
            :bioregion, Symbol("$(GCM)_initial_coral_cover");
            ylabel="Initial coral cover", xlabel="Bellwether Reefs"
        );
        save("../figs/$(GCM)/initial_cover_violin_evenness.png", initial_cover_violin_evenness)

        dhw_evenness_cor_violin = grouped_violin_plots(
            filtered_bior_evenness_no_zeros,
            Symbol("$(GCM)_target_reefs_bior_evenness_cat"),
            :bioregion, Symbol("$(GCM)_dhw_evenness_cor");
            ylabel="Taxa evenness - DHW correlation", xlabel="Bellwether Reefs"
        );
        save("../figs/$(GCM)/dhw_evenness_cor_violin_evenness.png", dhw_evenness_cor_violin)

        for fg in 1:6

            fg_dhw_cover_cor_violin = grouped_violin_plots(
                filtered_bior_evenness_no_zeros,
                Symbol("$(GCM)_target_reefs_bior_evenness_cat"),
                :bioregion, Symbol("$(GCM)_dhw_species$(fg)_cover_cor");
                ylabel="Functional Group-$(fg) DHW-cover correlation", xlabel="Bellwether Reefs"
            );
            save("../figs/$(GCM)/dhw_fg$(fg)_cover_cor_violin_evenness.png", fg_dhw_cover_cor_violin)
        end

        mean_COTS_mortality_violin_evenness = grouped_violin_plots(
            filtered_bior_evenness_no_zeros,
            Symbol("$(GCM)_target_reefs_bior_evenness_cat"),
            :bioregion, Symbol("$(GCM)_mean_cots_mortality");
            ylabel="mean CoTS mortality", xlabel="Bellwether Reefs"
        );
        save("../figs/$(GCM)/mean_COTS_mortality_violin_evenness.png", mean_COTS_mortality_violin_evenness)


        # timeseries_bior_evenness = timeseries_plot(evenness_no_zeros, :target_reefs_bior_evenness_cat, (1,50), "Scaled Taxa Evenness", 0.05)
        # save("../figs/$(GCM)/bior_evenness_timeseries_plot.png", timeseries_bior_evenness)

        # lagged_ts_bior_evenness = lagged_timeseries_plot(evenness_no_zeros, :target_reefs_bior_evenness_cat, (1,50), "Year", "Scaled Taxa Evenness", 0.05, 4)
        # save("../figs/$(GCM)/bior_evenness_lagged_timeseries.png", lagged_ts_bior_evenness)

        # combined_ts_bior_evenness = lagged_timeseries_plot_combined(evenness_no_zeros, :target_reefs_bior_evenness_cat, (1,50), "Year", "Scaled taxa evenness", 0.1, 4)
        # save("../figs/$(GCM)/bior_evenness_combined_timeseries.png", combined_ts_bior_evenness)

        bioregion_grouped_timeseries_evenness = grouped_timeseries_plots(
            evenness_no_zeros,
            "$(GCM)_target_reefs_bior_evenness_cat",
            "$(GCM)_lag4_bior_evenness",
            :bioregion,
            (1,50), 0;
            ylabel="Scaled taxa evenness",
            xlabel="Year"
        )
        save("../figs/$(GCM)/bioregion_evenness_timeseries.png", bioregion_grouped_timeseries_evenness)
        bioregion_grouped_timeseries_lagged_evenness = grouped_timeseries_plots(
            evenness_no_zeros,
            "$(GCM)_target_reefs_bior_evenness_cat",
            "$(GCM)_lag4_bior_evenness",
            :bioregion,
            (1,50), 4;
            ylabel="Scaled taxa evenness",
            xlabel="Year"
        )
        save("../figs/$(GCM)/bioregion_evenness_timeseries_lagged.png", bioregion_grouped_timeseries_lagged_evenness)

        bior_evenness_dhw = extract_timeseries(dhw_ts, filtered_bior_evenness_no_zeros, ["bioregion", "$(GCM)_lag4_bior_evenness","$(GCM)_target_reefs_bior_evenness_cat", "management_area"])
        bior_evenness_dhw_plot = grouped_timeseries_plots(
            bior_evenness_dhw,
            "$(GCM)_target_reefs_bior_evenness_cat",
            "$(GCM)_lag4_bior_evenness",
            :bioregion,
            (1,50), 0;
            ylabel="Degree Heating Week (°C)",
            xlabel="Year"
        )
        save("../figs/$(GCM)/evenness_bellwether_dhw_timeseries.png", bior_evenness_dhw_plot)

        # # Mean_DHW plot
        # mean_dhw_evenness_plot = explore_regression_scatter_plots(
        #     filtered_bior_evenness_no_zeros.mean_dhw,
        #     filtered_bior_evenness_no_zeros.lag4_bior_evenness .^2;
        #     color = categorical(filtered_bior_evenness_no_zeros.target_reefs_bior_evenness).refs,
        #     xlab="mean DHW",
        #     ylab="(correlation at lag 5) ^2",
        # )
        # save("../figs/$(GCM)/mean_dhw_evenness_correlation_scatter.png", mean_dhw_evenness_plot)

        # initial_cover_plot = explore_regression_scatter_plots(
        #     filtered_bior_evenness_no_zeros.initial_coral_cover,
        #     filtered_bior_evenness_no_zeros.lag4_bior_evenness
        # )

        # summarise_scores_evenness = combine(DataFrames.groupby(all_reefs_no_zeros, :bioregion), :lag4_bior_evenness => mean, :total_strength => mean, :ReefMod_area_m2 => mean, :conn_score => mean)
        # highest_bioregion = summarise_scores_evenness.bioregion[argmax(summarise_scores_evenness.lag4_bior_evenness_mean)]
        # lowest_bioregion = summarise_scores_evenness.bioregion[argmin(summarise_scores_evenness[summarise_scores_evenness.lag4_bior_evenness_mean .> 0.0, :lag4_bior_evenness_mean])]

        # FN_filtered_evenness = filtered_bior_evenness_no_zeros[filtered_bior_evenness_no_zeros.management_area_short .== "FarNorthern", :]
        # evenness_FN_map = plot_map_bellwether_reefs(FN_filtered_evenness, :bioregion, :target_reefs_bior_evenness_cat)[1]
        # save("../figs/$(GCM)/evenness_FarNorthern_map.png", evenness_FN_map)

        # CC_filtered_evenness = filtered_bior_evenness_no_zeros[filtered_bior_evenness_no_zeros.management_area_short .== "Cairns-Cooktown", :]
        # evenness_CC_map = plot_map_bellwether_reefs(CC_filtered_evenness, :bioregion, :target_reefs_bior_evenness_cat)[1]
        # save("../figs/$(GCM)/evenness_CairnsCooktown_map.png", evenness_CC_map)

        # TSV_filtered_evenness = filtered_bior_evenness_no_zeros[filtered_bior_evenness_no_zeros.management_area_short .== "Townsville-Whitsunday", :]
        # evenness_TSV_map = plot_map_bellwether_reefs(TSV_filtered_evenness, :bioregion, :target_reefs_bior_evenness_cat)[1]
        # save("../figs/$(GCM)/evenness_TownsvilleWhitsunday_map.png", evenness_TSV_map)

        # MCap_filtered_evenness = filtered_bior_evenness_no_zeros[filtered_bior_evenness_no_zeros.management_area_short .== "Mackay-Capricorn", :]
        # evenness_MCap_map = plot_map_bellwether_reefs(MCap_filtered_evenness, :bioregion, :target_reefs_bior_evenness_cat)[1]
        # save("../figs/$(GCM)/evenness_MackayCapricorn_map.png", evenness_MCap_map)
    end
end


# Creating ECS selection figure
ECS_values = DataFrame(GCM = ["EC-Earth3-Veg", "UKESM1-0-LL", "GFDL-ESM4", "CNRM-ESM2-1"], ECS = [4.33, 5.36, 3.9, 4.79])
likely = Point2f[(0.9, 2.5), (0.95, 2.5), (0.95, 4), (0.9, 4), (0.9, 2.5)]
very_likely = Point2f[(0.9, 2), (0.95, 2), (0.95, 5), (0.9, 5), (0.9, 2)]

low_impact = Point2f[(1.1, first(ECS_values[ECS_values.GCM .== "GFDL-ESM4", :ECS])), (1.11, first(ECS_values[ECS_values.GCM .== "GFDL-ESM4", :ECS]))]
mid_impact = Point2f[
    (1.1, first(ECS_values[ECS_values.GCM .== "EC-Earth3-Veg", :ECS])),
    (1.11, first(ECS_values[ECS_values.GCM .== "EC-Earth3-Veg", :ECS])),
    (1.11, first(ECS_values[ECS_values.GCM .== "CNRM-ESM2-1", :ECS])),
    (1.1, first(ECS_values[ECS_values.GCM .== "CNRM-ESM2-1", :ECS]))]
high_impact = Point2f[(1.1, first(ECS_values[ECS_values.GCM .== "UKESM1-0-LL", :ECS])), (1.11, first(ECS_values[ECS_values.GCM .== "UKESM1-0-LL", :ECS]))]

f = Figure(size = (800, 800))
ax = Axis(
    f[1,1];
    xticks = ([1.0],[""]),
    xticksvisible=false,
    limits= ((0.75, 1.25), (1,  6)),
    width=700,
    height=700,
    ylabel="Equilibrium Climate Sensitivity (ECS °C)",
    ylabelsize=16
)



scatter!(fill(0.925, length(ECS_values.ECS)), ECS_values.ECS; markersize = 15, color=:black)
text!(
    [1,1,1,1, 0.8, 0.8, 1.12, 1.12, 1.12],
    vcat(
        ECS_values.ECS .- 0.05,
        [
            3.25,
            2.2,
            first(ECS_values[ECS_values.GCM .== "GFDL-ESM4", :ECS]) - 0.05,
            mean(ECS_values[ECS_values.GCM .∈ [["EC-Earth3-Veg", "CNRM-ESM2-1"]], :ECS]) - 0.05,
            first(ECS_values[ECS_values.GCM .== "UKESM1-0-LL", :ECS]) - 0.05
        ]);
    text=vcat(ECS_values.GCM, ["likely\n(low confidence)", "very likely\n(high confidence)", "low impact", "mid impact", "high impact"]),
    fontsize=16
)
poly!(likely; color=(:red, 0.2), strokecolor=:black, strokewidth=2)
poly!(very_likely; color=(:red, 0.1), linestyle=(:dot, :loose), strokecolor=:black, strokewidth=2)
lines!.([low_impact, mid_impact, high_impact]; color=:black)

save("../figs/ECS_plot.png", f)
