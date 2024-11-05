using GLMakie
import GeoDataFrames as GDF

no_crash_gbr
no_crash_matrix = rel_cover[:, collect(getAxis("sites", rel_cover).val) .∈ [no_crash_gbr.RME_UNIQUE_ID]]

filtered_bior = no_crash_gbr[no_crash_gbr.bioregion .∈ [unique(no_crash_gbr[no_crash_gbr.target_reefs_bior, :bioregion])], :]
filtered_bior_matrix = rel_cover[:, collect(getAxis("sites", rel_cover).val) .∈ [filtered_bior.RME_UNIQUE_ID]]

target_order = combine(groupby(target_order,:management_area), sdf -> sort(sdf,:bioregion))
#target_order = target_order[target_order.so_to_si .< 15, :]
target_order.target_reefs_bior = ifelse.(target_order.target_reefs_bior, "target", "non-target")

test = DataFrame(no_crash_matrix.data, collect(getAxis("sites", no_crash_matrix).val))

function waterfall_categorical(data_df, order_df, color_by, xlab, ylab, zlab, alpha=0.7)
    data_matrix = Matrix(DataFrames.select(data_df, order_df.RME_UNIQUE_ID))

    if size(unique(order_df[:, color_by]),1) == 2
        alph_palette = colorscheme_alpha(ColorSchemes.tableau_20, Base.alpha)[2:3]
        palette = ColorSchemes.tableau_20.colors[2:3]
    elseif size(unique(order_df[:, color_by]),1) <= 20
        alph_palette = colorscheme_alpha(ColorSchemes.tableau_20, alpha)
        palette = ColorSchemes.tableau_20.colors
    else
        alph_palette = colorscheme_alpha(ColorSchemes.flag_ec, alpha)
        palette = ColorSchemes.flag_ec.colors
    end

    color_indices = groupindices(DataFrames.groupby(order_df, color_by))
    color_names = unique(DataFrame(indices=color_indices, names=order_df[:, color_by]))

    unique_names = color_names.names
    legend_entries = []
    for name in eachrow(color_names)
        col = palette[name.indices]
        LE = PolyElement(; color=col)
        push!(legend_entries, [LE])
    end

    xi = 1:1:78

    fig = Figure()
    ax = Axis3(fig[1,1];
        xticks=0:10:size(data_matrix, 1), yticks=0:200:size(target_order, 1),
        elevation=0.4, azimuth=-4,
        #limits=((0,80), (0), (0,1)),
        #xgridvisible=false, ygridvisible=false, zgridvisible=false,
        (Symbol.([:x, :y, :z], "spinecolor_", [2 3]) .=> :transparent)...,
        xlabel=xlab, ylabel=ylab, zlabel=zlab,
        xreversed=true)

    for (ind, z) in enumerate(eachcol(data_matrix))
        band!(Point3.(xi, ind, 0), Point3.(xi, ind, z), color=alph_palette[color_indices][ind])
        lines!(xi, fill(ind, length(xi)), z, color=palette[color_indices][ind])
    end

    Legend(fig[2, 1], legend_entries, unique_names, nbanks=3, tellheight=true, tellwidth=false, orientation=:horizontal, labelsize=10)

    display(fig)

    return fig
end

function waterfall_all_continuous(data_df, order_df, color_by, continuous_y, yticks, xlab, ylab, zlab, alpha=0.7)
    data_matrix = Matrix(DataFrames.select(df, order_df.RME_UNIQUE_ID))

    if size(unique(order_df[:, color_by]),1) == 2
        alph_palette = colorscheme_alpha(ColorSchemes.tableau_20, Base.alpha)[2:3]
        palette = ColorSchemes.tableau_20.colors[2:3]
    elseif size(unique(order_df[:, color_by]),1) <= 20
        alph_palette = colorscheme_alpha(ColorSchemes.tableau_20, alpha)
        palette = ColorSchemes.tableau_20.colors
    else
        alph_palette = colorscheme_alpha(ColorSchemes.flag_ec, alpha)
        palette = ColorSchemes.flag_ec.colors
    end

    yvals = order_df[:, continuous_y]

    color_indices = groupindices(DataFrames.groupby(order_df, color_by))
    color_names = unique(DataFrame(indices=color_indices, names=order_df[:, color_by]))

    unique_names = color_names.names
    legend_entries = []
    for name in eachrow(color_names)
        col = palette[name.indices]
        LE = PolyElement(; color=col)
        push!(legend_entries, [LE])
    end

    fig = Figure()
    ax = Axis3(fig[1,1];
        xticks=0:10:size(data_matrix, 1), yticks=yticks,
        #limits=((0,80), (0), (0,1)),
        #xgridvisible=false, ygridvisible=false, zgridvisible=false,
        (Symbol.([:x, :y, :z], "spinecolor_", [2 3]) .=> :transparent)...,
        xlabel=xlab, ylabel=ylab, zlabel=zlab,
        xreversed=true
    )

    for (ind, z) in enumerate(eachcol(data_matrix))
        band!(Point3.(xi, yvals[ind], 0), Point3.(xi, yvals[ind], z), color=alph_palette[color_indices][ind])
        lines!(xi, fill(yvals[ind], length(xi)), z, color=palette[color_indices][ind])
    end

    Legend(fig[2, 1], legend_entries, unique_names, nbanks=3, tellheight=true, tellwidth=false, orientation=:horizontal, labelsize=10)

    display(fig)
end

# # test lines
# df = DataFrame(tac_sites.data, collect(getAxis("sites", tac_sites).val))

# rel_cover = mapcols(relative_site_cover, df[2:79,:])
# rel_cover.year = [string(i+1) for i in 1:size(rel_cover,1)]
# select!(rel_cover, :year, Not(:year))

# data = permutedims(rel_cover, 1, "UNIQUE_ID")

# data = leftjoin(data, context_layers[:, [:UNIQUE_ID, :target_reefs]], on=:UNIQUE_ID)
# filtered = data[(data.UNIQUE_ID .∈ [filtered_bior.UNIQUE_ID]),:]

# data = leftjoin(data, context_layers[:, [:UNIQUE_ID, :management_area]], on=:UNIQUE_ID)

# data.target_reefs = ifelse.(data.target_reefs, "target", "non_target")

region_colors = Dict(
    "Far Northern Management Area" => colorscheme_alpha(ColorSchemes.Greens_9, 0.4),
    "Cairns/Cooktown Management Area" => colorscheme_alpha(ColorSchemes.lapaz10, 0.4),
    "Townsville/Whitsunday Management Area" => colorscheme_alpha(ColorSchemes.buda25, 0.4),
    "Mackay/Capricorn Management Area" => colorscheme_alpha(ColorSchemes.bilbao25, 0.4),
    "NA" => colorscheme_alpha(ColorSchemes.fes10, 0.4)
)

function waterfall_categorical(data_df, order_df, color_by, xlab, ylab, zlab, alpha=0.7)
    data_df = DataFrames.select(data_df, order_df.RME_UNIQUE_ID)

    fig = Figure()
    ax = Axis3(fig[1,1];
        xticks=0:10:78, yticks=0:200:size(target_order, 1),
        elevation=0.4, azimuth=-4,
        #limits=((0,80), (0), (0,1)),
        #xgridvisible=false, ygridvisible=false, zgridvisible=false,
        (Symbol.([:x, :y, :z], "spinecolor_", [2 3]) .=> :transparent)...,
        xlabel=xlab, ylabel=ylab, zlabel=zlab,
        xreversed=true
    )

    # if size(unique(order_df[:, color_by]),1) == 2
    #     alph_palette = colorscheme_alpha(ColorSchemes.tableau_20, Base.alpha)[2:3]
    #     palette = ColorSchemes.tableau_20.colors[2:3]
    # elseif size(unique(order_df[:, color_by]),1) <= 20
    #     alph_palette = colorscheme_alpha(ColorSchemes.tableau_20, alpha)
    #     palette = ColorSchemes.tableau_20.colors
    # else
    #     alph_palette = colorscheme_alpha(ColorSchemes.flag_ec, alpha)
    #     palette = ColorSchemes.flag_ec.colors
    # end

    xi = 1:1:78

# color_names = unique(DataFrame(indices=color_indices, names=order_df[:, color_by]))

    # unique_names = color_names.names
    # legend_entries = []
    # for name in eachrow(color_names)
    #     col = palette[name.indices]
    #     LE = PolyElement(; color=col)
    #     push!(legend_entries, [LE])
    # end
    unique_names = []
    legend_entries = []

    for region in unique(order_df.management_area)
        region_df = order_df[order_df.management_area .== region, :]

        color_indices = groupindices(DataFrames.groupby(region_df, color_by))
        alph_palette = region_colors[region]
        color_names = unique(DataFrame(indices=color_indices, names=region_df[:, color_by]))
        push!(unique_names, color_names.names)
        for name in eachrow(color_names)
            col=alph_palette[name.indices]
            LE = PolyElement(; color=col)
            push!(legend_entries, [LE])
        end

        data_matrix = Matrix(DataFrames.select(data_df, region_df.RME_UNIQUE_ID))

        for (ind_in_region, z) in enumerate(eachcol(data_matrix))
            ind = findall(x -> x==z, eachcol(data_df))
            band!(Point3.(xi, ind[1], 0), Point3.(xi, ind[1], z), color=alph_palette[color_indices][ind_in_region])
            lines!(xi, fill(ind[1], length(xi)), z, color=alph_palette[color_indices][ind_in_region])
        end
    end

    unique_names = vcat(unique_names...)
    Legend(fig[2, 1], legend_entries, unique_names, nbanks=3, tellheight=true, tellwidth=false, orientation=:horizontal, labelsize=5)

    display(fig)

    return fig
end
