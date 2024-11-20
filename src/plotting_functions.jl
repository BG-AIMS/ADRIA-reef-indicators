"""
File includes the plotting functions for ADRIA-reef-indicators analysis.
"""

function _convert_plottable(gdf::Union{DataFrame,DataFrameRow}, geom_col::Symbol)
    local plottable
    try
        if gdf isa DataFrame
            plottable = GeoMakie.geo2basic(AG.forceto.(gdf[!, geom_col], AG.wkbPolygon))
        else
            plottable = GeoMakie.geo2basic(AG.forceto(gdf[geom_col], AG.wkbPolygon))
        end
    catch
        # Column is already in a plottable form, or some unrelated error occurred
        if gdf isa DataFrame
            plottable = gdf[:, geom_col]
        else
            plottable = [gdf[geom_col]]
        end
    end

    return plottable
end

"""
    plot_map(gdf::DataFrame; geom_col::Symbol=:geometry, color_by::Symbol)

Convenience plot function.

# Arguments
- `gdf` : GeoDataFrame
- `color_by` : Column name holding factor to color reefs by (e.g. :management_area)
- `geom_col` : Column name holding geometries to plot
"""
function plot_map(gdf::Union{DataFrame,DataFrameRow}; geom_col::Symbol=:geometry)
    f = Figure(; size=(600, 900))
    ga = GeoAxis(
        f[1, 1];
        dest="+proj=latlong +datum=WGS84",
        xlabel="Longitude",
        ylabel="Latitude",
        xticklabelpad=15,
        yticklabelpad=40,
        xticklabelsize=10,
        yticklabelsize=10,
        aspect=AxisAspect(0.75),
        xgridwidth=0.5,
        ygridwidth=0.5,
    )

    plottable = _convert_plottable(gdf, geom_col)
    poly!(ga, plottable)

    display(f)

    return f, ga
end

function plot_map!(ga::GeoAxis, gdf::DataFrame; geom_col=:geometry, color=nothing)::Nothing

    plottable = _convert_plottable(gdf, geom_col)
    if !isnothing(color)
        poly!(ga, plottable, color=color)
    else
        poly!(ga, plottable)
    end

    # Set figure limits explicitly
    xlims!(ga)
    ylims!(ga)

    return nothing
end

function plot_map!(gdf::DataFrame; geom_col=:geometry, color=nothing)::Nothing
    return plot_map!(current_axis(), gdf; geom_col=geom_col, color=color)
end

function plot_map(gdf::Union{DataFrame,DataFrameRow}, color_by::Symbol; geom_col::Symbol=:geometry)
    f = Figure(; size=(600, 900))
    ga = GeoAxis(
        f[1, 1];
        dest="+proj=latlong +datum=WGS84",
        xlabel="Longitude",
        ylabel="Latitude",
        xticklabelpad=15,
        yticklabelpad=40,
        xticklabelsize=10,
        yticklabelsize=10,
        aspect=AxisAspect(0.75),
        xgridwidth=0.5,
        ygridwidth=0.5,
    )

    plottable = _convert_plottable(gdf, geom_col)

    # Define the unique colors and names for each level of factor color_by.
    # Use a different color palette for factors with high numbers of levels
    # (this palette is not as good for visualisation).
    if size(unique(gdf[:, color_by]),1) <= 20
        palette = ColorSchemes.tableau_20.colors
    else
        palette = ColorSchemes.flag_ec.colors
    end

    color_indices = groupindices(DataFrames.groupby(gdf, color_by))
    names = unique(DataFrame(indices=color_indices, names=gdf[:, color_by]))

    # Create the unique legend entries for each level of color_by
    unique_names = names.names
    legend_entries = []
    for name in eachrow(names)
        col = palette[name.indices]
        LE = PolyElement(; color=col)
        push!(legend_entries, [LE])
    end

    polys = poly!(ga, plottable, color=palette[color_indices])

    Legend(f[2, 1], legend_entries, unique_names, nbanks=3, tellheight=true,
    tellwidth=false, orientation=:horizontal, labelsize=10)

    display(f)
    return f, ga
end

function colorscheme_alpha(cscheme::ColorScheme, alpha = 0.5)
    ncolors = length(cscheme)

    return ColorScheme([RGBA(get(cscheme, k), alpha) for k in range(0, 1, length=ncolors)])
end

"""
    plot_lines(
    df::DataFrame,
    color_by::Symbol,
    y_start_col,
    y_end_col,
    x_lab::String,
    y_lab::String,
    alpha=1
    )

Convenience plot function for plotting timeseries lines with categorical colours.

# Arguments
- `df` : DataFrame with wide format (each unique line {reef} in a plot occupies it's own row in df)
- `color_by` : Column name holding factor to color reefs by (e.g. :management_area)
- `y_start_col` : Column name holding the first timepoint of the series (often 1 or 2)
- `y_end_col` : Column name holding the last timepoint of the series (often 79)
- `x_lab` : Label for X Axis
- `y_lab` : Label for Y Axis
- `alpha` : Alpha value for transparency (0-1)
"""
function plot_lines(
    df::DataFrame,
    color_by::Symbol,
    y_start_col,
    y_end_col,
    x_lab::String,
    y_lab::String,
    alpha=1
    )

    f = Figure()
    ax = Axis(f[1,1]; xlabel=x_lab, ylabel=y_lab)

    # Define the unique colors and names for each level of factor color_by.
    # Use a different color palette for factors with high numbers of levels
    # (this palette is not as good for visualisation).
    if size(unique(df[:, color_by]),1) == 2
        alph_palette = colorscheme_alpha(ColorSchemes.tableau_20, alpha)[2:3]
        palette = ColorSchemes.tableau_20.colors[2:3]
    elseif size(unique(df[:, color_by]),1) <= 20
        alph_palette = colorscheme_alpha(ColorSchemes.tableau_20, alpha)
        palette = ColorSchemes.tableau_20.colors
    else
        alph_palette = colorscheme_alpha(ColorSchemes.flag_ec, alpha)
        palette = ColorSchemes.flag_ec.colors
    end

    color_indices = groupindices(DataFrames.groupby(df, color_by))
    color_names = unique(DataFrame(indices=color_indices, names=df[:, color_by]))

    # Create the unique legend entries for each level of color_by
    unique_names = color_names.names
    legend_entries = []
    for name in eachrow(color_names)
        col = palette[name.indices]
        LE = PolyElement(; color=col)
        push!(legend_entries, [LE])
    end

    lines = series!(ax, Matrix(DataFrames.select(df, DataFrames.Between(Symbol(y_start_col), Symbol(y_end_col)))), color=alph_palette[color_indices])

    Legend(f[2, 1], legend_entries, unique_names, nbanks=3, tellheight=true,
    tellwidth=false, orientation=:horizontal, labelsize=10)

    display(f)
    return f, ax
end

# Checking how prominent the reefs with correlation at later lags are at earlier lags
using Colors, Random

# For horizintal and vertical text alignment:
justifyme(θ) = (0≤θ<π/2 || 3π/2<θ≤2π) ? :left : (π/2<θ<3π/2) ? :right : :center
justifymeV(θ) = π/4≤θ≤3π/4 ? :bottom : 5π/4<θ≤7π/4 ? :top : :center

# Radar plot function:
function radarplot(ax::Axis, v, val_labels; p_grid = maximum(v) * (1.0:length(val_labels)) / length(val_labels), title = "", labels = eachindex(v), labelsize = 1, points=true, spokeswidth= 1.5, spokescolor=:salmon, fillalpha=0.2, linewidth=1.5)
    # Axis attributes
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xminorgridvisible = false
    ax.yminorgridvisible = false
    ax.leftspinevisible = false
    ax.rightspinevisible = false
    ax.bottomspinevisible = false
    ax.topspinevisible = false
    ax.xminorticksvisible = false
    ax.yminorticksvisible = false
    ax.xticksvisible = false
    ax.yticksvisible = false
    ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
    ax.aspect = DataAspect()
    ax.title = title
    #
    l = length(v)
    rad = (0:(l-1)) * 2π / l
    # Point coordinates
    x = v .* cos.(rad) .- 0.05
    y = v .* sin.(rad) .- 0.05
    # if p_grid != 0
        # Coordinates for radial grid
        xa = maximum(p_grid) * cos.(rad) * 1.1
        ya = maximum(p_grid) * sin.(rad) * 1.1
        # Coordinates for polar grid text
        radC = (rad[Int(round(l / 2))] + rad[1 + Int(round(l / 2))]) / 2.0
        xc = p_grid * cos(radC)
        yc = p_grid * sin(radC)
        for i in p_grid
            poly!(ax, Circle(Point2f(0, 0), i), color = :transparent, strokewidth=1, strokecolor=ax.xgridcolor)
        end
        text!(ax, xc, yc, text=val_labels, fontsize = 12, align = (:center, :baseline), color=ax.xlabelcolor)
        arrows!(ax, zeros(l), zeros(l), xa, ya, color=ax.xgridcolor, linestyle=:solid, arrowhead=' ')
        if length(labels) == l
            for i in eachindex(rad)
                text!(ax, xa[i], ya[i], text=string(labels[i]), fontsize = labelsize, markerspace = :data, align = (justifyme(rad[i]), justifymeV(rad[i])), color=ax.xlabelcolor)
            end
        elseif length(labels) > 1
            printstyled("WARNING! Labels omitted:  they don't match with the points ($(length(labels)) vs $l).\n", bold=true, color=:yellow)
        end
    # end
    pp = scatter!(ax, [(x[i], y[i]) for i in eachindex(x)], color=RGB{Float32}(0.4, 0.4, 0.4))
    cc = to_value(pp.color)
    m_color = RGBA{Float32}(comp1(cc), comp2(cc), comp3(cc), fillalpha)
    s_color = RGB{Float32}(comp1(cc), comp2(cc), comp3(cc))
    pp.color = m_color
    pp.strokecolor = s_color
    pp.strokewidth= linewidth
    arrows!(ax, zeros(l), zeros(l), x, y, color=spokescolor, linewidth=spokeswidth, arrowhead=' ')
    if points
        scatter!(ax, x, y)
    end
    ax
end

function radarplot!(ax, v)

    l = length(v)
    rad = (0:(l-1)) * 2π / l

    x = v .* cos.(rad) .- 0.15
    y = v .* sin.(rad) .- 0.15

    pp = scatter!(ax, [(x[i], y[i]) for i in eachindex(x)], color=RGB{Float32}(0.4, 0.4, 0.4))

    ax
end

"""
    radarplot_df(df, cols, val_labels, labels; spokeswidth = 0, labelsize = 1)

Plots a radar plot intended to visualise each of the bellwether reefs and check whether their
correlation occurs at all lags or only some lags.

# Arguments
- `df` : dataframe containing a row for each reef and columns for each relevant lag. Each
lag column has corresponding integer values if a reef is a bellwether reef at that lag, e.g.
if reef-1 is a bellwether reef at lag5 the value of df[reef-1, lag5] = 5.
- `cols` : vector of column names containing the lag values for each reef.
- `val_labels` : labels for the values of each ring in the plot.
- `labels` : Labels for each reef. (RME_UNIQUE_ID or similar identifier)
"""
function radarplot_df(df, cols, val_labels, labels; spokeswidth = 0, labelsize = 1)
    fig = Figure()
    ax = Axis(fig[1,1])
    max_vals = [maximum(df[:, col]) for col in cols]
    initial_column = df[:, cols[argmax(max_vals)]]

    f = radarplot(ax, initial_column, val_labels; labels = labels, spokeswidth = spokeswidth, labelsize = labelsize)

    for col in cols
        f = radarplot!(ax, df[:, col])
    end

    display(fig)

    return fig
end

"""
    bellwether_reef_numbers_plot(df, variable_col, value_col, xlabel="lags", ylabel="number of bellwether reefs")

Intended to plot the numbers of bellwether reefs seen in each analysis level on a scatter plot.

# Arguments
- `df` : A long format dataframe with column lags containing Integers for each relevant lag,
a variable column with the analysis level names for each lag, a value column with the numbers of bellwether reefs.
- `variable_col` : name of column containing the analysis level strings for colouring points.
- `value_col` : name of column containing the numbers of bellwether reefs for each lag/analysis level.
"""
function bellwether_reef_numbers_plot(df, variable_col, value_col, xlabel="lags", ylabel="number of bellwether reefs")
    f = Figure()
    ax = Axis(f[1,1]; xlabel=xlabel, ylabel=ylabel)

    alph_palette = colorscheme_alpha(ColorSchemes.tableau_20, 0.7);
    palette = ColorSchemes.tableau_20.colors;

    color_indices = groupindices(DataFrames.groupby(df, variable_col))
    color_names = unique(DataFrame(indices=color_indices, names=df[:, variable_col]))

    # Create the unique legend entries for each level of variable_col
    unique_names = color_names.names
    legend_entries = []
    for name in eachrow(color_names)
        col = palette[name.indices]
        LE = PolyElement(; color=col)
        push!(legend_entries, [LE])
    end

    points = scatter!(df.lags, df[:, value_col], color=alph_palette[color_indices], markersize=10)

    Legend(f[2, 1], legend_entries, unique_names, nbanks=3, tellheight=true,
    tellwidth=false, orientation=:horizontal, labelsize=10)

    display(f)
end

function basic_target_reef_boxplot(
    categories, values;
    xlabel="Bellwether Reefs",
    ylabel="values",
    xticks=unique(categories),
    title="",
    method=nothing
)
    fig = Figure()
    target_reefs_axis = Axis(
        fig[1,1];
        xlabel = xlabel,
        ylabel = ylabel,
        xticks = (unique(categories.refs), xticks),
        title=title
    )

    colors = [Makie.wong_colors(); Makie.wong_colors()];
    legend_entries = []
    for (i, col) in enumerate(unique(colors[indexin(categories, unique(categories))]))
        LE = PolyElement(; color=col)
        push!(legend_entries, [LE])
    end

    if method == "rainclouds"
        f = rainclouds!(target_reefs_axis, categories.refs, values; color=colors[indexin(categories, unique(categories))], markersize=3, jitter_width=0.1, side_nudge=0.27)
    else
        f = boxplot!(target_reefs_axis, categories.refs, values; color=colors[indexin(categories, unique(categories))])
    end

    Legend(fig[1,2], legend_entries, unique(categories))

    display(fig)

    return fig
end

function basic_target_reef_boxplot(
    categories, values, yticks;
    xlabel="Bellwether Reefs",
    ylabel="values",
    xticks=["non-bellwether reefs", "bellwether reefs"],
    title=""
)
    fig = Figure()
    target_reefs_axis = Axis(
        fig[1,1];
        xlabel = xlabel,
        ylabel = ylabel,
        xticks = ([1.0,length(xticks)],xticks),
        yticks = yticks,
        title=title
    )
    f = boxplot!(target_reefs_axis, categories, values)

    display(fig)

    return fig
end

function timeseries_plot(data, target_reefs_col, length, ylabel, alpha)
    target_reefs = Matrix(
        DataFrames.select(
            data[data[:, target_reefs_col] .== "bellwether", :],
            DataFrames.Between(Symbol(length[1]), Symbol(length[2]))
        )
    )'
    target_reefs = median(target_reefs, dims=2)
    non_target_reefs = Matrix(
        DataFrames.select(
            data[data[:, target_reefs_col] .== "non-bellwether", :],
            DataFrames.Between(Symbol(length[1]), Symbol(length[2]))
        )
    )'
    non_target_reefs = median(non_target_reefs, dims=2)

    f, ax = plot_lines(data, target_reefs_col, length[1], length[2], "Year", ylabel, alpha)
    series!(target_reefs', solid_color=:red)
    series!(non_target_reefs', solid_color=:blue)

    return f
end

function skipmissing_median(x)
    new_median = Vector{Union{Missing, Float64}}(missing, size(x, 1))
    for (i, row) in enumerate(eachrow(x))
        if any(ismissing.(row))
            continue
        end
        new_median[i] = median(row)
    end

    return new_median
end

function lagged_timeseries_plot(data, target_reefs_col, length_t, xlab, ylab, alpha, lag)
    target_reefs = Matrix(
        DataFrames.select(
            data[data[:, target_reefs_col] .== "bellwether", :],
            DataFrames.Between(Symbol(length_t[1]), Symbol(length_t[2]))
        )
    )'
    target_reefs_buffer = fill(missing, lag, size(target_reefs, 2))
    target_reefs = vcat(target_reefs_buffer, target_reefs)

    target_reefs_median = skipmissing_median(target_reefs)

    non_target_reefs = Matrix(
        DataFrames.select(
            data[data[:, target_reefs_col] .== "non-bellwether", :],
            DataFrames.Between(Symbol(length_t[1]), Symbol(length_t[2]))
        )
    )'
    non_target_reefs_median = median(non_target_reefs, dims=2)


    f = Figure()
    ax = Axis(f[1,1]; xlabel=xlab, ylabel=ylab)

    alph_palette = colorscheme_alpha(ColorSchemes.tableau_20, alpha)[2:3];
    palette = ColorSchemes.tableau_20.colors[2:3]

    legend_entries = [
        PolyElement(;color = palette[1]),
        PolyElement(;color = palette[2])
    ]

    lines = series!(ax, target_reefs', solid_color=alph_palette[2])
    series!(non_target_reefs', solid_color=alph_palette[1])
    series!(target_reefs_median', solid_color=:red)
    series!(non_target_reefs_median', solid_color=:blue)

    Legend(f[2, 1], legend_entries, ["non-bellwether reefs", "bellwether reefs"], nbanks=3, tellheight=true,
    tellwidth=false, orientation=:horizontal, labelsize=10)

    return f
end

function lagged_timeseries_plot_combined(data, target_reefs_col, length_t, xlab, ylab, alpha, lag)
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size = (1000, 700))
    ga = Axis(f[1, 1], xlabel=xlab, ylabel=ylab)
    gb = Axis(f[1, 2], xlabel=xlab, ylabel=ylab)

    alph_palette = colorscheme_alpha(ColorSchemes.tableau_20, alpha)[2:3];
    palette = ColorSchemes.tableau_20.colors[2:3];

    legend_entries = [
        [PolyElement(;color = palette[1]), LineElement(;color=:blue)],
        [PolyElement(;color = palette[2]), LineElement(;color=:red)]
    ]

    # Prepare plotting matrices for non-lagged plot
    target_reefs = Matrix(
        DataFrames.select(
            data[data[:, target_reefs_col] .== "bellwether", :],
            DataFrames.Between(Symbol(length_t[1]), Symbol(length_t[2]))
        )
    )'
    target_reefs_median = skipmissing_median(target_reefs)

    non_target_reefs = Matrix(
        DataFrames.select(
            data[data[:, target_reefs_col] .== "non-bellwether", :],
            DataFrames.Between(Symbol(length_t[1]), Symbol(length_t[2]))
        )
    )'
    non_target_reefs_median = median(non_target_reefs, dims=2)

    # Plot non-lagged timeseries lines
    line_a = series!(ga, target_reefs', solid_color=(palette[2], 0.3))
    series!(ga, non_target_reefs', solid_color=alph_palette[1])
    series!(ga, target_reefs_median', solid_color=:red)
    series!(ga, non_target_reefs_median', solid_color=:blue)

    # Prepare bellwether reefs lagged timeseries data
    target_reefs = Matrix(
        DataFrames.select(
            data[data[:, target_reefs_col] .== "bellwether", :],
            DataFrames.Between(Symbol(length_t[1]), Symbol(length_t[2]))
        )
    )'
    target_reefs_buffer = fill(missing, lag, size(target_reefs, 2))
    target_reefs = vcat(target_reefs_buffer, target_reefs)
    target_reefs_median = skipmissing_median(target_reefs)

    # Plot lagged timeseries lines
    line_b = series!(gb, target_reefs', solid_color=(palette[2], 0.3))
    series!(gb, non_target_reefs', solid_color=alph_palette[1])
    series!(gb, target_reefs_median', solid_color=:red)
    series!(gb, non_target_reefs_median', solid_color=:blue)

    Legend(f[2, 1:2], legend_entries, ["non-bellwether reefs", "bellwether reefs"], nbanks=3, tellheight=true,
    tellwidth=false, orientation=:horizontal, labelsize=10, position=(0.8,0.3))

    Label(f[1,1,TopLeft()], "a",
        fontsize = 26,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right
    )
    Label(f[1,2,TopLeft()], "b",
        fontsize = 26,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right
    )

    return f
end

function explore_regression_scatter_plots(x, y; color=:blue, xlab="", ylab="", title="")

    fig = Figure()
    ax = Axis(
        fig[1,1];
        xlabel = xlab,
        ylabel = ylab,
        title = title
    )
    f = scatter!(
        x, y; color=color
    )

    return fig
end

function combined_bellwether_reef_boxplot(
    dataframe, bellwether_reefs_col;
    xlabel="Bellwether Reefs",
    xticks=unique(dataframe[:, bellwether_reefs_col]),
    title="",
    method=nothing
)
    dataframe = sort(dataframe, bellwether_reefs_col)
    categories = categorical(dataframe[:, bellwether_reefs_col])

    fig = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size = (1000, 700))
    ga = Axis(
        fig[1,1];
        xlabel = xlabel,
        ylabel = "mean DHW",
        xticks = (unique(categories.refs), xticks),
        title=title,
        width=350,
        height=200
    )
    gb = Axis(
        fig[1,2];
        xlabel = xlabel,
        ylabel = "Source to sink ratio",
        xticks = (unique(categories.refs), xticks),
        title=title,
        width=350,
        height=200
    )
    gc = Axis(
        fig[1,3];
        xlabel = xlabel,
        ylabel = "Connectivity eigenvector centrality",
        xticks = (unique(categories.refs), xticks),
        title=title,
        width=350,
        height=200
    )
    gd = Axis(
        fig[2,1];
        xlabel = xlabel,
        ylabel = "Total connectivity strength",
        xticks = (unique(categories.refs), xticks),
        title=title,
        width=350,
        height=200
    )
    ge = Axis(
        fig[2,2];
        xlabel = xlabel,
        ylabel = "Initial coral cover",
        xticks = (unique(categories.refs), xticks),
        title=title,
        width=350,
        height=200
    )
    resize_to_layout!(fig)

    colors = [Makie.wong_colors(); Makie.wong_colors()];
    legend_entries = []
    for (i, col) in enumerate(unique(colors[indexin(categories, unique(categories))]))
        LE = PolyElement(; color=col)
        push!(legend_entries, [LE])
    end

    f = rainclouds!(ga, categories.refs, dataframe.mean_dhw; color=colors[indexin(categories, unique(categories))], markersize=3, jitter_width=0.1, side_nudge=0.27)
    f = rainclouds!(gb, categories.refs, dataframe.so_to_si; color=colors[indexin(categories, unique(categories))], markersize=3, jitter_width=0.1, side_nudge=0.27)
    f = rainclouds!(gc, categories.refs, dataframe.conn_score; color=colors[indexin(categories, unique(categories))], markersize=3, jitter_width=0.1, side_nudge=0.27)
    f = rainclouds!(gd, categories.refs, dataframe.total_strength; color=colors[indexin(categories, unique(categories))], markersize=3, jitter_width=0.1, side_nudge=0.27)
    f = rainclouds!(ge, categories.refs, dataframe.initial_coral_cover; color=colors[indexin(categories, unique(categories))], markersize=3, jitter_width=0.1, side_nudge=0.27)

    Label(fig[1,1,TopLeft()], "a",
        fontsize = 26,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right
    )
    Label(fig[1,2,TopLeft()], "b",
        fontsize = 26,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right
    )
    Label(fig[1,3,TopLeft()], "c",
        fontsize = 26,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right
    )
    Label(fig[2,1,TopLeft()], "d",
        fontsize = 26,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right
    )
    Label(fig[2,2,TopLeft()], "e",
        fontsize = 26,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :right
    )

    resize_to_layout!(fig)
    Legend(fig[2, 3], legend_entries, unique(categories), nbanks=2, padding = 15)

    display(fig)

    return fig
end

function bioregion_grouped_boxplots(
    dataframe, bellwether_reefs_col, grouping, variable, ncol;
    xlabel="Bellwether Reefs",
    ylabel="Value",
    xticks=unique(dataframe[:, bellwether_reefs_col]),
    title="",
    method=nothing
)
    dataframe = sort(dataframe, [bellwether_reefs_col, :management_area]; rev=true)
    categories = categorical(dataframe[:, bellwether_reefs_col])

    gdf = DataFrames.groupby(dataframe, grouping)

    if length(gdf) < 27
        labels = string.(collect('a':'z'))[1:length(gdf)]
    else
        error("number of bioregions > 26 (number of letters for labelling)")
    end

    fig = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size = (2130, 1500))

    colors = [Makie.wong_colors(); Makie.wong_colors()];
    legend_entries = []
    for (i, col) in enumerate(unique(colors[indexin(categories, unique(categories))]))
        LE = MarkerElement(; color=col, marker=:circle)
        push!(legend_entries, [LE])
    end

    min_var, max_var = extrema(dataframe[:, variable])
    var_difference = max_var - min_var
    # if variable ∈ [:mean_dhw, :so_to_si]
    #     limits = (nothing, (floor(Int64, min_var - 2), ceil(Int64, max_var + 2)))
    #     yticks = round.(Int64, min_var:2:max_var)
    # elseif variable == :initial_coral_cover
    #     limits = (nothing, (floor(Int64, min_var - 5), ceil(Int64, max_var + 5)))
    #     yticks = round.(Int64, min_var:10:max_var)
    # elseif variable == :total_strength
    #     limits = (nothing, (round(min_var - 0.3, digits=2), round(max_var + 0.3, digits=2)))
    #     yticks = round.(min_var:0.5:max_var, digits=1)
    # elseif variable ∈ [:dhw_cover_cor, :dhw_evenness_cor]
    #     limits = (nothing, (-1.1,1.1))
    #     yticks = (-1:0.5:1)
    # elseif variable == :conn_score
    #     limits = (nothing, nothing)
    #     yticks = Makie.automatic
    # end

    if var_difference < 4
        limits = (nothing, (min_var - (var_difference/5), max_var + (var_difference/5)))
        yticks = round.(min_var:(var_difference / 4):max_var; digits=3)
    else
        limits = (nothing, (floor(Int64, min_var - 1), ceil(Int64, max_var + 1)))
        yticks = round.(min_var:(var_difference / 4):max_var; digits=1)
    end

    if length(gdf) < 10
        xsize, ysize = 180, 120
        label_size = 12
    else
        xsize, ysize = 150, 150
        label_size = 12
    end

    management_areas = Dict(
        "Far Northern Management Area" => (:red, 0.7),
        "Cairns/Cooktown Management Area" => (:green, 0.7),
        "Townsville/Whitsunday Management Area" => (:blue, 0.7),
        "Mackay/Capricorn Management Area" => (:Orange, 0.7)
    )

    for (xi, groupdf) in enumerate(gdf)
        categories = categorical(groupdf[:, bellwether_reefs_col])

        bellwether_var_med = median(groupdf[groupdf[:, bellwether_reefs_col] .== "bellwether", variable])
        non_bellwether_var_75 = quantile(groupdf[groupdf[:, bellwether_reefs_col] .== "non-bellwether", variable], 0.75)
        non_bellwether_var_25 = quantile(groupdf[groupdf[:, bellwether_reefs_col] .== "non-bellwether", variable], 0.25)

        if bellwether_var_med < non_bellwether_var_25
            background_color = (:royalblue, 0.2)
        elseif bellwether_var_med > non_bellwether_var_75
            background_color = (:goldenrod1, 0.2)
        else
            background_color = :white
        end

        if xi == 1
            ax = Axis(
                fig[fldmod1(xi, ncol)...];
                backgroundcolor=background_color,
                xlabel = xlabel,
                xlabelsize = label_size+1,
                ylabelsize = label_size+1,
                ylabel = ylabel,
                xticks = (unique(categories.refs), xticks),
                xticklabelsize=label_size,
                limits= limits,
                yticks = yticks,
                yticklabelsize=label_size,
                title=title,
                width=xsize,
                height=ysize
            )
        else
            ax = Axis(
                fig[fldmod1(xi, ncol)...];
                backgroundcolor=background_color,
                xticklabelsize=label_size,
                yticklabelsize=label_size,
                xticks = (unique(categories.refs), ["",""]),
                xticksvisible=false,
                limits= limits,
                yticks = yticks,
                title=title,
                width=xsize,
                height=ysize
            )
        end
        f = rainclouds!(
            ax,
            categories.refs,
            groupdf[:, variable];
            color=colors[indexin(categories, unique(categories))],
            jitter_width=0.27,
            markersize=3.5,
            cloudwidth=0.75,
            gap=0.01,
            side_nudge=0.32,
        )

        if variable == :so_to_si
            hlines!(1; color=(:gray, 0.7))
        end

        bioregion_label_color = management_areas[first(unique(groupdf.management_area))]
        Label(
            fig[fldmod1(xi, ncol)..., TopLeft()],
            labels[xi],
            color = bioregion_label_color,
            fontsize = 15,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right
        )
    end

    last_figure = fldmod1(length(gdf), ncol)
    if last_figure[2] < ncol
        Legend(fig[last_figure[1], last_figure[2]+1:ncol], legend_entries, unique(categories), nbanks=1)
    else
        Legend(fig[last_figure[1] + 1, 1:ncol], legend_entries, unique(categories), nbanks=1)
    end
    resize_to_layout!(fig)

    display(fig)

    return fig
end


function bioregion_grouped_lagged_timeseries(
    dataframe, bellwether_reefs_col, length_t, grouping, lag, ncol;
    xlabel="Years",
    ylabel="Value",
    title=""
)
    dataframe = sort(dataframe, [bellwether_reefs_col, :management_area]; rev=true)
    categories = categorical(dataframe[:, bellwether_reefs_col])

    gdf = DataFrames.groupby(dataframe, grouping)

    if length(gdf) < 27
        labels = string.(collect('a':'z'))[1:length(gdf)]
    else
        error("number of bioregions > 26 (number of letters for labelling)")
    end

    fig = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98), size = (2130, 1500))

    colors = [Makie.wong_colors(); Makie.wong_colors()];
    legend_entries = [
        [PolyElement(;color = colors[1]), LineElement(;color=:blue)],
        [PolyElement(;color = colors[2]), LineElement(;color=:red)]
    ]

    min_var, max_var = extrema(
        Matrix(DataFrames.select(dataframe, DataFrames.Between(
            Symbol(length_t[1]), Symbol(length_t[2]))
        ))
    )
    var_difference = max_var - min_var
    if bellwether_reefs_col == :target_reefs_bior_evenness_cat
        limits = (nothing, (min_var, max_var))
        yticks = round.(min_var:(var_difference / 3):max_var; digits=2)
    else
        limits = (nothing, (floor(Int64, min_var - 2), ceil(Int64, max_var + 2)))
        yticks = round.(Int64, min_var:(var_difference / 3):max_var)
    end


    if ncol < 5
        xsize, ysize = 180, 110
        label_size = 12
    else
        xsize, ysize = 150, 140
        label_size = 12
    end

    management_areas = Dict(
        "Far Northern Management Area" => (:red, 0.7),
        "Cairns/Cooktown Management Area" => (:green, 0.7),
        "Townsville/Whitsunday Management Area" => (:blue, 0.7),
        "Mackay/Capricorn Management Area" => (:Orange, 0.7)
    )

    for (xi, groupdf) in enumerate(gdf)
        categories = categorical(groupdf[:, bellwether_reefs_col])

        target_reefs = Matrix(
        DataFrames.select(
            groupdf[groupdf[:, bellwether_reefs_col] .== "bellwether", :],
            DataFrames.Between(Symbol(length_t[1]), Symbol(length_t[2]))
            )
        )'
        if lag > 0
            target_reefs_buffer = fill(missing, lag, size(target_reefs, 2))
            target_reefs = vcat(target_reefs_buffer, target_reefs)
        end
        target_reefs_median = skipmissing_median(target_reefs)

        non_target_reefs = Matrix(
            DataFrames.select(
                groupdf[groupdf[:, bellwether_reefs_col] .== "non-bellwether", :],
                DataFrames.Between(Symbol(length_t[1]), Symbol(length_t[2]))
            )
        )'
        non_target_reefs_median = median(non_target_reefs, dims=2)

        if xi == 1
            ax = Axis(
                fig[fldmod1(xi, ncol)...];
                xlabel = xlabel,
                xlabelsize = label_size+1,
                ylabelsize = label_size+1,
                ylabel = ylabel,
                xticklabelsize=label_size,
                limits= limits,
                yticks = yticks,
                yticklabelsize=label_size,
                title=title,
                width=xsize,
                height=ysize
            )
        else
            ax = Axis(
                fig[fldmod1(xi, ncol)...];
                xticklabelsize=label_size,
                limits= limits,
                yticks = yticks,
                yticklabelsize=label_size,
                title=title,
                width=xsize,
                height=ysize
            )
        end

        lines = series!(ax, target_reefs', solid_color=(colors[2], 0.3))
        series!(non_target_reefs', solid_color=(colors[1], 0.15))
        series!(target_reefs_median', solid_color=:red)
        series!(non_target_reefs_median', solid_color=:blue)

        bioregion_label_color = management_areas[first(unique(groupdf.management_area))]
        Label(
            fig[fldmod1(xi, ncol)..., TopLeft()],
            labels[xi],
            color = bioregion_label_color,
            fontsize = 15,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right
        )
    end

    last_figure = fldmod1(length(gdf), ncol)
    if last_figure[2] < ncol
        Legend(fig[last_figure[1], last_figure[2]+1:ncol], legend_entries, unique(categories), nbanks=1)
    else
        Legend(fig[last_figure[1] + 1, 1:ncol], legend_entries, unique(categories), nbanks=1)
    end
    resize_to_layout!(fig)

    display(fig)

    return fig, gdf.keymap
end

function plot_map_bellwether_reefs(gdf::Union{DataFrame,DataFrameRow}, color_by::Symbol, bellwether_reefs_col::Symbol; geom_col::Symbol=:geometry)
    if first(unique(gdf.management_area_short)) ∈ ["FarNorthern", "Cairns-Cooktown"]
        f = Figure(; size=(700, 1000))
    else
        f = Figure(; size=(1000, 800))
    end
    ga = GeoAxis(
        f[1, 1];
        dest="+proj=latlong +datum=WGS84",
        xlabel="Longitude",
        ylabel="Latitude",
        xticklabelpad=15,
        yticklabelpad=40,
        xticklabelsize=10,
        yticklabelsize=10,
    )

    plottable = _convert_plottable(gdf, geom_col)

    # Define the unique colors and names for each level of factor color_by.
    # Use a different color palette for factors with high numbers of levels
    # (this palette is not as good for visualisation).
    if size(unique(gdf[:, color_by]),1) <= 20
        palette = ColorSchemes.tableau_20.colors
        alph_palette = colorscheme_alpha(ColorSchemes.tableau_20, 0.6)
    else
        palette = ColorSchemes.flag_ec.colors
        alph_palette = colorscheme_alpha(ColorSchemes.flag_ec, 0.6)
    end

    color_indices = groupindices(DataFrames.groupby(gdf, color_by))
    names = unique(DataFrame(indices=color_indices, names=gdf[:, color_by]))

    # Create the unique legend entries for each level of color_by
    unique_names = names.names
    legend_entries = []
    for name in eachrow(names)
        col = palette[name.indices]
        LE = PolyElement(; color=col)
        push!(legend_entries, [LE])
    end
    bellwether_indices = gdf[:, bellwether_reefs_col] .== "bellwether"
    non_bellwether_indices = gdf[:, bellwether_reefs_col] .== "non-bellwether"

    non_bell_polys = poly!(ga, plottable[non_bellwether_indices], color=alph_palette[color_indices[non_bellwether_indices]])
    bell_polys = poly!(ga, plottable[bellwether_indices], color=palette[color_indices[bellwether_indices]], strokewidth=0.5)

    Legend(f[2, 1], legend_entries, unique_names, nbanks=3, tellheight=true,
    tellwidth=false, orientation=:horizontal, labelsize=9, patchsize=(11,11))
    resize_to_layout!(f)

    display(f)
    return f, ga
end
