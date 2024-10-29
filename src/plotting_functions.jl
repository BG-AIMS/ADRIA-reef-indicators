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
