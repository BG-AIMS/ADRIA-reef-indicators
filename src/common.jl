using Glob

using
    CSV,
    Dates,
    DataFrames,
    YAXArrays

using
    ColorSchemes,
    GeoMakie,
    GLMakie

using
    Statistics,
    Bootstrap,
    LinearAlgebra

import GeoDataFrames as GDF
import GeoFormatTypes as GFT
import ArchGDAL as AG
import GeoInterface as GI

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

    color_indices = groupindices(groupby(gdf, color_by))
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
    if size(unique(df[:, color_by]),1) <= 20
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

"""
    find_intersections(
        x::DataFrame,
        y::DataFrame,
        x_id::Symbol,
        y_id::Symbol,
        y_geom_col::Symbol=:geometry;
        proportion::Bool=false,
    )::DataFrame

Find the areas of `y` that intersect with each polygon in `x`.
`rel_areas` contains corresponding `y_id` for each intersecting polygon in x (can then be
joined to `x`).

If `proportion = true`: polygons of `y` are only chosen if the intersection with `x` is >
50% the area of `x`.

# Arguments
- `x` : The target GeoDataFrame to compare with
- `y` : GeoDataFrame containing polygons to match against `x`
- `xid` : Column name holding unique IDs for x geometries (referred to as GBRMPA_ID in rel_areas)
- `yid` : Column name holding variable of interest for y geometries
- `y_geom_col` : Column name holding geometries in y
- `proportion` : Only select y polygons if the intersection with x polygon is > 50% of x polygon area
                 (default: `false`).
"""
function find_intersections(
    x::DataFrame,
    y::DataFrame,
    x_id::Symbol,
    y_id::Symbol,
    y_geom_col::Symbol=:geometry;
    proportion::Bool=false
)::DataFrame
    rel_areas = DataFrame(
        [Vector{Any}(missing, size(x, 1)) for _ in 1:2],
        [:GBRMPA_ID, :area_ID]
    )

    for (x_i, reef_poly) in enumerate(eachrow(x))
        intersecting = DataFrame(
            [Vector{Any}(missing, size(y, 1)) for _ in 1:3],
            [:GBRMPA_ID, :area_ID, :inter_area]
        )

        for (y_i, interest_area) in enumerate(eachrow(y))
            if AG.intersects(reef_poly.geometry, interest_area[y_geom_col])
                inter_area = AG.intersection(
                    reef_poly.geometry, interest_area[y_geom_col]
                )

                inter_area = AG.geomarea(inter_area)
                if proportion
                    prop_area = inter_area / AG.geomarea(reef_poly.geometry)

                    if prop_area >= 0.5
                        data = [reef_poly[x_id], interest_area[y_id], inter_area]

                    else
                        data = [missing, missing, missing]
                    end
                else
                    data = [reef_poly[x_id], interest_area[y_id], inter_area]
                end
            else
                data = [reef_poly[x_id], missing, missing]
            end

            intersecting[y_i, :] = data
        end

        if all(ismissing, intersecting.area_ID)
            x_data = [intersecting[1, x_id], intersecting[1, :area_ID]]
        else
            dropmissing!(intersecting)
            max_inter_area = argmax(intersecting.inter_area)
            x_data = [intersecting[max_inter_area, x_id], intersecting[max_inter_area, :area_ID]]
        end

        rel_areas[x_i, :] = x_data
    end

    return rel_areas
end

DATE_FORMAT = "YYYY-mm-dd-THH-MM-SS"

"""
    _get_file_timestamp(file_path, dt_length)::DateTime

Extract the timestamp from a given file name.
"""
function _get_file_timestamp(file_path, dt_length)::DateTime
    # Get name of file without extension
    filename = splitext(basename(file_path))[1]

    local fn_timestamp
    try
        fn_timestamp = Dates.DateTime(filename[end-(dt_length-1):end], DATE_FORMAT)
    catch err
        if !(err isa ArgumentError)
            rethrow(err)
        end

        # Otherwise, some unexpected date format was encountered so we assign an
        # very early date/time.
        fn_timestamp = Dates.DateTime("1900-01-01-T00-00-00", DATE_FORMAT)
    end

    # Return datetime stamp
    return fn_timestamp
end

"""
    find_latest_file(
        target_dir::String;
        prefix::String="rrap_canonical",
        ext::String="gpkg"
    )::String

Identify the latest output file in a directory based on the timestamp included in the
file name (default: `YYYY-mm-dd-THH-MM-SS`). Intended to find the latest output file for
input into the next script.

# Arguments
- `target_dir` : Target directory
- `prefix` : prefix of target file
- `ext` : the file extension

# Returns
Path to latest output file.
"""
function find_latest_file(
    target_dir::String;
    prefix::String="rrap_canonical",
    ext::String="gpkg"
)::String
    # Get list of files matching the given pattern
    candidate_files = glob("$(prefix)*.$(ext)", target_dir)

    timestamps = map(f -> _get_file_timestamp(f, length(DATE_FORMAT)), candidate_files)
    latest = candidate_files[argmax(timestamps)]

    return latest
end

"""
    relative_site_cover(x)

Calculate the cover of a site at each timestep of x relative to the site's initial cover (x[1]) (standardise all series to start at 1).
"""
function relative_site_cover(x)
    init = x[1]
    for (index, step) in enumerate(x)
        x[index] = x[index]/init
    end

    return x
end

"""
    mapslices_toFloat64(func, data::YAXArray, dim::Symbol)

Apply a function across grouped slices of a YAXArray.
E.g. To reduce from 3 to 2 dimensions, use a summarising `func` and the dimension to be
reduced as `dim`.
E.g. To apply function over timesteps for each reef, use a non-summarising `func` and
timesteps as `dim`.

# Arguments
- `func` : Function to apply to array slices
- `data` : YAXArray for slicing (2D or 3D)
- `dim` : Target dimension symbol (if reducing use unwanted dim, if transforming across each series use timesteps)

# Returns
YAXArray after applying function for use in further ADRIA analysis/viz functions.
"""
function mapslices_toFloat64(func, data::YAXArray, dim::Symbol)
    init_array = mapslices(func, data, dims=[dim])
    dat = convert.(Float64, init_array.data)
    new_array::YAXArray{Float64} = YAXArray(dims(init_array), dat)

    return new_array
end

"""
    cross_correlation(
        x::AbstractVector{<:Real},
        y::AbstractVector{<:Real},
        lags::AbstractVector{<:Integer}
    )

Function to calculate the normalised cross correlation of two vectors x and y with time series
lags. If `x` is ahead of `y` then a positive lag will result in positive correlation. If `y`
is ahead of `x`, then a negative lag will result in positive correlation.
E.g. If testing for x reef to be ahead of y reef, test for correlation at positive lag.

Based on StatsBase https://github.com/JuliaStats/StatsBase.jl/blob/60fb5cd400c31d75efd5cdb7e4edd5088d4b1229/src/signalcorr.jl#L400
and https://paulbourke.net/miscellaneous/correlate/

# Arguments
- `x` : Vector of interest to test for being ahead or behind `y`
- `y` : Vector to test lags of `x` against
- `lags` : Vector of lags to apply to vector. Positive lags test for `x` leading `y`, negative lags test for `y` leading `x`.

# Returns
`r` : Vector of correlation values for each lag in `lags`.
"""
function cross_correlation(
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    lags::AbstractVector{<:Integer}
    )

    r = Vector{Float64}()
    lx = length(x)
    m = length(lags)

    zx::Vector{Float64} = x
    zy::Vector{Float64} = y

    for k = 1 : m  # foreach lag value
        l = lags[k]

        if l >= 0
           sub_x = zx[1:lx-l]
           sub_y = zy[1+l:lx]
        else
           sub_x = zx[1-l:lx]
           sub_y = zy[1:lx+l]
        end

        sc = sqrt(dot(sub_x, sub_x) * dot(sub_y, sub_y))

        push!(r, dot(sub_x, sub_y) / sc)
    end

   return r
end

"""
    cross_correlation(
        x::AbstractVector{<:Real},
        y::AbstractVector{<:Real},
        lags::AbstractVector{<:Integer},
        demean::Bool
    )

Function to calculate the normalised cross correlation of two vectors x and y with time series
lags. If `x` is ahead of `y` then a positive lag will result in positive correlation. If `y`
is ahead of `x`, then a negative lag will result in positive correlation.
E.g. If testing for x reef to be ahead of y reef, test for correlation at positive lag.

Based on StatsBase https://github.com/JuliaStats/StatsBase.jl/blob/60fb5cd400c31d75efd5cdb7e4edd5088d4b1229/src/signalcorr.jl#L400
and https://paulbourke.net/miscellaneous/correlate/

# Arguments
- `x` : Vector of interest to test for being ahead or behind `y`
- `y` : Vector to test lags of `x` against
- `lags` : Vector of lags to apply to vector. Positive lags test for `x` leading `y`, negative lags test for `y` leading `x`.
- `demean` : Subtract the mean of each vector from each element of `x` and `y`. If demean is intended include it as true, otherwise do not include `demean` argument.

# Returns
`r` : Vector of correlation values for each lag in `lags`.
"""
function cross_correlation(
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    lags::AbstractVector{<:Integer},
    demean::Bool
    )

    r = Vector{Float64}()
    lx = length(x)
    m = length(lags)

    if demean
        zx::Vector{T} = x .- mean(x)
    else
        throw("`demean` must be true if included. Intended use for applying mean subtraction to `x` and `y`.")
    end

    if demean
        zy::Vector{S} = y .- mean(y)
    end

    for k = 1 : m  # foreach lag value
        l = lags[k]

        if l >= 0
           sub_x = zx[1:lx-l]
           sub_y = zy[1+l:lx]
        else
           sub_x = zx[1-l:lx]
           sub_y = zy[1:lx+l]
        end

        sc = sqrt(dot(sub_x, sub_x) * dot(sub_y, sub_y))

        push!(r, dot(sub_x, sub_y) / sc)
    end

   return r
end

"""
    lagged_cluster_analysis(
        region_rel_cover::YAXArray,
        region_clusters::Vector{Int64},
        lags::AbstractVector{<:Integer}
    )::DataFrame

Function to perform lagged cross correlation analysis across a number of clusters for a target
region. Uses `mapslices_toFloat64()` and `cross_correlation()` functions.

# Arguments
- `region_rel_cover` : YAXArray of cover trajectories for each reef in the target region.
- `region_clusters` : Vector of cluster categories assigned to each reef.
- `lags` : Vector of lags to apply with `cross_correlation()`. Positive lags test for `x` leading `y`, negative lags test for `y` leading `x`.

# Returns
DataFrame of correlation values for each reef in each cluster with columns for the reef UNIQUE_ID, cluster category and correlation values for each lag step.
"""
function lagged_cluster_analysis(
        region_rel_cover::YAXArray,
        region_clusters::Vector{Int64},
        lags::AbstractVector{<:Integer}
    )::DataFrame

    lags_symbols = [Symbol("lag" * string(lags[a])) for a in eachindex(lags)]
    cross_cor = DataFrame(
        [Vector{Any}() for _ in 1:(length(lags)+2)],
        [:UNIQUE_ID; :t_cluster ; lags_symbols]
    )

    for cluster in 1:n_clusters
        target_cluster = region_rel_cover[:, (region_clusters .== cluster)]

        cluster_median = mapslices_toFloat64(median, target_cluster, :sites)

        for (ind, reef) in enumerate(eachcol(target_cluster))
            reef_name = target_cluster.sites[ind]
            correlation = cross_correlation(reef, cluster_median, lags)

            push!(cross_cor, [reef_name; cluster; correlation])
        end
    end

    return cross_cor
end

"""
    lagged_region_analysis(
        region_rel_cover::YAXArray,
        reg::String,
        lags::AbstractVector{<:Integer}
    )::DataFrame

Function to perform lagged cross correlation analysis across a a target region.
Uses `mapslices_toFloat64()` and `cross_correlation()` functions.

# Arguments
- `region_rel_cover` : YAXArray of cover trajectories for each reef in the target region.
- `reg` : String of region name.
- `lags` : Vector of lags to apply with `cross_correlation()`. Positive lags test for `x` leading `y`, negative lags test for `y` leading `x`.

# Returns
DataFrame of correlation values for each reef in target_region with columns for the reef UNIQUE_ID, region and correlation values for each lag step.
"""
function lagged_region_analysis(
        region_rel_cover::YAXArray,
        reg::String,
        lags::AbstractVector{<:Integer}
    )::DataFrame

    lags_symbols = [Symbol("lag" * string(lags[a])) for a in eachindex(lags)]
    cross_cor = DataFrame(
        [Vector{Any}() for _ in 1:(length(lags)+2)],
        [:UNIQUE_ID; :region; lags_symbols]
    )

    reg_median = mapslices_toFloat64(median, region_rel_cover, :sites)

    for (ind, reef) in enumerate(eachcol(region_rel_cover))
        reef_name = region_rel_cover.sites[ind]
        correlation = cross_correlation(reef, reg_median, lags)

        push!(cross_cor, [reef_name; reg; correlation])
    end

    return cross_cor
end
