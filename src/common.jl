using Glob

using
    CSV,
    Dates,
    DataFrames,
    YAXArrays,
    DimensionalData,
    DataStructures

using
    ColorSchemes,
    GeoMakie,
    GLMakie,
    GraphMakie

using
    Statistics,
    Bootstrap,
    LinearAlgebra

using Graphs, SimpleWeightedGraphs
import Graphs.Parallel
import GeoDataFrames as GDF
import GeoFormatTypes as GFT
import ArchGDAL as AG
import GeoInterface as GI

gbr_domain_path = "c:/Users/bgrier/Documents/Projects/ADRIA_Domains/rme_ml_2024_01_08/"
conn_path = joinpath(gbr_domain_path, "data_files/con_bin")

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

# """
#     mapslices_toFloat64(func, data::YAXArray, dim::Symbol)

# Apply a function across grouped slices of a YAXArray.
# E.g. To reduce from 3 to 2 dimensions, use a summarising `func` and the dimension to be
# reduced as `dim`.
# E.g. To apply function over timesteps for each reef, use a non-summarising `func` and
# timesteps as `dim`.

# # Arguments
# - `func` : Function to apply to array slices
# - `data` : YAXArray for slicing (2D or 3D)
# - `dim` : Target dimension symbol (if reducing use unwanted dim, if transforming across each series use timesteps)

# # Returns
# YAXArray after applying function for use in further ADRIA analysis/viz functions.
# """
# function mapslices_toFloat64(func, data::YAXArray, dim::Symbol)
#     init_array = mapslices(func, data, dims=[dim])
#     dat = convert.(Float64, init_array.data)
#     new_array::YAXArray{Float64} = YAXArray(dims(init_array), dat)

#     return new_array
# end

"""
    cross_correlation(
        x::AbstractVector{<:Real},
        y::AbstractVector{<:Real},
        lags::AbstractVector{<:Integer}
    )

Calculate the normalised cross correlation of two vectors x and y with time series
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
Vector of correlation values for each lag in `lags`.
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

Calculate the normalised cross correlation of two vectors x and y with time series
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
Vector of correlation values for each lag in `lags`.
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
        zx::Vector{Float64} = x .- mean(x)
    else
        throw("`demean` must be true if included. Intended use for applying mean subtraction to `x` and `y`.")
    end

    if demean
        zy::Vector{Float64} = y .- mean(y)
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

Perform lagged cross correlation analysis across a number of clusters for a target
region. Uses `mapslices_toFloat64()` and `cross_correlation()` functions.

# Arguments
- `region_rel_cover` : YAXArray of cover trajectories for each reef in the target region.
- `region_clusters` : Vector of cluster categories assigned to each reef.
- `lags` : Vector of lags to apply with `cross_correlation()`. Positive lags test for `x` leading `y`, negative lags test for `y` leading `x`.

# Returns
DataFrame of correlation values for each reef in each cluster with columns for the reef RME_UNIQUE_ID, cluster category and correlation values for each lag step.
"""
function lagged_cluster_analysis(
        region_rel_cover::YAXArray,
        region_clusters::Vector{Int64},
        lags::AbstractVector{<:Integer}
    )::DataFrame

    lags_symbols = [Symbol("lag" * string(lags[a])) for a in eachindex(lags)]
    cross_cor = DataFrame(
        [Vector{Any}() for _ in 1:(length(lags)+2)],
        [:RME_UNIQUE_ID; :t_cluster ; lags_symbols]
    )

    for cluster in 1:n_clusters
        target_cluster = region_rel_cover[:, (region_clusters .== cluster)]

        cluster_median = Float64.(mapslices(median, target_cluster, dims=[:sites]))

        for (ind, reef) in enumerate(eachcol(target_cluster))
            reef_name = target_cluster.sites[ind]
            correlation = cross_correlation(reef, cluster_median, lags, true)

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

Perform lagged cross correlation analysis across a a target region.
Uses `mapslices_toFloat64()` and `cross_correlation()` functions.

# Arguments
- `region_rel_cover` : YAXArray of cover trajectories for each reef in the target region.
- `reg` : String of region name.
- `lags` : Vector of lags to apply with `cross_correlation()`. Positive lags test for `x` leading `y`, negative lags test for `y` leading `x`.

# Returns
DataFrame of correlation values for each reef in target_region with columns for the reef RME_UNIQUE_ID, region and correlation values for each lag step.
"""
function lagged_region_analysis(
        region_rel_cover::YAXArray,
        reg::String,
        lags::AbstractVector{<:Integer}
    )::DataFrame

    lags_symbols = [Symbol("lag" * string(lags[a])) for a in eachindex(lags)]
    cross_cor = DataFrame(
        [Vector{Any}() for _ in 1:(length(lags)+2)],
        [:RME_UNIQUE_ID; :region; lags_symbols]
    )

    reg_median = Float64.(mapslices(median, region_rel_cover, dims=[:sites]))

    for (ind, reef) in enumerate(eachcol(region_rel_cover))
        reef_name = region_rel_cover.sites[ind]
        correlation = cross_correlation(reef, reg_median, lags, true)

        push!(cross_cor, [reef_name; reg; correlation])
    end

    return cross_cor
end

"""
    subregion_analysis(
        subregions::Vector{String},
        rel_cover::YAXArray,
        context_layers::DataFrame,
        category::Symbol
    )::DataFrame

Apply lagged_region_analysis across a list of subregions such as closest_ports or bioregions.
Uses `lagged_region_analysis()` function, lags must be specified.

# Arguments
- `subregions` : Vector of subregion names to apply analysis to. Must be found in `context_layers.category` column.
- `rel_cover` : YAXArray of cover trajectories for each reef.
- `context_layers` : DataFrame containing columns `RME_UNIQUE_ID` and category column matching `subregions`.
- `category` : Column name in `context_layers` that contains `subregions` for each reef.
- `lags` : Vector of lags to apply with `cross_correlation()`. Positive lags test for `x` leading `y`, negative lags test for `y` leading `x`.

# Returns
DataFrame of correlation values for each reef and surrounding reefs in 'subregion'. Contains RME_UNIQUE_ID, region and correlation values for each lag step.
"""
function subregion_analysis(
    subregions::Vector{String},
    rel_cover::YAXArray,
    context_layers::DataFrame,
    category::Symbol,
    lags::AbstractVector{<:Integer}
    )::DataFrame

    lagged_analysis_sub = DataFrame()
    for subregion in subregions
        subregion_reefs = context_layers[(context_layers[:, category] .== subregion), :RME_UNIQUE_ID]
        subregion_cover = rel_cover[:, rel_cover.sites .∈ [subregion_reefs]]

        subregion_lagged_analysis = lagged_region_analysis(subregion_cover, subregion, lags)
        lagged_analysis_sub = vcat(lagged_analysis_sub, subregion_lagged_analysis)
    end

    return lagged_analysis_sub
end

"""
    DataCube(data::AbstractArray; kwargs...)::YAXArray

Constructor for YAXArray. When used with `axes_names`, the axes labels will be UnitRanges
from 1 up to that axis length.

# Arguments
- `data` : Array of data to be used when building the YAXArray
- `axes_names` :
"""
function DataCube(data::AbstractArray; kwargs...)::YAXArray
    return YAXArray(Tuple(Dim{name}(val) for (name, val) in kwargs), data)
end
function DataCube(data::AbstractArray, axes_names::Tuple)::YAXArray
    return DataCube(data; NamedTuple{axes_names}(1:len for len in size(data))...)
end

function load_connectivity(
    conn_path::String, loc_ids::Vector{String}
)::Tuple{YAXArray, YAXArray}
    conn_files = glob("*CONNECT_ACRO*", conn_path)
    if isempty(conn_files)
        ArgumentError("No CONNECT_ACRO data files found in: $(conn_path)")
    end

    n_locs = length(loc_ids)
    tmp_mat = zeros(n_locs, n_locs, length(conn_files))
    for (i, fn) in enumerate(conn_files)
        # File pattern used is "CONNECT_ACRO_[YEAR]_[DAY].bin"
        # We use a clunky regex approach to identify the year.
        # tmp = replace(split(fn, r"(?=CONNECT_ACRO_[0-9,4]+)")[2], "CONNECT_ACRO_"=>"")
        # year_id = split(tmp, "_")[1]

        # Turns out, there's only data for each year so just read in directly
        # Have to read in binary data - read first two values as Int32, and the rest
        # as Float32. Then reshape into a square (n_locs * n_locs) matrix.
        data = IOBuffer(read(fn))
        x = read(data, Int32)
        y = read(data, Int32)

        ds = Vector{Float32}(undef, x * y)
        tmp_mat[:, :, i] .= reshape(read!(data, ds), (n_locs, n_locs))
    end

    # Mean/stdev over all years
    mean_conn_data::Matrix{Float64} = dropdims(mean(tmp_mat; dims=3); dims=3)
    stdev_conn_data::Matrix{Float64} = dropdims(std(tmp_mat; dims=3); dims=3)

    mean_conn = DataCube(
        mean_conn_data;
        Source=loc_ids,
        Sink=loc_ids,
    )

    stdev_conn = DataCube(
        stdev_conn_data;
        Source=loc_ids,
        Sink=loc_ids,
    )

    return mean_conn, stdev_conn
end

canonical_reefs = find_latest_file("../canonical-reefs/output/")
canonical_reefs = GDF.read(canonical_reefs)

#MANAGEMENT_AREAS = unique(canonical_reefs.management_area)

# Get location indices
REGION_REEFS = DataStructures.OrderedDict(
    "Far Northern Management Area"=>Int64[],
    "Cairns/Cooktown Management Area"=>Int64[],
    "Townsville/Whitsunday Management Area"=>Int64[],
    "Mackay/Capricorn Management Area"=>Int64[],
)

for mgmt_area in collect(keys(REGION_REEFS))
    if mgmt_area == "NA"
        continue
    end

    REGION_REEFS[mgmt_area] = findall(canonical_reefs.management_area .== mgmt_area)
end

"""
    connectivity_scoring(
        conn_matrix::YAXArray;
        gdf::DataFrame=nothing,
        context_layer::Symbol=nothing,
        conn_col_name::Symbol=nothing,
        by_layer::Bool=false
    )::DataFrame

Calculate eigenvector_centrality connectivity scores for GBR reefs at different spatial scales.
When by_layer=true and a DataFrame is given for gdf and a layer symbol is given for context_layer,
then reefs are subset into categories via context_layer groups and then eigenvectors are calculated.

# Arguments
- `conn_matrix` : YAXArray containing the connectivity values. Diagonal should be set to 0 prior to use.
- `gdf` : DataFrame containing context_layer.
- `context_layer` : Categorical column with levels for grouping.
- `by_layer` : Calculate by grouping reefs by context_layer (vs whole GBR connectivity).
"""
function connectivity_scoring(
    conn_matrix::YAXArray;
    gdf=nothing,
    context_layer=nothing,
    conn_col_name=nothing,
    by_layer::Bool=false
)::DataFrame
    RME_UNIQUE_ID = collect(getAxis("Source", conn_matrix).val)
    connectivity_scores = DataFrame(
        [RME_UNIQUE_ID, Vector{Union{Missing, Float64}}(missing, 3806)],
        [:RME_UNIQUE_ID, :conn_score]
    )

    if by_layer
        for level in unique(gdf[:, context_layer])
            level_reefs = gdf[gdf[:, context_layer] .== level, :RME_UNIQUE_ID]
            level_matrix = conn_matrix[conn_matrix.Source .∈ [level_reefs], conn_matrix.Sink .∈ [level_reefs]]

            g = SimpleWeightedDiGraph(level_matrix)
            conn_scores = Dict(zip(collect(level_matrix.Source), eigenvector_centrality(g)))

            for (k, v) in conn_scores
                connectivity_scores[connectivity_scores.RME_UNIQUE_ID .== k, :conn_score] .= v
            end
        end

        rename!(connectivity_scores, :conn_score => conn_col_name)
    else
        g = SimpleWeightedDiGraph(conn_matrix)
        conn_score = eigenvector_centrality(g)
        connectivity_scores.conn_score = conn_score
    end

    return connectivity_scores
end

"""
    weight_by_context(
        gdf::DataFrame,
        target_col::Symbol,
        context_layer::Symbol,
        new_col_name::Symbol
    )::DataFrame

Rank the values in the target_col by their numerical order within each level of context_layer category.

# Arguments
- `gdf` : DataFrame containing all columns.
- `target_col` : Column containing values for ranking.
- `context_layer` : Categorical column with levels for ranking.
- `new_col_name` : Column name for new ranked values.
"""
function weight_by_context(
    gdf::DataFrame,
    target_col::Symbol,
    context_layer::Symbol,
    new_col_name::Symbol
    )::DataFrame

    gdf[:, new_col_name] .= 0.0
    for level in unique(gdf[:, context_layer])
        gdf_level = gdf[gdf[:, context_layer] .== level, :]

        max_target = maximum(gdf_level[:, target_col])
        gdf[gdf[:, context_layer] .== level, new_col_name] = gdf_level[:, target_col] ./ max_target
    end

    return gdf
end
