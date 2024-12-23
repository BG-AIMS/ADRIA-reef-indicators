"""
Functions used in 3_analysis_context_layers to load data directly from RME instead of using ADRIA.RMEDomain as
RME uses an additional species etc.
"""

"""
    _get_relevant_files(fn_path, ident)

Filter files found in given path down to those that have the provided identifier.

# Arguments
- `fn_path` : directory of files
- `ident` : keep files with `ident` in their filename.
"""
function _get_relevant_files(fn_path::String, ident::String)
    valid_files = filter(isfile, readdir(fn_path; join=true))
    return filter(x -> occursin(ident, x), valid_files)
end

"""
    load_DHW(::Type{RMEDomain}, data_path::String, rcp::String, timeframe=(2022, 2100))::YAXArray

Loads ReefMod DHW data as a datacube.

# Arguments
- `RMEDomain`
- `data_path` : path to ReefMod data
- `rcp` : RCP identifier
- `timeframe` : range of years to represent.

# Returns
YAXArray[timesteps, locs, scenarios]
"""
function load_DHW(data_path::String, rcp::String, timeframe=(2022, 2099))::YAXArray
    dhw_path = joinpath(data_path, "dhw_csv")
    rcp_files = _get_relevant_files(dhw_path, rcp)
    rcp_files = filter(x -> occursin("SSP", x), rcp_files)
    if isempty(rcp_files)
        ArgumentError("No DHW data files found in: $(dhw_path)")
    end

    first_file = CSV.read(rcp_files[1], DataFrame)
    loc_ids = String.(first_file[:, 1])

    data_tf = parse.(Int64, names(first_file[:, 2:end]))
    tf_start = findall(timeframe[1] .∈ data_tf)[1]
    tf_end = findall(timeframe[2] .∈ data_tf)[1]

    d1 = first_file[:, (tf_start+1):(tf_end+1)]
    data_shape = reverse(size(d1))
    data_cube = zeros(data_shape..., length(rcp_files))
    data_cube[:, :, 1] .= Matrix(d1)'

    local tf_start
    local tf_end
    keep_ds = fill(true, length(rcp_files))
    for (i, rcp_data_fn) in enumerate(rcp_files[2:end])
        d = CSV.read(rcp_data_fn, DataFrame)[:, 2:end]

        if size(d, 1) == 0
            @info "Empty file?" rcp_data_fn
            continue
        end

        data_tf = parse.(Int64, names(d))
        try
            tf_start = findall(timeframe[1] .∈ data_tf)[1]
            tf_end = findall(timeframe[2] .∈ data_tf)[1]
        catch err
            if !(err isa BoundsError)
                rethrow(err)
            end

            @warn "Building DHW: Could not find matching time frame, skipping $rcp_data_fn"
            keep_ds[i] = false  # mark scenario for removal
            continue
        end

        data_cube[:, :, i+1] .= Matrix(d[:, tf_start:tf_end])'
    end

    # Only return valid scenarios
    return DataCube(
        data_cube[:, :, keep_ds];
        timesteps=timeframe[1]:timeframe[2],
        locs=loc_ids,
        scenarios=rcp_files[keep_ds],
    )
end

"""
    load_DHW(::Type{RMEDomain}, data_path::String, rcp::String, timeframe=(2022, 2100))::YAXArray

Loads ReefMod DHW data as a datacube.

# Arguments
- `RMEDomain`
- `data_path` : path to ReefMod data
- `rcp` : RCP identifier
- `timeframe` : range of years to represent.

# Returns
YAXArray[timesteps, locs, scenarios]
"""
function load_DHW(data_path::String, rcp::String, GCM::String, timeframe=(2022, 2099))::YAXArray
    dhw_path = joinpath(data_path, "dhw_csv")
    rcp_files = _get_relevant_files(dhw_path, rcp)
    rcp_files = filter(x -> occursin(GCM, x), rcp_files)
    if GCM != "GFDL-ESM4"
        rcp_files = filter(x -> occursin("SSP", x), rcp_files)
    end
    if isempty(rcp_files)
        ArgumentError("No DHW data files found in: $(dhw_path)")
    end

    first_file = CSV.read(rcp_files[1], DataFrame)
    loc_ids = String.(first_file[:, 1])

    data_tf = parse.(Int64, names(first_file[:, 2:end]))
    tf_start = findall(timeframe[1] .∈ data_tf)[1]
    tf_end = findall(timeframe[2] .∈ data_tf)[1]

    d1 = first_file[:, (tf_start+1):(tf_end+1)]
    data_shape = reverse(size(d1))
    data_cube = zeros(data_shape..., length(rcp_files))
    data_cube[:, :, 1] .= Matrix(d1)'

    local tf_start
    local tf_end
    keep_ds = fill(true, length(rcp_files))
    for (i, rcp_data_fn) in enumerate(rcp_files[2:end])
        d = CSV.read(rcp_data_fn, DataFrame)[:, 2:end]

        if size(d, 1) == 0
            @info "Empty file?" rcp_data_fn
            continue
        end

        data_tf = parse.(Int64, names(d))
        try
            tf_start = findall(timeframe[1] .∈ data_tf)[1]
            tf_end = findall(timeframe[2] .∈ data_tf)[1]
        catch err
            if !(err isa BoundsError)
                rethrow(err)
            end

            @warn "Building DHW: Could not find matching time frame, skipping $rcp_data_fn"
            keep_ds[i] = false  # mark scenario for removal
            continue
        end

        data_cube[:, :, i+1] .= Matrix(d[:, tf_start:tf_end])'
    end

    # Only return valid scenarios
    return DataCube(
        data_cube[:, :, keep_ds];
        timesteps=timeframe[1]:timeframe[2],
        locs=loc_ids,
        scenarios=rcp_files[keep_ds],
    )
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
