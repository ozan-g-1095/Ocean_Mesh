using NCDatasets
using JLD2
using Plots
using LaTeXStrings
using ProgressMeter
using Printf

const R_E = 6.371e6      # Earth's radius in meters
const Ω = 2π/(24*60*60)  # Earth's angular velocity in rad/s

"""
    z_filtered = filter_topo(lon, lat, z; L=1e5)

Filter the topography `z` on a (`lat`, `lon`) grid by averaging over circles of diameter `L`.

Land points (``z > 0``) are ignored.

See also [`circle_average`](@ref).
"""
function filter_topo(lon, lat, z; L=1e5, ocean_only=false)
    z_filtered = copy(z)
    @showprogress for i in axes(z, 1), j in axes(z, 2)
        if z[i, j] > 0 && ocean_only
            continue # skip land points
        end
        z_filtered[i, j] = circle_average(lon, lat, z, i, j; D=L)
    end
    return z_filtered
end

"""
    z_avg = circle_average(lon, lat, z, i0, j0; D)

Compute the average of `z` within a circle of diameter `D` centered at `(i0, j0)`.

The average is computed only over ocean points (``z ≤ 0``).

See also [`index_square_ring`](@ref) and [`sphere_distance`](@ref).
"""
function circle_average(lon, lat, z, i0, j0; D)
    lon0 = lon[i0]    # center longitude
    lat0 = lat[j0]    # center latitude
    z_tot = z[i0, j0] # initialize with the value at the center
    count = 1         # number of points within the circle
    n = 1             # starting index to add to (i0, j0)
    ocean_count = 1   # number of ocean points in the ring
    while ocean_count > 0
        ocean_count = 0
        for I in index_square_ring(i0, j0, n, length(lon), length(lat))
            if z[I] <= 0
                if sphere_distance(lon0, lat0, lon[I[1]], lat[I[2]]) <= D/2
                    z_tot += z[I]
                    count += 1
                    ocean_count += 1
                end
            end
        end
        n += 1
    end
    return z_tot/count
end

"""
    indices = index_square_ring(i0, j0, n, nlon, nlat)

Return a vector of indices in a square ring of size `n` centered at `(i0, j0)`.

The indices wrap around `nlon` and `nlat`.
"""
function index_square_ring(i0, j0, n, nlon, nlat)
    # list of all indices in a square ring of size n centered at (i0, j0)
    indices = Vector{CartesianIndex}(undef, 8n)
    idx = 1
    for i in (i0-n):(i0+n)
        if j0 - n ≥ 1
            indices[idx] = CartesianIndex(mod1(i, nlon), j0 - n)  # left edge
            idx += 1
        end
        if j0 + n ≤ nlat 
            indices[idx] = CartesianIndex(mod1(i, nlon), j0 + n)  # right edge
            idx += 1
        end
    end
    for j in (j0-n+1):(j0+n-1) # avoid double counting corners
        if j ≥ 1 && j ≤ nlat
            indices[idx] = CartesianIndex(mod1(i0 - n, nlon), j)  # bottom edge
            idx += 1
            indices[idx] = CartesianIndex(mod1(i0 + n, nlon), j)  # top edge
            idx += 1
        end
    end 
    indices = indices[1:idx-1]  # trim unused space if ring overlapped pole
    return indices
end

"""
    d = sphere_distance(lon1, lat1, lon2, lat2)

Compute the distance between two points on the sphere given their longitudes and latitudes.
"""
function sphere_distance(lon1, lat1, lon2, lat2)
    dlon = abs(lon2 - lon1)
    if dlon > 180
        dlon = 360 - dlon  # handle wrap-around at 180 degrees
    end
    dlat = abs(lat2 - lat1)
    dx = R_E*deg2rad(dlon)*cos(deg2rad(lat1))
    dy = R_E*deg2rad(dlat)
    return √(dx^2 + dy^2)
end

function plot_H(lon, lat, z, fname; vmax=6000)
    H = -z
    H[H .< 0] .= NaN
    heatmap(lon, lat, H', dpi=300,
            xlabel="Longitude (°)", ylabel="Latitude (°)",
            cb_title=L"$H$ (m)", aspect_ratio=:equal,
            clims=(0, vmax), 
            # c=cgrad(:magma, rev=true),
            c=:seaborn_rocket_gradient,
            grid=false, xlimits=(lon[1], lon[end]),
            ylimits=(lat[1], lat[end]), tickdirection=:out)
    savefig(fname)
    @info "Saved '$fname'"
end
function plot_H_panama(lon, lat, z, fname; vmax=6000)
    lon_min, lon_max = -95.0, -65.0
    lat_min, lat_max = 3.0, 18.0
    lon_p = lon[lon .>= lon_min .&& lon .<= lon_max]
    lat_p = lat[lat .>= lat_min .&& lat .<= lat_max]
    z_p = z[lon .>= lon_min .&& lon .<= lon_max, lat .>= lat_min .&& lat .<= lat_max]
    plot_H(lon_p, lat_p, z_p, fname; vmax=vmax)
end

function fill_panama!(lon, lat, z)
    lon_mask = lon .>= -80 .&& lon .<= -75
    lat_mask = lat .>= 8 .&& lat .<= 9
    z[lon_mask, lat_mask] .= NaN
    lon_mask = lon .>= -82 .&& lon .<= -76
    lat_mask = lat .>= 9 .&& lat .<= 10
    z[lon_mask, lat_mask] .= NaN
    lon_mask = lon .>= -83 .&& lon .<= -82
    lat_mask = lat .>= 10 .&& lat .<= 11
    z[lon_mask, lat_mask] .= NaN
    lon_mask = lon .>= -84 .&& lon .<= -83
    lat_mask = lat .>= 10 .&& lat .<= 13
    z[lon_mask, lat_mask] .= NaN
    return z
end

# load data
# name = "blacksea"; res = "5min"
name = "global"; res = "60min"
ifile = joinpath(@__DIR__, "etopo_$(name)_$(res).jld2")
@info "Loading data from '$ifile'..."
d = jldopen(ifile)
lon = d["lon"]
lat = d["lat"]
z = d["z"]
close(d)
@info "Dimensons" length(lon) length(lat) size(z)
# plot_H_panama(lon, lat, z, joinpath(@__DIR__, "H_$(name)_panama.png"))

# filter
L_land = 250e3 # 250 km
L_ocean = 500e3 # 500 km
z = filter_topo(lon, lat, z; L=L_land,  ocean_only=false) 
z = filter_topo(lon, lat, z; L=L_ocean, ocean_only=true) 
ofile = joinpath(@__DIR__, @sprintf("etopo_%s_%s_%dkm_%dkm.jld2", name, res, L_land/1e3, L_ocean/1e3))
jldsave(ofile; lon, lat, z)
@info "Saved '$ofile'"

# plot
ifile = joinpath(@__DIR__, @sprintf("etopo_%s_%s_%dkm_%dkm.jld2", name, res, L_land/1e3, L_ocean/1e3))
@info "Loading filtered data from '$ifile'..."
d = jldopen(ifile)
lon = d["lon"]
lat = d["lat"]
z = d["z"]
close(d)
plot_H(lon, lat, z, joinpath(@__DIR__, "H_$(name)_filtered.png"))

# fill Panama region
plot_H_panama(lon, lat, z, joinpath(@__DIR__, "H_$(name)_filtered_panama.png"))
z = fill_panama!(lon, lat, z)
plot_H(lon, lat, z, joinpath(@__DIR__, "H_$(name)_filtered_panama_filled.png"))
plot_H_panama(lon, lat, z, joinpath(@__DIR__, "H_$(name)_filtered_panama_filled_panama.png"))
ofile = joinpath(@__DIR__, @sprintf("etopo_%s_%s_%dkm_%dkm_panama_filled.jld2", name, res, L_land/1e3, L_ocean/1e3))
jldsave(ofile; lon, lat, z)
@info "Saved '$ofile'"