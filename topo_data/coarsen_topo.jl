using NCDatasets
using Plots
using LaTeXStrings
using JLD2

# data from https://www.ncei.noaa.gov/products/etopo-global-relief-model
# 60 arc-second resolution ice surface elevation

d = NCDataset("ETOPO_2022_v1_60s_N90W180_surface.nc")
lon = d["lon"][:]
lat = d["lat"][:]
z = Float64.(d["z"][:, :])
close(d)

# d = NCDataset("ETOPO_2022_v1_60s_N90W180_geoid.nc")
# lon_geoid = d["lon"][:]
# lat_geoid = d["lat"][:]
# z_geoid = Float64.(d["z"][:, :])
# close(d)

# @assert lon_geoid == lon
# @assert lat_geoid == lat

# # convert EGM2008 to WGS84
# z .+= z_geoid

"""
    lon, lat, z = coarsen_topo(lon, lat, z; arcmins=60)

Coarse topography data to a resolution of `arcmins` arc-minutes (default 60).
"""
function coarsen_topo(lon, lat, z; arcmins=60)
    lon = lon[1:arcmins:end]
    lat = lat[1:arcmins:end]
    z = z[1:arcmins:end, 1:arcmins:end]
    return lon, lat, z
end

function plot_H(lon, lat, z, fname; vmax=6000)
    H = -z
    H[H .< 0] .= NaN
    heatmap(lon, lat, H', dpi=300,
            xlabel="Longitude (°)", ylabel="Latitude (°)",
            cb_title=L"$H$ (m)", aspect_ratio=:equal,
            clims=(0, vmax), c=:seaborn_rocket_gradient,
            grid=false, xlimits=(lon[1], lon[end]),
            ylimits=(lat[1], lat[end]), tickdirection=:out)
    savefig(fname)
    @info "Saved '$fname'"
end

# # coarsen global data
# arcmins = 60 # 60 arc-minutes (1 degree)
# lon_global, lat_global, z_global = coarsen_topo(lon, lat, z; arcmins)
# plot_H(lon_global, lat_global, z_global, "H.png")
# jldsave("etopo_global_$(arcmins)min.jld2"; lon=lon_global, lat=lat_global, z=z_global)
# @info "Saved 'etopo_global_$(arcmins)min.jld2'"

# # coarsen and zoom in to Black Sea region 
# arcmins = 5 # 5 arc-minutes (1/12th degree)
# lon_box, lat_box, z_box = coarsen_topo(lon, lat, z; arcmins=5)
# lon_min, lon_max = 26.2, 42.25
# lat_min, lat_max = 40.13, 47.71
# lon_mask = lon_min .≤ lon_box .≤ lon_max
# lat_mask = lat_min .≤ lat_box .≤ lat_max
# lon_box = lon_box[lon_mask]
# lat_box = lat_box[lat_mask]
# z_box = z_box[lon_mask, lat_mask]
# plot_H(lon_box, lat_box, z_box, "H_blacksea.png")
# jldsave("etopo_blacksea_$(arcmins)min.jld2"; lon=lon_box, lat=lat_box, z=z_box)
# @info "Saved 'etopo_blacksea_$(arcmins)min.jld2'"