using NCDatasets
using Gmsh: gmsh
using Shapefile
using Plots

gmsh.initialize()
gmsh.model.add("Caspian_Sea")

R = 6378000

shp = Shapefile.Table("lines.shp")

geoms = Shapefile.shapes(shp)
# println(geoms)

# x = Float32[]
# y = Float32[]

# for i in (1:length(geoms))
#     point = geoms[i].points
#     x_cor = point[1].x
#     push!(x, x_cor)
#     y_cor = point[1].y
#     push!(y, y_cor)
# end

# # scatter(x, y, xlabel="Longitude", ylabel="Latitude", title="Point Scatter")

# for i in (1:length(x))
#     lat_point = deg2rad(y[i])
#     long_point = deg2rad(x[i])

#     x_cor = R * cos(lat_point) * cos(long_point)
#     y_cor = R * cos(lat_point) * sin(long_point)
#     z_cor = R * sin(lat_point)
#     gmsh.model.occ.addPoint(x_cor, y_cor, z_cor)

# end




# for feature in shp
#     geoms = Shapefile.shapes(shp)
#     println("whaaa")
#     if geoms isa Shapefile.Point
#         println(geoms.x, ", ", geoms.y)
#     end
# end

# for feature in shp
#     println(feature.getproperty)
#     # println(feature.geometry)
# end



ds = NCDataset("topo_25.1_coarsened1024.nc", "r")

lon = ds["lon"][:]
lat = ds["lat"][:]
z   = ds["z"][:,:]
# z[z .> 0] .= 0


for i in (argmin(abs.(lon .- 15)):  argmin(abs.(lon .- 45)) ) # longitude
    for j in (argmin(abs.(lat .- 30)):  argmin(abs.(lat .- 45))) # latitude
        depth = z[i,j]
    
        if (abs(depth) < 2)
            println("depth ", z[i,j])
            
            lat_point = deg2rad(lat[j])
            long_point = deg2rad(lon[i])

            x_cor = R * cos(lat_point) * cos(long_point)
            y_cor = R * cos(lat_point) * sin(long_point)
            z_cor = R * sin(lat_point)
            gmsh.model.occ.addPoint(x_cor, y_cor, z_cor)
        end
    end
end


gmsh.model.occ.synchronize()

gmsh.model.mesh.generate(0)
gmsh.write("Caspian_Sea.msh")
gmsh.finalize()