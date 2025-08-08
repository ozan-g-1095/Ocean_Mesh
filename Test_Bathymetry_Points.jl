using Interpolations
using Plots
gr()

using NCDatasets

using ProgressMeter
using Interpolations
using JLD2

using Contour
using LinearAlgebra

using Gmsh: gmsh
R = 6378000
gmsh.initialize()
gmsh.model.add("Test_Bathymetry_Points")

ds = NCDataset("topo_data/ETOPO_2022_v1_60s_N90W180_surface.nc", "r")

lon = ds["lon"][:]
lat = ds["lat"][:]
z   = ds["z"][:,:]
close(ds)

# d = jldopen("topo_data/etopo_blacksea_5min_50km.jld2")
# lon = d["lon"]
# lat = d["lat"]
# z = d["z"]
# close(d)


lon_min =26.27#-179
lon_max =42.25 #179

lat_min= 40.13#-79
lat_max= 47.70#79

displacement_lon = abs(lon[1] - lon[2])
displacement_lat = abs(lat[1] - lat[2])

min_dist = sqrt((displacement_lat/4)^2 + (displacement_lon/4)^2 )


function find_biliniear_interpolation()
    @showprogress 1 ("Processing longitude loop...") for i in (argmin(abs.(lon .- lon_min)): argmin(abs.(lon .- lon_max))) # longitude
        for j in (argmin(abs.(lat .- lat_min)):argmin(abs.(lat .- lat_max))) # latitude
            if z[i,j]/50<=0
                x_max = abs((argmin(abs.(lon .- lon_min)) - argmin(abs.(lon .- lon_max))))
                y_max = abs((argmin(abs.(lat .- lat_min)) -argmin(abs.(lat .- lat_max))))
                gmsh.model.occ.addPoint(i/x_max, j/x_max, z[i,j]/x_max)
            end
        end
    end
    gmsh.model.occ.synchronize()
end


find_biliniear_interpolation()

xmin, xmax = 12.9, 13.8
ymin, ymax = 8.15, 8.65
lc = 0.0              # mesh size at points (0 lets background settings control)

# addRectangle(x, y, z, dx, dy; tag?)
surf_tag = gmsh.model.occ.addRectangle(xmin, ymin, 0.0, xmax - xmin, ymax - ymin)

gmsh.model.occ.synchronize()





# find_next_point(0,0)

# gmsh.model.mesh.field.add("MathEval", 3)
# gmsh.model.mesh.field.setString(3, "F", "($h*3)") # for circular basin / ( ( 1 + ( 2*x )^2 + ( (2*y)^2 )  ) ^ (1/2) )

# gmsh.model.mesh.field.setAsBackgroundMesh(3)


# gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
# gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
# gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

# gmsh.option.setNumber("Mesh.Algorithm", 5)


gmsh.model.occ.removeAllDuplicates()
# gmsh.option.setNumber("Mesh.SaveAll", 1)

gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(3)
gmsh.write("Test_Bathymetry_Points.msh")
gmsh.finalize()