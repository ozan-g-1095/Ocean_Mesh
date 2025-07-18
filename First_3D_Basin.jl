using Gmsh: gmsh
using Roots
using Plots

gmsh.initialize()
gmsh.model.add("Basin_3D")

point_no = 30
h = 0.1

function f_prime(x,y)
    ((1-2x)*y*(1-y)*10, (1-2y)*x*(1-x)*10)
end 

function f(x,y)
    x*(1-x)*y*(1-y)*10
end 

points = Tuple{Float64, Float64}[]

function find_next_point(xi, yi)
    while (1-yi > h)
        if (1-xi > h) 
            push!(points,(xi,yi))
            depth = f(xi, yi)
            gmsh.model.geo.addPoint(xi, yi, f(xi,yi), h)
            xi = h/ (sqrt( 1 + (f_prime(xi,yi)[1])^2 ))+ xi
        end
            println("Here are the different x points", xi , " and their y: ", yi)
        if (1-xi < h)
            yi = yi +h
            xi = 0 
        end
    end 
end

find_next_point(0,0)

# println(x_points)

xs = getindex.(points, 1)
ys = getindex.(points, 2)

# Plot the points
plot(xs, ys, marker=:circle, label="", xlabel="x", ylabel="y", title="Tuple Points")


# connect as lines
for i in 1:length(points)-1
    gmsh.model.geo.addLine(i, i+1)
end
gmsh.model.geo.addLine(length(points), 1)

# # loop together to define region
# gmsh.model.geo.addCurveLoop(1:length(points))
# gmsh.model.geo.addPlaneSurface([1], 1)
# gmsh.model.geo.synchronize()
# extrude the surface to create a volume (for fun)
# gmsh.model.geo.extrude([(2, 1)], 0, 0, 15)

# # # generate
gmsh.model.geo.synchronize()
# gmsh.model.mesh.generate(3)
#  # gmsh.model.mesh.generate(3) # if you want to do the extrusion
gmsh.write("Basin_3D.msh")
gmsh.finalize()

