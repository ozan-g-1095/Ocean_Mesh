using Gmsh: gmsh
using Roots
using Plots

gmsh.initialize()
gmsh.model.add("Concave2_2D")

point_no = 30
h = 0.1

function f_prime(x)
    2x 
end 

function f(x)
    x^2
end 

points = Tuple{Float64, Float64}[]

function find_next_point(xi)
    for i in 1:point_no
        push!(points,(xi,0))
        yf = f(xi)
        gmsh.model.geo.addPoint(xi, yf, 0, h)
        xi = h/ (sqrt( 1 + (f_prime(xi))^2 )) + xi
        
        println("Here are the different x points", xi , " and their y: ", yf)

    end 
end

find_next_point(-1)

# println(x_points)

xs = getindex.(points, 1)
ys = getindex.(points, 2)

# Plot the points
plot(xs, ys, marker=:circle, label="", xlabel="x", ylabel="y", title="Tuple Points")


# # connect as lines
# for i in 1:point_no-1
#     gmsh.model.geo.addLine(i, i+1)
# end
# gmsh.model.geo.addLine(point_no, 1)

# # loop together to define region
# gmsh.model.geo.addCurveLoop(1:point_no)
# gmsh.model.geo.addPlaneSurface([1], 1)
# # gmsh.model.geo.synchronize()
# # extrude the surface to create a volume (for fun)
# gmsh.model.geo.extrude([(2, 1)], 0, 0, 15)

# # # # generate
# gmsh.model.geo.synchronize()
# gmsh.model.mesh.generate(3)
# #  # gmsh.model.mesh.generate(3) # if you want to do the extrusion
# gmsh.write("Concave2_2D.msh")
# gmsh.finalize()

