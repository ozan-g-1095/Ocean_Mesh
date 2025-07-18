using Gmsh: gmsh
using Roots

h = 1 # resolution

gmsh.initialize()
gmsh.model.add("Concave2D")

using Roots
point_no = 53


function Find_Points()
    # Define constants
    Px = 5.0
    Py = 25
    A = 1

    # Define the equation as a function of t
    f(t) = t^2 - 2t*Px + t^4 - 2t^2*Py - (A^2 - Px^2 - Py^2)

    for i in 1:point_no
        gmsh.model.geo.addPoint(Px, Py, 0, h)

        roots = find_zeros(f, Px-10, Px)

        Px = roots[1]
        Py = Px^2

        println("Roots of the equation:", roots)
    end 

end 

Find_Points()

# connect as lines
for i in 1:point_no-1
    gmsh.model.geo.addLine(i, i+1)
end
gmsh.model.geo.addLine(point_no, 1)

# loop together to define region
gmsh.model.geo.addCurveLoop(1:point_no)
gmsh.model.geo.addPlaneSurface([1], 1)
# gmsh.model.geo.synchronize()
# extrude the surface to create a volume (for fun)
gmsh.model.geo.extrude([(2, 1)], 0, 0, 15)

# # # generate
gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(3)
#  # gmsh.model.mesh.generate(3) # if you want to do the extrusion
gmsh.write("Concave2D.msh")
gmsh.finalize()