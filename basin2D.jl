using Gmsh: gmsh
using Roots

h = 0.1 # resolution

gmsh.initialize()
gmsh.model.add("basin2D")

# surface "mesh"
# x = -1:0.0:1
nx = 45 #length(x)

function Find_Points()
    # Define parameters
    A  = 0.1
    Px = 1
    Py = 0
    C = (A^2 - (Px^2) - (Py^2) - 1) / 2

    # Define the function
    f(x) = cos(x) * (Px) + sin(x) * (Py) - (C)
    # # depth function
    # H = 1 .- x.^2

    # add points at z = -H
    for i in 1:nx
        gmsh.model.geo.addPoint(Px, Py, 0, h)

        println(" x position: ", Px, " y position: ", Py, "index ", i)
        x_sol = find_zeros(f, 0, 4pi)

        if (Px == 1 && Py == 0)
            Px = -cos(x_sol[1])
            Py = sin(x_sol[1])
        else
            Px = -cos(x_sol[2])
            Py = -sin(x_sol[2])
        end

    end
end

Find_Points()

# connect as lines
for i in 1:nx-1
    gmsh.model.geo.addLine(i, i+1)
end
gmsh.model.geo.addLine(nx, 1)

# loop together to define region
gmsh.model.geo.addCurveLoop(1:nx)
gmsh.model.geo.addPlaneSurface([1], 1)
# gmsh.model.geo.synchronize()
# # # # extrude the surface to create a volume (for fun)
# # # gmsh.model.geo.extrude([(2, 1)], 0, 1, 0)

# # # generate
gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)
#  # gmsh.model.mesh.generate(3) # if you want to do the extrusion
gmsh.write("basin2D.msh")
gmsh.finalize()