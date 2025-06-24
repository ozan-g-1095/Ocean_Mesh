using Gmsh: gmsh

gmsh.initialize()
gmsh.model.add("toy_column")

# resolution
h = 1

# bottom triangle
gmsh.model.occ.addPoint(0.0, 0.0, -0.8, 1)
gmsh.model.occ.addPoint(1.0, 0.0, -1.0, 2)
gmsh.model.occ.addPoint(0.5, 0.5, -1.0, 3)
gmsh.model.occ.addPoint(0.0, 1.0, -1.0, 4)

# triangle one
gmsh.model.occ.addLine(1, 2, 1)
gmsh.model.occ.addLine(2, 3, 2)
gmsh.model.occ.addLine(3, 1, 3)
gmsh.model.occ.addCurveLoop([1, 2, 3], 1)
gmsh.model.occ.addPlaneSurface([1], 1)

# triangle two
# gmsh.model.occ.addLine(1, 3, 4) # this line was already made
gmsh.model.occ.addLine(3, 4, 5)
gmsh.model.occ.addLine(4, 1, 6)
gmsh.model.occ.addCurveLoop([-3, 5, 6], 2) # flip 3 to keep orientation
gmsh.model.occ.addPlaneSurface([2], 2)

# extrude triangle one to surface
gmsh.model.occ.extrude([(2, 1)], 0, 0, 1)

# extrude triangle two to surface
gmsh.model.occ.extrude([(2, 2)], 0, 0, 1)

# fuse the two extruded volumes
gmsh.model.occ.fuse([(3, 1)], [(3, 2)])

# box at surface
gmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3)

# cut the volume off at z = 0
gmsh.model.occ.cut([(3, 1)], [(3, 3)])

# set mesh resolution
gmsh.model.occ.mesh.setSize(gmsh.model.occ.getEntities(0), h)

# generate
gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(3)
gmsh.write("toy_column.msh")
gmsh.fltk.run()
gmsh.finalize()
