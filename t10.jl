# ------------------------------------------------------------------------------
#
#  Gmsh Julia tutorial 10
#
#  Mesh size fields
#
# ------------------------------------------------------------------------------

# In addition to specifying target mesh sizes at the points of the geometry (see
# `t1.jl') or using a background mesh (see `t7.jl'), you can use general mesh
# size "Fields".

using Gmsh: gmsh


gmsh.initialize(append!(["gmsh"], ARGS))

gmsh.model.add("t10")

# Let's create a simple rectangular geometry:

r = 1  # radius

# Define points

disc = gmsh.model.occ.addDisk(0,0,0,r,r)

# p0 = gmsh.model.occ.addPoint(0.0, 0.0, 0.0)  
# p1 = gmsh.model.occ.addPoint(r, 0.0, 0.0)   
# p2 = gmsh.model.occ.addPoint(0.0, r, 0.0)    
# p3 = gmsh.model.occ.addPoint(-r, 0.0, 0.0)   
# p4 = gmsh.model.occ.addPoint(0.0, -r, 0.0)   


# c1 = gmsh.model.occ.addCircleArc(p1, p0, p2)
# c2 = gmsh.model.occ.addCircleArc(p2, p0, p3)
# c3 = gmsh.model.occ.addCircleArc(p3, p0, p4)
# c4 = gmsh.model.occ.addCircleArc(p4, p0, p1)

# lc = .15
# gmsh.model.geo.addPoint(0.0, 0.0, 0, lc, 1)
# gmsh.model.geo.addPoint(1, 0.0, 0, lc, 2)
# gmsh.model.geo.addPoint(1, 1, 0, lc, 3)
# gmsh.model.geo.addPoint(0, 1, 0, lc, 4)
# # gmsh.model.geo.addPoint(0.2, .5, 0, lc, 5)

# gmsh.model.geo.addLine(1, 2, 1)
# gmsh.model.geo.addLine(2, 3, 2)
# gmsh.model.geo.addLine(3, 4, 3)
# gmsh.model.geo.addLine(4, 1, 4)

# gmsh.model.occ.addCurveLoop([c1, c2, c3, c4], 5)
# gmsh.model.occ.addPlaneSurface([5], 6)

gmsh.model.occ.synchronize()

# Say we would like to obtain mesh elements with size lc/30 near curve 2 and
# point 5, and size lc elsewhere. To achieve this, we can use two fields:
# "Distance", and "Threshold". We first define a Distance field (`Field[1]') on
# points 5 and on curve 2. This field returns the distance to point 5 and to
# (100 equidistant points on) curve 2.
# gmsh.model.mesh.field.add("Distance", 1)
# gmsh.model.mesh.field.setNumbers(1, "PointsList", [5])
# gmsh.model.mesh.field.setNumbers(1, "CurvesList", [2])
# gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

# We then define a `Threshold' field, which uses the return value of the
# `Distance' field 1 in order to define a simple change in element size
# depending on the computed distances
#
# SizeMax -                     /------------------
#                              /
#                             /
#                            /
# SizeMin -o----------------/
#          |                |    |
#        Point         DistMin  DistMax
# gmsh.model.mesh.field.add("Threshold", 2)
# gmsh.model.mesh.field.setNumber(2, "InField", 1)
# gmsh.model.mesh.field.setNumber(2, "SizeMin", lc / 30)
# gmsh.model.mesh.field.setNumber(2, "SizeMax", lc)
# gmsh.model.mesh.field.setNumber(2, "DistMin", 0.15)
# gmsh.model.mesh.field.setNumber(2, "DistMax", 0.5)

# Say we want to modulate the mesh element sizes using a mathematical function
# of the spatial coordinates. We can do this with the MathEval field:

h = 0.1

# find_next_point(0,0)

gmsh.model.mesh.field.add("MathEval", 3)
# gmsh.model.mesh.field.setString(3, "F", "0.1 * ((x +1)/2 ) ")
# gmsh.model.mesh.field.setString(3, "F", "$h/ ( ( 1 + ( (1-2*x)*y*(1-y)*10 )^2 + ( (1-2*y)*x*(1-x)*10 )^2 )  ) ^ (1/2)" ) this is for the square basin model
gmsh.model.mesh.field.setString(3, "F", "($h * 3 )/ ( ( 1 + ( 2*x )^2 + ( (2*y)^2 )  ) ^ (1/2) )") # for circular basin
# gmsh.model.mesh.field.setString(3, "F", "$h/ ( ( 1 + ( (3*x)^2 - (9*x)^8 + (7*x)^6   )^2 + ( (0)^2 )  ) ^ (1/2) )")

# We could also combine MathEval with values coming from other fields. For
# example, let's define a `Distance' field around point 1
# gmsh.model.mesh.field.add("Distance", 4)
# gmsh.model.mesh.field.setNumbers(4, "PointsList", [1])

# We can then create a `MathEval' field with a function that depends on the
# return value of the `Distance' field 4, i.e., depending on the distance to
# point 1 (here using a cubic law, with minimum element size = lc / 100)
# gmsh.model.mesh.field.add("MathEval", 5)
# gmsh.model.mesh.field.setString(5, "F", "F4^3 + " * string(lc / 100))

# We could also use a `Box' field to impose a step change in element sizes
# inside a box
# gmsh.model.mesh.field.add("Box", 6)
# gmsh.model.mesh.field.setNumber(6, "VIn", lc / 15)
# gmsh.model.mesh.field.setNumber(6, "VOut", lc)
# gmsh.model.mesh.field.setNumber(6, "XMin", 0.3)
# gmsh.model.mesh.field.setNumber(6, "XMax", 0.6)
# gmsh.model.mesh.field.setNumber(6, "YMin", 0.3)
# gmsh.model.mesh.field.setNumber(6, "YMax", 0.6)

# Many other types of fields are available: see the reference manual for a
# complete list. You can also create fields directly in the graphical user
# interface by selecting `Define->Size fields' in the `Mesh' module.

# Let's use the minimum of all the fields as the mesh size field:
# gmsh.model.mesh.field.add("Min", 7)
# gmsh.model.mesh.field.setNumbers(7, "FieldsList", [2, 3, 5, 6])

gmsh.model.mesh.field.setAsBackgroundMesh(3)

# The API also allows to set a global mesh size callback, which is called each
# time the mesh size is queried
# function meshSizeCallback(dim, tag, x, y, z, lc)
#     return min(lc, 0.02 * x + 0.01)
# end

# gmsh.model.mesh.setSizeCallback(meshSizeCallback)

# To determine the size of mesh elements, Gmsh locally computes the minimum of
#
# 1) the size of the model bounding box;
# 2) if `Mesh.MeshSizeFromPoints' is set, the mesh size specified at geometrical
#    points;
# 3) if `Mesh.MeshSizeFromCurvature' is positive, the mesh size based on
#    curvature (the value specifying the number of elements per 2 * pi rad);
# 4) the background mesh size field;
# 5) any per-entity mesh size constraint;
#
# The value can then be further modified by the mesh size callback, if any,
# before being constrained in the interval [`Mesh.MeshSizeMin',
# `Mesh.MeshSizeMax'] and multiplied by `Mesh.MeshSizeFactor'.  In addition,
# boundary mesh sizes are interpolated inside surfaces and/or volumes depending
# on the value of `Mesh.MeshSizeExtendFromBoundary' (which is set by default).
#
# When the element size is fully specified by a mesh size field (as it is in
# this example), it is thus often desirable to set

gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

# This will prevent over-refinement due to small mesh sizes on the boundary.

# Finally, while the default "Frontal-Delaunay" 2D meshing algorithm
# (Mesh.Algorithm = 6) usually leads to the highest quality meshes, the
# "Delaunay" algorithm (Mesh.Algorithm = 5) will handle complex mesh size fields
# better - in particular size fields with large element size gradients:




gmsh.option.setNumber("Mesh.Algorithm", 5)

# points = gmsh.model.getEntities(0)


# println(points)

gmsh.model.mesh.generate(2)
gmsh.write("t10.msh")

# Launch the GUI to see the results:
# if !("-nopopup" in ARGS)
#     gmsh.fltk.run()
# end

element_type, element_tags, node_Tags =gmsh.model.mesh.getElements(2)
# println("the elements", convert(node_Tags, Int))
gmsh.finalize()