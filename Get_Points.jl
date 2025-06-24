using Gmsh: gmsh

gmsh.initialize()

gmsh.open("t10.msh")       

h = 0.1
function f(x,y)
    x*(1-x)*y*(1-y)*10
end 

nodeTags, coords, _ = gmsh.model.mesh.getNodes()
nodes = reshape(coords, 3, length(nodeTags))' 

length_of = size(nodes, 1)
tag_to_point_coords = Dict(nodeTags[i] => nodes[i,:] for i in 1:length_of) # generates a dictionary
# where each tag represents the coordinates of the points that the tags represent

elementType = gmsh.model.mesh.getElementType("Triangle", 1)
elementTags, elementNodeTags = gmsh.model.mesh.getElementsByType(elementType)
triangles = reshape(elementNodeTags, 3, :)' 

for i in 1:size(triangles, 1)  
    idx = triangles[i, :]   
    idx[:] = sort(Vector(idx))
    triangles[i, :] = idx
end


tag_to_nodes_triangles = Dict( elementTags[i] => triangles[i,:] for i in eachindex(elementTags))

edgeSets = Dict{Tuple{Int, Int}, Int}() 
nodeSets = Set{Int}()
loopTags = Int[]


function check_if_in_set(set, tag1, tag2)
    return (haskey(set, (tag1, tag2)))
end

gmsh.model.add("Reconstruct")

for i in 1:size(triangles, 1)  
    idx = triangles[i, :]   
    # idx[:] = sort(Vector(idx))
    # println(idx)
    triangle_coords = nodes[idx, 1:2]    
    for k in (1:3)
        check = idx[k]
        if !(check in nodeSets)
            push!(nodeSets, idx[k])
            x_cor, y_cor = triangle_coords[k, :]
            z_cor = f(x_cor, y_cor)
            gmsh.model.occ.addPoint(x_cor, y_cor, z_cor, h, idx[k])
        end
    end 

    # println(tag_to_nodes_triangles[elementTags[i]])
    node1 = idx[1]
    node2 = idx[2]
    node3 = idx[3]


    if !(check_if_in_set(edgeSets, node1, node2))
        edgeSets[(node1, node2)] = gmsh.model.occ.addLine(node1, node2)
    end
    if !(check_if_in_set(edgeSets, node2, node3))
        edgeSets[(node2, node3)] = gmsh.model.occ.addLine(node2, node3)
    end
    if !(check_if_in_set(edgeSets, node3, node1))
        edgeSets[(node3, node1)] = gmsh.model.occ.addLine(node3, node1)
    end 
    
end
println(nodeSets)

for i in (1:length(elementTags))
    points = tag_to_nodes_triangles[elementTags[i]] 
    point1 = points[1]
    point2 = points[2]
    point3 = points[3]

    edge1  = edgeSets[(point1, point2)]
    edge2  = edgeSets[(point2, point3)]
    edge3  = edgeSets[(point3, point1)]

    curve_loop = gmsh.model.occ.addCurveLoop([edge1, edge2, edge3])
    push!(loopTags, curve_loop)


end

volumeTags = Int[]
surfaceTags = Int[]

for i in (1:length(loopTags))
    surface = gmsh.model.occ.addPlaneSurface([loopTags[i]])
    push!(surfaceTags, surface)
end

c1=gmsh.model.occ.addPoint(0,0,0)
c2=gmsh.model.occ.addPoint(1,0,0)
c3=gmsh.model.occ.addPoint(1,1,0)
c4=gmsh.model.occ.addPoint(0,1,0)

l1=gmsh.model.occ.addLine(c1, c2)
l2=gmsh.model.occ.addLine(c2, c3)
l3=gmsh.model.occ.addLine(c3, c4)
l4=gmsh.model.occ.addLine(c4, c1)

curve1 = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
surface1 = gmsh.model.occ.addPlaneSurface([curve1])
push!(surfaceTags, surface1)
gmsh.model.occ.synchronize()
gmsh.model.occ.removeAllDuplicates()
shell = gmsh.model.occ.addSurfaceLoop(Int32.(surfaceTags)) 

vol   = gmsh.model.occ.addVolume([shell])        

gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(3)
gmsh.write("Reconstruct.msh")
gmsh.finalize()