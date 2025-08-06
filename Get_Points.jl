using Gmsh: gmsh
using JLD2
using NCDatasets
using DataStructures

gmsh.initialize()

# ds = NCDataset("topo_25.1_coarsened1024.nc", "r")

# lon = ds["lon"][:]
# lat = ds["lat"][:]
# z   = ds["z"][:,:]
# close(ds)

d = jldopen("topo_data/etopo_blacksea_7min_50km.jld2")
lon = d["lon"]
lat = d["lat"]
z = d["z"]
close(d)

# gmsh.open("t10.msh")      
gmsh.open("Boundary_Fixed_Part2.msh") 

h_global = 0.5          
h_0    = 0.1       

R = 1 #6378000 # Earth radius in meters

α = 0.5 


function f(x,y)
    100 - (20*x)^2 - (30*y)^2 #((x-(3.14*10^6))/(10^3))^2 - ((y-(4.53*10^6))/(10^3))^2
end 

nodeTags, coords, _ = gmsh.model.mesh.getNodes()
nodes = reshape(coords, 3, length(nodeTags))' 

length_of = size(nodes, 1)
tag_to_point_coords = Dict(nodeTags[i] => nodes[i,:] for i in 1:length_of) # generates a dictionary
# where each tag represents the coordinates of the points that the tags represent

# println(tag_to_point_coords)

elementType = gmsh.model.mesh.getElementType("Triangle", 1)
elementTags, elementNodeTags = gmsh.model.mesh.getElementsByType(elementType)
triangles = reshape(elementNodeTags, 3, :)'

lineType = gmsh.model.mesh.getElementType("Line", 1)
lineTags, lineNodeTags = gmsh.model.mesh.getElementsByType(lineType)
lines = reshape(lineNodeTags, 2, :)'
tag_to_nodes_lines = Dict(lineTags[i] => lines[i,:] for i in eachindex(lineTags))

# curve_loops = gmsh.model.getBoundary([(2, surfaceTag)], oriented=false, recursive=false)
curve_entities = gmsh.model.getEntities(1)  # 1 = dimension for curves
curve_tags = [c[2] for c in curve_entities]

dim = 1  # dimension for curves
# tag = curve_tags[1]

line_ID_length = size(triangles, 1) * size(triangles, 2) 
line_tag_IDs = collect(1:line_ID_length)

free_line_tag_IDs = line_tag_IDs#setdiff(line_tag_IDs, curve_tags) 

# Returns: (lowerDimTags, higherDimTags)
# for tag in curve_tags
#     points, pairs = gmsh.model.getAdjacencies(dim, tag)
# end


A = gmsh.model.getEntities(2)

elementTypesTTriangle = Vector{Int}()
elementTagsTriangle = Vector{Vector{Int}}()
nodeTagsElementTriangle = Vector{Vector{Int}}()

nodes_of_triangles_in_each_loop = Vector{Vector{Int}}()

for elem in A
    elementTypesForLoop, elementTagsForLoop, nodeTagsElementForLoop = gmsh.model.mesh.getElements(elem[1], elem[2]) # elem = (dim, tag)
    push!(nodeTagsElementTriangle, nodeTagsElementForLoop[1])
    push!(elementTagsTriangle, elementTagsForLoop[1])
end



for i in (1:size(triangles, 1))
    idx = triangles[i, :]   
    idx[:] = sort(Vector(idx))
    triangles[i, :] = idx
end

tag_to_nodes_triangles = Dict( elementTags[i] => triangles[i,:] for i in eachindex(elementTags))

edgeSets = Dict{Tuple{Int, Int}, Int}() 
boundaryEdgeSet = Int[]
edges    = Set{Int}()  
nodeSets = Set{Int}()
mesh_to_occ_point = Dict{Int, Int}()
loopTags = Vector{Vector{Int}}()
boundary = Dict{Int, Tuple{Float64, Float64}}()

# @load "boundary_line_orders" ordered_surface_lines
@load "cell_coefficients" cell_coefficients

@load "x_min_normalize" x_min_normalize
@load "y_min_normalize" y_min_normalize
@load "x_max_normalize" x_max_normalize
@load "z_max" z_max
@load "L" L
@load "loops_fixed_boundary" loops_fixed_boundary

# @load "curve_loops" curve_loops

function check_if_in_set(set, tag1, tag2)
    if ((haskey(set, (tag1, tag2))) || (haskey(set, (tag2, tag1))))
        return true
    else
        return false
    end
end

function get_line_no(dict, tag1, tag2)
    for (k, val) in dict
        if (val[1] == tag1 && val[2] == tag2) || (val[1] == tag2 && val[2] == tag1)
            return k
        end
    end
    return 0 # error("Line connecting $tag1 and $tag2 not found in dictionary.")
end


function check_if_bonudary_edge(set, node1, node2)
    for i in (1:size(set, 1))
        if (set[i][2] == node1 && set[i][3] == node2) ||
           (set[i][2] == node2 && set[i][3] == node1)
            return set[i][1]
        end
    end
    return 0
end

function bilinear_depth(coef, x, y)
    return coef[1] + coef[2]*x + coef[3]*y + coef[4]*x*y
end


"""
    edges, boundary_indices, emap = all_edges(t)

Find all unique edges (no_edges x 2 array) in the triangulation `t` 
Second output is indices to the boundary edges.
Third output emap (nt x 3 array) is a mapping from local triangle edges
to the global edge list, i.e., emap[it,k] is the global edge number
for local edge k (1,2,3) in triangle it.
"""
function all_edges(t)
    etag = vcat(t[:,[1,2]], t[:,[2,3]], t[:,[3,1]])
    etag = hcat(sort(etag, dims=2), 1:3*size(t,1))
    etag = sortslices(etag, dims=1)
    dup = all(etag[2:end,1:2] - etag[1:end-1,1:2] .== 0, dims=2)[:]
    keep = .![false;dup]
    edges = etag[keep,1:2]
    emap = cumsum(keep)
    invpermute!(emap, etag[:,3])
    emap = reshape(emap,:,3)
    dup = [dup;false]
    dup = dup[keep]
    bndix = findall(.!dup)
    return edges, bndix, emap
end

"""
    e = boundary_nodes(t)

Find all boundary nodes in the triangulation `t`.
"""
function boundary_nodes(t)
    edges, boundary_indices, _ = all_edges(t)
    return unique(edges[boundary_indices,:][:])
end


function get_correct_edge(dict, tag1, tag2)
    if (haskey(dict, (tag1, tag2)))
        return dict[(tag1, tag2)]
    elseif (haskey(dict, (tag2, tag1)))
        return dict[(tag2, tag1)]
    end
end

unique_edges, bndix, emap = all_edges(triangles)
boundary_matrix = unique_edges[bndix, :]

each_node_and_neighbor_nodes = Dict{Int, Vector{Int}}()
unique_edge_boundary = Tuple{Int, Int}[]

function find_connectivities(unique_edge_sets, bndix)
    for elem in bndix
        # println(unique_edge_sets[elem])
        nodes = unique_edge_sets[elem, :]
        push!(unique_edge_boundary, (nodes[1], nodes[2]))
    end

    for line in unique_edge_boundary
        push!(get!(each_node_and_neighbor_nodes, line[1], Int[]), line[2])
        push!(get!(each_node_and_neighbor_nodes, line[2], Int[]), line[1])
    end
end


n = length(loops_fixed_boundary) 
loop_line_tags_ordered = [Int[] for _ in 1:n]

function create_boundary_loop(dict)
    for i in (1:length(loops_fixed_boundary))
        # println("dict: ", dict)
        k , v = first(dict)
        # println("This is k ", k)
        # println("This is v ", v)
        initial = k
        current = v[1]
        # println("initial ", initial)
        # println("current ", current)
        # delete!(dict, k)

        line_no = get_correct_edge(boundary_node_IDs_line_ID, initial, current)
        # println("line_no: ", line_no)
        edge_new = gmsh.model.occ.addLine(initial, current, line_no)
        push!(loop_line_tags_ordered[i], edge_new)
        push!(edges, edge_new)
        edgeSets[(initial, current)] = edge_new
        # for i in 1:length(dict)
        while current != initial
            # println("current ", current)
            prev = current
            v = dict[prev]
            # println("prev ", prev, " new ", v)
            if v[1] == k
                current = v[2]
            else
                current = v[1]
            end
            k = prev
            # println("current to be deleted ", current)
            delete!(dict, prev)

            line_no = get_correct_edge(boundary_node_IDs_line_ID, prev, current)
            # println(" LINE NO ", line_no)
            edge_new = gmsh.model.occ.addLine(prev, current, line_no)
            push!(loop_line_tags_ordered[i], edge_new)
            push!(edges, edge_new)
            edgeSets[(prev, current)] = edge_new

            # if current == initial 
            #     # println("DONE HERE")
            #     break
            # end

        end
        # println("dict: ", dict)
        delete!(dict, current)
    end
end

gmsh.finalize()


gmsh.initialize()
gmsh.model.add("Reconstruct")

number_oooo = 0

pointtagsno = SortedSet()

boundary_node_IDs_line_ID = Dict{Tuple{Int, Int}, Int}()

function add_nodes_and_lines()
    for i in (1:size(triangles, 1))
        # println(i)
        idx = triangles[i, :] 

        # println("idx ", idx)

        node1 = idx[1]
        node2 = idx[2]
        node3 = idx[3]

        triangle_coords = zeros(3, 3)
        node1_coords = tag_to_point_coords[node1]  
        triangle_coords[1, :] = node1_coords
        node2_coords = tag_to_point_coords[node2]
        triangle_coords[2, :] = node2_coords
        node3_coords = tag_to_point_coords[node3]
        triangle_coords[3, :] = node3_coords

        # println("triangle coords ", node1_coords, " ", node2_coords, " ", node3_coords)

        for k in (1:3)
            check = idx[k]
            if !(check in nodeSets)
                # println(nodeSets, " this is the node set ")
                push!(nodeSets, idx[k])
                x_cor, y_cor = triangle_coords[k, 1:2]

                lon_pos = rad2deg(((x_cor *(x_max_normalize-x_min_normalize) ) + x_min_normalize ) / R)# rad2deg(x_cor/R)
                lat_pos = rad2deg(((y_cor *(x_max_normalize-x_min_normalize) ) + y_min_normalize ) / R)#rad2deg(y_cor/R)
                
                p = argmin(abs.(lon .- lon_pos))
                j = argmin(abs.(lat .- lat_pos))

                # println("lat_pos: ", lat_pos, " lon_pos: ", lon_pos)
                # println("p: ", p, " j: ", j)

                z_cor_calc = bilinear_depth(cell_coefficients[p, j], x_cor, y_cor)
                # println("cell_coefficients : ", cell_coefficients[p,j])
                # println("Calculated z depth ", z_cor_calc)

                if idx[k] in boundary_matrix
                    z_cor_calc = 0.0
                end

                if -0.01< z_cor_calc< 0
                    z_cor_calc -= 0.01
                elseif z_cor_calc > 0 
                    z_cor_calc = -0.01/5
                    # println("zcor ", z_cor_calc)
                # elseif -0.05 < z_cor_calc < 0 
                #     println("zcor very sMAAAAL", z_cor_calc)
                #     z_cor_calc = -rand(0.05:0.01:0.15)
                end

                push!(pointtagsno, idx[k])
                gmsh.model.occ.addPoint(x_cor, y_cor, z_cor_calc, h_0/5, idx[k])
                # gmsh.model.occ.addPoint((x_cor-x_min_normalize)/(x_max_normalize - x_min_normalize), (y_cor-y_min_normalize)/(x_max_normalize - x_min_normalize), z_cor_calc, h_0, idx[k])

            end
        end 

            if !( (node1, node2) in unique_edge_boundary ||  (node2, node1) in unique_edge_boundary)
                if !(check_if_in_set(edgeSets, node1, node2)) 
                    edge_new = gmsh.model.occ.addLine(node1, node2, popfirst!(free_line_tag_IDs) )
                    # println("Lets see what edge new is ", edge_new)
                    push!(edges, edge_new)
                    edgeSets[(node1, node2)]  = edge_new
                end
            else  
                boundary_node_IDs_line_ID[(node1, node2)] = popfirst!(free_line_tag_IDs)
            end

            if !( (node2, node3) in unique_edge_boundary ||  (node3, node2) in unique_edge_boundary)
                if !(check_if_in_set(edgeSets, node2, node3)) 
                    edge_new = gmsh.model.occ.addLine(node2, node3, popfirst!(free_line_tag_IDs) )
                    # println("Lets see what edge new is ", edge_new)
                    push!(edges, edge_new)
                    edgeSets[(node2, node3)]  = edge_new
                end
            else  
                boundary_node_IDs_line_ID[(node2, node3)] = popfirst!(free_line_tag_IDs)
            end

            if !( (node3, node1) in unique_edge_boundary ||  (node1, node3) in unique_edge_boundary)
                if !(check_if_in_set(edgeSets, node3, node1)) 
                    edge_new = gmsh.model.occ.addLine(node3, node1, popfirst!(free_line_tag_IDs) )
                    # println("Lets see what edge new is ", edge_new)
                    push!(edges, edge_new)
                    edgeSets[(node3, node1)]  = edge_new
                end
            else  
                boundary_node_IDs_line_ID[(node3, node1)] = popfirst!(free_line_tag_IDs)
            end
    end
end


find_connectivities(unique_edges, bndix)
add_nodes_and_lines()


create_boundary_loop(each_node_and_neighbor_nodes)
collect(pointtagsno)


for idx in (1:length(elementTagsTriangle))

    elementTagsForLoop = elementTagsTriangle[idx]

    push!(loopTags, Int[])

    for i in (1:length(elementTagsForLoop))
        points = tag_to_nodes_triangles[elementTagsForLoop[i]]
        # println("Here are the element tags number ", i, " ",  Int.(points) )

        point1 = points[1]
        point2 = points[2]
        point3 = points[3]

        edge1  = get_correct_edge(edgeSets, point1, point2)
        edge2  = get_correct_edge(edgeSets, point2, point3)
        edge3  = get_correct_edge(edgeSets, point3, point1)

        # println("edge1 ", edge1, " edge2 ", edge2, " edge3 ", edge3)

        curve_loop = gmsh.model.occ.addCurveLoop([edge1, edge2, edge3])
        push!(loopTags[idx], curve_loop)
    end

end

n = length(loop_line_tags_ordered) 
surfaceTags = Vector{Vector{Int}}() # Int[]
volumeTags = Vector{Vector{Int}}()

for i in (1:length(loopTags))
    surfaces_of_loop = Int[]
    for elem in loopTags[i]
        # println("elem: ", elem)
        surface = gmsh.model.occ.addPlaneSurface([elem])
        push!(surfaces_of_loop, surface)
    end
    push!(surfaceTags, surfaces_of_loop)
end

println("Hello it has come here")

top_surfaces = Int[]

println("HERE??")

for loop in loop_line_tags_ordered
    curve1 = gmsh.model.occ.addCurveLoop(loop)
    surface1 = gmsh.model.occ.addPlaneSurface([curve1])
    push!(top_surfaces, surface1)
end

for (idx, elem) in enumerate(top_surfaces)
    # println("idx: ", idx)
    # println("surfaceTags[idx]: ", surfaceTags[idx])
    push!(surfaceTags[idx], elem)
end

volumes = Int[]
for surfaceTagsLoop in surfaceTags
    shell = gmsh.model.occ.addSurfaceLoop(Int.(surfaceTagsLoop)) # [surface1] 
    vol   = gmsh.model.occ.addVolume([shell])   
    push!(volumes, vol)  
end

gmsh.model.occ.removeAllDuplicates()
gmsh.model.occ.synchronize()


for i in 1:length(top_surfaces)
    pop!(surfaceTags[i])
end

gmsh.model.addPhysicalGroup(0, Int32.(collect(nodeSets)), 1, "bot")
gmsh.model.addPhysicalGroup(1, Int32.(edges), 2, "bot") #Int32.(collect(edges))
gmsh.model.addPhysicalGroup(2, top_surfaces, 3, "sfc")
gmsh.model.addPhysicalGroup(2, vcat(surfaceTags...), 4, "bot")
gmsh.model.addPhysicalGroup(3, volumes, 5, "int")

gmsh.model.mesh.generate(3)
gmsh.write("Reconstruct.msh")
gmsh.write("Reconstruct.vtk")
gmsh.finalize()