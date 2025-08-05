include("utils.jl")

using Gmsh: gmsh
using JLD2

gmsh.initialize()
gmsh.open("Boundary_Fixed.msh")

h_0 = 0.1

elementType = gmsh.model.mesh.getElementType("Triangle", 1)
elementTags, elementNodeTags = gmsh.model.mesh.getElementsByType(elementType)
triangles = reshape(elementNodeTags, 3, :)'


# entities = gmsh.model.get_entities()
# dim = entities[end][1]
# tag = entities[end][2]
# _, _, elements = gmsh.model.mesh.get_elements(dim, tag)
# elements = Matrix(reshape(Int.(elements[1]), 3, :)')

# nodeTags, coords, _ = gmsh.model.mesh.getNodes()
# nodes = reshape(coords, 3, length(nodeTags))' 

# nodeTags, coords, _ = gmsh.model.mesh.getNodes()


@load "ordered_surface_lines" ordered_surface_lines

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

function boundary_nodes(t)
    edges, boundary_indices, _ = all_edges(t)
    return unique(edges[boundary_indices,:][:])
end


function get_line_tag_from_points(p1::Int, p2::Int)
    dim = 1  # Lines are 1D entities
    tags = gmsh.model.getEntities(dim)
    println("tags ", tags)
    for (dim_tag, line_tag) in tags
        node_tags, _ = gmsh.model.getAdjacencies(dim, line_tag)
        if length(node_tags) == 2
            if (p1 == node_tags[1] && p2 == node_tags[2]) || (p1 == node_tags[2] && p2 == node_tags[1])
                return line_tag
            end
        end
    end
    error("No line found connecting points $p1 and $p2")
end


# line = get_line_tag_from_points(510, 515)

# gmsh.model.occ.remove([(0, 200)], true)

# elementTypes, elementTags, nodeTags_for_lines = gmsh.model.mesh.getElements(1)

unique_edges, bndix, emap = all_edges(triangles)

edges_on_the_boundary = Tuple{Int, Int, Bool}[] 

for elem in bndix
    push!(edges_on_the_boundary, (unique_edges[elem, :][1],unique_edges[elem, :][2], false) )
end

edges_on_the_boundary_ordered = Vector{Vector{NTuple{4, Int64}}}()


function order_the_surface_lines(lines_on_surface, surface_lines_ordered)
    line_numbers = collect(1:length(lines_on_surface)*2)
    global loop_no = 1
    while !isempty(lines_on_surface)
        initial = lines_on_surface[1][1]
        current = lines_on_surface[1][2]
        current_index = 1

        line_no = popfirst!(line_numbers)
        push!(surface_lines_ordered, [(Int64(line_no), initial, current, lines_on_surface[1][3])])
        println((Int64(line_no), initial, current))

        while (initial != current)
            global loop_formed = false
            
            deleteat!(lines_on_surface, current_index) # it deletes the chosen point so that the algorithm can march forward and find the 
            # next grid that shares the chosen current's tag. 
            has_found_connectivity = false
            for (idx,elem) in enumerate(lines_on_surface)
               
                if elem[1] == current
                    println("Helllloooo")
                    
                    line_no = popfirst!(line_numbers) #gmsh.model.occ.addLine(elem[1], elem[2])
                    
                    push!(surface_lines_ordered[loop_no], (Int64(line_no), elem[1], elem[2], elem[3]))
                
                    current = elem[2]
                    current_index = idx
                    has_found_connectivity = true
                    break
                    
                elseif elem[2] == current
                    println("Helllloooo")
                    line_no = popfirst!(line_numbers)# gmsh.model.occ.addLine(elem[2], elem[1])
                    
                    push!(surface_lines_ordered[loop_no], (Int64(line_no), elem[2], elem[1], elem[3]))
                
                    current = elem[1]
                    current_index = idx
                    has_found_connectivity = true  
                    break
                end
            end

            if !(has_found_connectivity)
                for i in (length(surface_lines_ordered[loop_no]):-1:1)
                    # gmsh.model.occ.remove(([1,surface_lines_ordered[loop_no][i][1]])) # recursive = true
                end
                deleteat!(surface_lines_ordered, loop_no)
                break
            else
                loop_formed = true
            end

        end
        if (loop_formed)
            global loop_no += 1
        end
    end
end


order_the_surface_lines(edges_on_the_boundary, edges_on_the_boundary_ordered)

# nodes_on_the_boundary = boundary_nodes(triangles)
# boundary_edges = Tuple{Int, Int, Bool}[]

# println(gmsh.model.mesh.getNodes())

nodeTags, coords, _ = gmsh.model.mesh.getNodes()
nodes = reshape(coords, 3, length(nodeTags))' 

length_of = size(nodes, 1)
tag_to_point_coords = Dict(nodeTags[i] => nodes[i,:] for i in 1:length_of) # generates a dictionary
# where each tag represents the coordinates of the points that the tags represent

# Use it like this:
# coord_map = get_mesh_node_coord_map() # replace with your node tag

points = Tuple{Float64,Float64,Float64,Float64,Int}[]

function get_all_edges(triangles::Matrix)
    edges = Tuple{Int, Int}[]
    for tri in eachrow(triangles)
        push!(edges, sort((tri[1], tri[2])))
        push!(edges, sort((tri[2], tri[3])))
        push!(edges, sort((tri[3], tri[1])))
    end
    return edges
end

function count_edge_occurrences(edges::Vector{Tuple{Int,Int}})
    edge_count = Dict{Tuple{Int,Int}, Int}()
    for e in edges
        edge_count[e] = get(edge_count, e, 0) + 1
    end
    return edge_count
end

function get_boundary_edges(edge_count::Dict{Tuple{Int,Int}, Int})
    return [e for (e, count) in edge_count if count == 1]
end

function get_boundary_nodes(boundary_edges::Vector{Tuple{Int,Int}})
    return unique(vcat([e[1] for e in boundary_edges], [e[2] for e in boundary_edges])...)
end

function find_boundary_node_tags(triangles::Matrix)
    edges = get_all_edges(triangles)
    edge_count = count_edge_occurrences(edges)
    boundary_edges = get_boundary_edges(edge_count)
    boundary_nodes = get_boundary_nodes(boundary_edges)
    return boundary_nodes, boundary_edges
end

# Use it
# boundary_nodes, boundary_edges = find_boundary_node_tags(triangles)

# Tag, coor , _ = gmsh.model.mesh.getNodes(0, 27, true)
# println("here is the tag ", Tag, " and the coordinates ", coor)


edges = Set{Tuple{Int, Int, Bool}}()
interior_edges = Set{Tuple{Int, Int, Bool}}()
edges_to_triangles = Dict{Tuple{Int, Int}, Vector{Int}}()

function create_edges(triangles_matrix, edges_set, interior_edge_set)
    for i in (1:size(triangles_matrix)[1])
        # println("i ", i)
        triangle = triangles_matrix[i,:]
        # println("tirangles, ", triangle)
        edge1 = (triangle[1], triangle[2], false)
        edge1_flipped = (triangle[2], triangle[1], false)
        edge2 = (triangle[2], triangle[3], false)
        edge2_flipped = (triangle[3], triangle[2], false)
        edge3 = (triangle[3], triangle[1], false)
        edge3_flipped = (triangle[1], triangle[3], false)

        if edge1 in edges_set 
            push!(interior_edge_set, edge1, false)
        elseif edge1_flipped in edges_set
            push!(interior_edge_set, edge1_flipped)
        else
            push!(edges_set, edge1)
            edges_to_triangles[(edge1[1], edge1[2])] = triangle
        end

        if edge2 in edges_set 
            push!(interior_edge_set, edge2)
        elseif edge2_flipped in edges_set
            push!(interior_edge_set, edge2_flipped)
        else
            push!(edges_set, edge2)
            edges_to_triangles[(edge2[1], edge2[2])] = triangle
            edges_to_triangles[(edge2[2], edge2[1])] = triangle
        end

        if edge3 in edges_set 
            push!(interior_edge_set, edge3)
        elseif edge3_flipped in edges_set
            push!(interior_edge_set, edge3_flipped)
        else
            push!(edges_set, edge3)
            edges_to_triangles[(edge3[1], edge3[2])] = triangle
            edges_to_triangles[(edge3[2], edge3[1])] = triangle
        end
    end
    return setdiff(edges_set, interior_edge_set)
end

boundary_edges_unordered = create_edges(triangles, edges, interior_edges)
boundary_edges_unordered_vector = sort(collect(boundary_edges_unordered))#Tuple{Int, Int, Bool}[]

nodes_on_the_boundary = Set()

for elem in boundary_edges_unordered_vector
    println(elem)
    push!(nodes_on_the_boundary, elem[1])
    push!(nodes_on_the_boundary, elem[2])
end


# for loop in ordered_surface_lines
#     for edge in loop
#         push!(nodes_on_the_boundary, edge[2])
#     end
# end


# Look up as needed

nodeTags, coords, actual_tags = gmsh.model.mesh.getNodes(-1)
point_tags  = last.(gmsh.model.getEntities(0))
function node_coord_map()
    # coords is flat: [x1,y1,z1, x2,y2,z2, ...]
    n = length(nodeTags)
    coord_of = Dict{Int, NTuple{3,Float64}}()
    @inbounds for i in 1:n
        coord_of[Int(nodeTags[i])] = (coords[3i-2], coords[3i-1], coords[3i])
    end
    return coord_of
end

# Build once
coord_of = node_coord_map()

boundary_edges_ordered = Vector{Vector{NTuple{4, Int}}}()


# for id in nodes_on_the_boundary
#     coor = coord_map[id]
#     push!(points, (coor[1], coor[2], coor[3], h_0, id))
# end

geometric_ID_to_node_ID = Dict{Int, Tuple}()

# for loop in ordered_surface_lines
#     for i in (1:length(loop))
#         node1_ID, node1_cor , _ = gmsh.model.mesh.getNodes(0, loop[i][2], false)
#         node2_ID, node2_cor , _ = gmsh.model.mesh.getNodes(0, loop[i][3], false)
#         geometric_ID_to_node_ID[loop[i][2]] = (node1_ID, node1_cor[1], node1_cor[2], node1_cor[3])
#         geometric_ID_to_node_ID[loop[i][3]] = (node2_ID, node2_cor[1], node2_cor[2], node2_cor[3])
#     end
# end

gmsh.finalize()

curve_loops = [Int[] for _ in 1:length(edges_on_the_boundary_ordered)]

added_points_to_the_mesh = Set()

function fix_boundary(curve_to_tri, bound_nodes, added_points)
    for (loop_index, loop) in enumerate(edges_on_the_boundary_ordered)
        index = 1
        try_reaching_the_node = false
        node_to_reach = nothing
        while index <= length(loop)
        # for i in (1:length(loop))
            first_node_index = nothing
            second_node_index = nothing

            node1_ID = loop[index][2]
            node1_cor_vector = tag_to_point_coords[node1_ID]
            node1_cor = (node1_cor_vector[1], node1_cor_vector[2], node1_cor_vector[3])

            node2_ID = loop[index][3]
            node2_cor_vector = tag_to_point_coords[node2_ID]
            node2_cor = (node2_cor_vector[1], node2_cor_vector[2], node2_cor_vector[3])

            # println("node1 ", node1_ID)
            edge = (node1_ID,node2_ID)
            edge_other_way = (node2_ID,node1_ID)
            # println("edge ", edge, " other way ", edge_other_way)

            tri = [1,2,3]
            if haskey(curve_to_tri, edge)
                tri = curve_to_tri[edge]
            elseif haskey(curve_to_tri, edge_other_way)
                tri = curve_to_tri[edge_other_way]
            end
            if tri != [31, 13249, 25]
                println("here is the triangle ", tri)
            end
            for i in 1:3
                if tri[i] == node1_ID
                    first_node_index = i 
                end
                if tri[i] == node2_ID
                    second_node_index = i 
                end 
            end
            third_node_index = 6 - first_node_index - second_node_index
            third_node = tri[third_node_index]
           
            if third_node in bound_nodes
                try_reaching_the_node = true
                node_to_reach = third_node
                

            else
                if !try_reaching_the_node
                    if !(node1_ID in added_points)
                        gmsh.model.occ.addPoint(node1_cor[1], node1_cor[2], node1_cor[3], h_0, node1_ID)
                        push!(added_points, node1_ID)
                    end
                    if !(node2_ID in added_points)
                        gmsh.model.occ.addPoint(node2_cor[1], node2_cor[2], node2_cor[3], h_0, node2_ID)
                        push!(added_points, node2_ID)
                    end
                    line_ID = gmsh.model.occ.addLine(node1_ID, node2_ID)
                    push!(curve_loops[loop_index], line_ID)

                end
                # if node1_ID == node_to_reach || node2_ID == node_to_reach 
                #     try_reaching_the_node = false
                # end
            end

            if try_reaching_the_node
                loop1_length = 0
                loop2_length = 0
                
                for j in (index:length(loop))
                    if loop[j][3] == node_to_reach
                        loop1_length += 1
                        break
                    else
                        println("node: ", loop[j][3])
                        loop1_length += 1
                    end
                end

                println("loop length1 ", loop1_length)

                for j in ( index+loop1_length : length(loop))
                    if loop[j][3] == node1_ID
                        
                        break
                    else
                        println("node: ", loop[j][3])
                        loop2_length += 1
                    end
                end

                println("loop length1 ", loop1_length, " loop length 2 ", loop2_length)

                if loop1_length < loop2_length
                    index += loop1_length

                    if !(node1_ID in added_points)
                        gmsh.model.occ.addPoint(node1_cor[1], node1_cor[2], node1_cor[3], h_0, node1_ID)
                        push!(added_points, node1_ID)
                    end
                    
                    # if !(node2_ID in added_points)
                    #     gmsh.model.occ.addPoint(node2_cor[1], node2_cor[2], node2_cor[3], h_0, node2_ID)
                    #     push!(added_points, node2_ID)
                    # end

                    third_node_coor_vector = tag_to_point_coords[third_node]
                    third_node_coor = (third_node_coor_vector[1], third_node_coor_vector[2], third_node_coor_vector[3])
                    gmsh.model.occ.addPoint(third_node_coor[1], third_node_coor[2], third_node_coor[3], h_0, third_node)
                    push!(added_points, third_node)

                    line_ID = gmsh.model.occ.addLine(node1_ID, third_node)
                    push!(curve_loops[loop_index], line_ID)

                else
                    if !(node1_ID in added_points)
                        gmsh.model.occ.addPoint(node1_cor[1], node1_cor[2], node1_cor[3], h_0, node1_ID)
                        push!(added_points, node1_ID)
                    end
                    
                    if !(node2_ID in added_points)
                        gmsh.model.occ.addPoint(node2_cor[1], node2_cor[2], node2_cor[3], h_0, node2_ID)
                        push!(added_points, node2_ID)
                    end

                    line_ID = gmsh.model.occ.addLine(node1_ID, node2_ID)
                    push!(curve_loops[loop_index], line_ID)
                    index += 1
                    
                end
                println("Index: ", index)
                try_reaching_the_node = false
            else
                index += 1
            end
        end
    end
end


gmsh.initialize()
gmsh.model.add("Boundary_Fixed_Part2.msh")

# create_loop(boundary_edges_unordered_vector, boundary_edges_ordered)

fix_boundary(edges_to_triangles, nodes_on_the_boundary, added_points_to_the_mesh)

gmsh.model.occ.synchronize()

list_of_loops = Int[]
list_of_surfaces = Int[]

function create_curve_loops(edge_number_of_loops)
    for loop in edge_number_of_loops
        loop_no = gmsh.model.occ.addCurveLoop(loop)
        push!(list_of_loops, loop_no)
    end
end

function create_surfaces_based_on_loops(loop_list, surfaces)
    for elem in loop_list
        surface = gmsh.model.occ.addPlaneSurface([elem])
        push!(surfaces, surface)
    end
end

create_curve_loops(curve_loops)
create_surfaces_based_on_loops(list_of_loops, list_of_surfaces)

boundary_edge_to_triangle = Dict{Tuple{Int, Int}, Vector{Int}}()

ordered_boundary_edges = Vector{Vector{NTuple{4, Int}}}()

function edge_to_triangle()
    for i in (1:size(triangles)[1])
        triangle = triangles[i,:]
        edge1 = (triangle[1], triangle[2])
        edge2 = (triangle[2], triangle[3])
        edge3 = (triangle[3], triangle[1])

        if (check_if_tuple_in_vector(boundary_edges, edge1))
            boundary_edge_to_triangle[edge1] = triangle
        end
        if (check_if_tuple_in_vector(boundary_edges, edge2))
            boundary_edge_to_triangle[edge2] = triangle
        end
        if (check_if_tuple_in_vector(boundary_edges, edge3))
            boundary_edge_to_triangle[edge3] = triangle
        end
    end
end

# edge_to_triangle()

function is_same(tuple_main, tuple)
    return tuple_main == tuple 
end



function get_node_coordinates(node_id)::Tuple{Float64, Float64, Float64}
    idx = findfirst(==(node_id), nodeTags)

    if isnothing(idx)
        error("Node ID $node_id not found in mesh.")
    end

    x = coords[3*(idx - 1) + 1]
    y = coords[3*(idx - 1) + 2]
    z = coords[3*(idx - 1) + 3]
    return (x, y, z)
end

added_nodes = Int[]

function point_exists(tag)
    return tag in added_nodes
end


@save "curve_loops" curve_loops

gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(3)
gmsh.write("Boundary_Fixed_Part2.msh")
gmsh.finalize()