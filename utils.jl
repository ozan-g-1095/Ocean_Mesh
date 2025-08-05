using ProgressMeter

function set_mesh_of_cells()
@showprogress 1 ("Processing each cell for mesh size...") for i in (argmin(abs.(lon .- lon_min)): argmin(abs.(lon .- lon_max))) # longitude
        for j in (argmin(abs.(lat .- lat_min)):argmin(abs.(lat .- lat_max))) # latitude

            x_min = (deg2rad(lon[i]) * R - x_min_normalize ) / (x_max_normalize - x_min_normalize)
            x_max = (deg2rad(lon[i+1]) * R -x_min_normalize) / (x_max_normalize - x_min_normalize)
            y_min = (deg2rad(lat[j]) * R - y_min_normalize) / (x_max_normalize - x_min_normalize)
            y_max = (deg2rad(lat[j+1]) * R -y_min_normalize) / (x_max_normalize - x_min_normalize)

            corner1 = (z[i, j] * α) /z_max # Q11

            corner2 = (z[i+1, j]*α)/z_max # Q21

            corner3 = (z[i, j+1]*α)/z_max # Q12

            corner4 = (z[i+1, j+1]*α)/z_max # Q22

            xs =  [x_min, x_max]
            ys =  [y_min, y_max]
            # xs = [deg2rad(lon[i]) * R, deg2rad(lon[i]+displacement_lon) * R] #the x axis is longitude
            # ys = [deg2rad(lat[j]) * R, deg2rad(lat[j]+displacement_lat) * R] #the y axis is latitude

            # println(" here is xs: ", xs ," and ys: ", ys)
            
            A = [corner1; corner3; corner2; corner4] # This is the matrix of corners needed to calculate the coefficients that represent the bilinear interpolation
            # for the specific grid space defined by the corners.

            coef = bilinear_coefficients(xs[1], xs[2], ys[1], ys[2], A)
            # println("Coefficients: ", coef, " for cell (", i, ",", j, ")")
            cell_coefficients[i, j] = (coef[1], coef[2], coef[3], coef[4])

            # println(typeof(coef))
            # println(coef)
            # cell_coefficients[i, j] = coef

            x_middle = (x_min + x_max)/2
            y_middle = (y_min + y_max)/2

            (Hx, Hy) = polynomial_fit_derivative(coef, x_middle, y_middle)

            # println("Hx ", Hx, " Hy ", Hy)
            desired_mesh_size_in_cell = (h)/(sqrt(1  + Hx^2 + Hy^2)) # *0.1
            # println("Desired mesh size in cell: ", desired_mesh_size_in_cell)

            println(desired_mesh_size_in_cell)
            # println("--------------")

            field_id = gmsh.model.mesh.field.add("Box")
            println(field_id)
            gmsh.model.mesh.field.setNumber(field_id, "VIn", desired_mesh_size_in_cell)
            gmsh.model.mesh.field.setNumber(field_id, "VOut", h)
            gmsh.model.mesh.field.setNumber(field_id, "XMin", x_min)
            gmsh.model.mesh.field.setNumber(field_id, "XMax", x_max)
            gmsh.model.mesh.field.setNumber(field_id, "YMin", y_min)
            gmsh.model.mesh.field.setNumber(field_id, "YMax", y_max)
            # gmsh.model.mesh.field.setNumber(field_id, "ZMin", z_min)  # set appropriately for your domain
            # gmsh.model.mesh.field.setNumber(field_id, "ZMax", z_max)

            # println("Hello")
            push!(field_ids, field_id)
        end
    end
    min_field = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(min_field, "FieldsList", field_ids)
    gmsh.model.mesh.field.setAsBackgroundMesh(min_field)
end

function cont()
    plt = Plots.contour(real_bathymetry_x, real_bathymetry_y, real_bathymetry_z', levels=21, fill=true, c=:balance)

    # plt = Plots.scatter(real_bathymetry_x, real_bathymetry_y, real_bathymetry_z)
    Plots.contour!(real_bathymetry_x, real_bathymetry_y, real_bathymetry_z'; levels=[0], fill=false,
        linewidth=2, linecolor=:red, colorbar=false, color =:balance) 
    savefig(plt, "bathymetry_contour.png")
    # display(plt)
    # f(x, y) = (3x + y^2) * abs(sin(x) + cos(y))
    # x = range(0, 5, length=100)
    # y = range(0, 3, length=50)
    # z = @. f(x', y)
    # plt = Plots.contour(x, y, z)
    # display(plt)
end

function create_surfaces(loop_list)
    for elem in curves_to_create_closed_loops
        surface = gmsh.model.occ.addPlaneSurface([elem])
        push!(surfaces, surface)
    end
end

function create_loops(surface_lines_ordered, loop_list)
    for elem in surface_lines_ordered
        curve_individual = gmsh.model.occ.addCurveLoop([line[1] for line in elem])
        push!(loop_list, curve_individual)
    end
end

function delete_small_loops(surface_lines_ordered::Vector{Vector{NTuple{4, Int64}}})
    for idx in length(surface_lines_ordered):-1:1
        elem = surface_lines_ordered[idx]
        if length(elem) < 40
            for i in (1:length(elem))
                println(elem[i])
                gmsh.model.occ.remove([(1, elem[i][1])], true)
            end
            deleteat!(surface_lines_ordered, idx)
        end
    end
end

function delete_small_lines_in_boundary(surface_lines_ordered)
    for loops in (length(surface_lines_ordered):-1:1)
        for line in (length(surface_lines_ordered[loops]):-1:1)

            if line == 1
                if surface_lines_ordered[loops][line][4] == 1 || surface_lines_ordered[loops][end][4] == 1
                    gmsh.model.occ.remove(([1, surface_lines_ordered[loops][line][1]]))
                    gmsh.model.occ.remove(([1, surface_lines_ordered[loops][end][1]]))

                    gmsh.model.occ.remove(([0, surface_lines_ordered[loops][line][2]]), true)
                    gmsh.model.occ.remove(([0, surface_lines_ordered[loops][end][3]]), true)

                    println("points: ", surface_lines_ordered[loops][end][2], " ", surface_lines_ordered[loops][line][3])
                    line_no = gmsh.model.occ.addLine(surface_lines_ordered[loops][end][2], surface_lines_ordered[loops][line][3])

                    surface_lines_ordered[loops][line-1] = (line_no, surface_lines_ordered[loops][end][2], surface_lines_ordered[loops][line][3], 0)
                    deleteat!(surface_lines_ordered[loops], line)
            
                end

            else

                if surface_lines_ordered[loops][line][4] == 1 || surface_lines_ordered[loops][line-1][4] == 1
                    gmsh.model.occ.remove(([1, surface_lines_ordered[loops][line][1]]))
                    gmsh.model.occ.remove(([1, surface_lines_ordered[loops][line-1][1]]))

                    gmsh.model.occ.remove(([0, surface_lines_ordered[loops][line][2]]), true)
                    gmsh.model.occ.remove(([0, surface_lines_ordered[loops][line-1][3]]), true)

                    println("points: ", surface_lines_ordered[loops][line-1][2], " ", surface_lines_ordered[loops][line][3])
                    line_no = gmsh.model.occ.addLine(surface_lines_ordered[loops][line-1][2], surface_lines_ordered[loops][line][3])

                    surface_lines_ordered[loops][line-1] = (line_no, surface_lines_ordered[loops][line-1][2], surface_lines_ordered[loops][line][3], 0)
                    deleteat!(surface_lines_ordered[loops], line)
            
                end
            end
        end
    end
end

function check_if_tuple_in_vector(vector, tuple::Tuple)
    return ((tuple[1], tuple[2]) in vector) || ((tuple[2], tuple[1]) in vector)
end

function create_loop(lines_on_surface, surface_lines_ordered)
    global loop_no = 1
    while !isempty(lines_on_surface)
        initial = lines_on_surface[1][1]
        current = lines_on_surface[1][2]
        current_index = 1

        line_no = gmsh.model.occ.addLine(initial, current)
        push!(surface_lines_ordered, [(Int64(line_no), initial, current, lines_on_surface[1][3])])
        println((Int64(line_no), initial, current))

        while (initial != current)
            global loop_formed = false
            
            deleteat!(lines_on_surface, current_index) # it deletes the chosen point so that the algorithm can march forward and find the 
            # next grid that shares the chosen current's tag. 
            has_found_connectivity = false
            for (idx,elem) in enumerate(lines_on_surface)
               
                if elem[1] == current
                    
                    line_no = gmsh.model.occ.addLine(elem[1], elem[2])
                    
                    push!(surface_lines_ordered[loop_no], (Int64(line_no), elem[1], elem[2], elem[3]))
                
                    current = elem[2]
                    current_index = idx
                    has_found_connectivity = true
                    break
                    
                elseif elem[2] == current
                    line_no = gmsh.model.occ.addLine(elem[2], elem[1])
                    
                    push!(surface_lines_ordered[loop_no], (Int64(line_no), elem[2], elem[1], elem[3]))
                
                    current = elem[1]
                    current_index = idx
                    has_found_connectivity = true  
                    break
                end
            end

            if !(has_found_connectivity)
                for i in (length(surface_lines_ordered[loop_no]):-1:1)
                    gmsh.model.occ.remove(([1,surface_lines_ordered[loop_no][i][1]])) # recursive = true
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