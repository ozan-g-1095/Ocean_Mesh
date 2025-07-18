using Interpolations
using Plots
gr()

using NCDatasets

using ProgressMeter
using Interpolations
using JLD2

using Contour
using LinearAlgebra

using Gmsh: gmsh
R = 6378000
gmsh.initialize()
gmsh.model.add("Biliniear_Correlation")

ds = NCDataset("topo_25.1_coarsened1024.nc", "r")

lon = ds["lon"][:]
lat = ds["lat"][:]
z   = ds["z"][:,:]
close(ds)

lon_min =26.8#-179
lon_max =30 #179

lat_min= 40.13#-79
lat_max=  41#79

displacement_lon = abs(lon[1] - lon[2])
displacement_lat = abs(lat[1] - lat[2])

min_dist = sqrt((displacement_lat/4)^2 + (displacement_lon/4)^2 )


x_min = deg2rad(lon[argmin(abs.(lon .- lon_min))]) * R
y_min = deg2rad(lat[argmin(abs.(lat .- lat_min))]) * R
x_max = deg2rad(lon[argmin(abs.(lon .- lon_max))]) * R
L = x_max - x_min

z_max = 0
for i in (argmin(abs.(lon .- lon_min)): argmin(abs.(lon .- lon_max))) 
    for j in (argmin(abs.(lat .- lat_min)):argmin(abs.(lat .- lat_max)))
        if z[i,j]> z_max
            global z_max = z[i,j]
        end
    end 
end

"""
This function calculates the bilinear coefficients for a given quadrilateral defined by its corner points (x1, x2, y1, y2) 
    and the function values at those corners (fQ).

y2        Q12--------Q22
           |          |
           |          |
y1        Q11--------Q21

           x1         x2

It returns the coefficients as a tuple (a1, a2, a3, a4) where:
    f = a1 + a2*x + a3*y + a4*x*y

"""

function bilinear_coefficients(x1, x2, y1, y2, fQ)    
    denom = (x2 - x1) * (y2 - y1)
    
    M = [x2*y2  -x2*y1  -x1*y2  x1*y1;
         -y2     y1      y2    -y1;
         -x2     x2      x1    -x1;
         1      -1      -1     1]
    
    a = (1/denom) *  M * fQ 
    return a
end

"""
This function evaluates a polynomial at given x and y coordinates using the coefficients provided.

"""
function polynomial_fit(coefficients, x, y)
    return coefficients[1] + coefficients[2]*x + coefficients[3]*y + coefficients[4] *x *y

end 


function polynomial_fit_derivative(coefficients, x, y)
    return  (coefficients[2]   , coefficients[3] )
    # returns (Hx, Hy) + coefficients[4] *y + coefficients[4] *x

end 

"""
This function finds the lattitude of the 0 contour defined by a vertical strip given by the longitude value

    |      |       |
    |      |       |
    |      |       |
    |      |       |
    |      |       |
  lon-1   lon    lon+1 

  It retuns 0 if the value does not fit between the min and max allowed values of latitude for a grid
"""
function find_0_contour_latitude_on_grid(coefficients, lon, ymin, ymax)
    x = lon
    y = -(coefficients[1] + coefficients[2]*x) / (coefficients[3] + (coefficients[4] *x))
    if ( ymin<= y <= ymax)
        return y
    else
        return 0
    end
end


"""
This function finds the longitude of the 0 contour defined by a horizontal strip given by the latitude value


    It retuns 0 if the value does not fit between the min and max allowed values of longitude for a grid

"""

function find_0_contour_longitude_on_grid(coefficients, lat, xmin, xmax)
    y = lat
    x = -(coefficients[1] + coefficients[3]*y) / (coefficients[2] + (coefficients[4] *y))
    if ( xmin<= x <= xmax)
        return x
    else
        return 0
    end
end

"""
This function returns true if a given values in the array all have the same signs (treats 0 as positive) false otherwise

Used to determine whether a 0 contour will pass in a given grid. Helps reduce computational cost by determining which cells
should be processed for coastline finding. 
"""
function all_same_sign(arr)
    arr = filter(!=(0), arr)
    s = sign.(arr)
    return all(s .== s[1])
end

x_points = Float32[]
y_points = Float32[]
z_points = Float32[]


surface_lines= Tuple{Int, Int, Bool}[]

"""
This function looks at the neighbors of a passed grid by values (i,j). It checks all 6 neighbors of the grid space to see which one it has a 
    connection with its lines. If the point on the edge of the grid is shared by one of the neighbors and it has already been placed by the 
    neighbor, it returns the existing tag of that point. If none of the neighbors have placed that point it creates the point and returns
    its tag. 

"""


function map_in_3D(lat, lon, i, j )
        for row in i-1:i+1
            for col in j-1:j+1
                if isapprox(lat, grid_space[row, col].lat1, atol=1e-6) && isapprox(lon, grid_space[row, col].lon1, atol=1e-6)
                    return (1,Int(grid_space[row, col].point1))
                elseif isapprox(lat, grid_space[row, col].lat2, atol=1e-6) && isapprox(lon, grid_space[row, col].lon2, atol=1e-6)
                    return (1,Int(grid_space[row, col].point2))
                elseif isapprox(lat, grid_space[row, col].lat3, atol=1e-6) && isapprox(lon, grid_space[row, col].lon3, atol=1e-6)
                    return (1,Int(grid_space[row, col].point3))
                elseif isapprox(lat, grid_space[row, col].lat4, atol=1e-6) && isapprox(lon, grid_space[row, col].lon4, atol=1e-6)
                    return (1,Int(grid_space[row, col].point4))
                end
            end
        end
        # If no point found, create a new point
        lat_point = deg2rad(lat)
        long_point = deg2rad(lon)

        x_cor = R * cos(lat_point) * cos(long_point)
        y_cor = R * cos(lat_point) * sin(long_point)
        z_cor = R * sin(lat_point)
         
        # no = gmsh.model.occ.addPoint(x_cor, y_cor, z_cor)
        no = gmsh.model.occ.addPoint(deg2rad(lon) * R, deg2rad(lat) * R, 0)
        push!(x_points, deg2rad(lon) * R)
        push!(y_points, deg2rad(lat) * R )
        push!(z_points, 0)

        return (0,Int(no))

end

boundary = Dict{Int, Tuple{Float64, Float64}}()
edgeSets = Dict{Tuple{Int, Int}, Int}() 
boundaryEdgeSet = Int[]

m = length(lon) 
n =  length(lat) 


"""
This data structure is used for each cell in the grid matrix that represent the lat and long combinations.

has_line : 0 indicates the cell doesn't have a line passing through it (no 0 contour), 1 indicates a 0 contour passes
latn: the latitude of point n.
lonn: the longitude of point n.
pointn: the tag of the given node n.

Since at most 4 points can be represented in a grid there are 4 points. The points are filled in an increasing order. First point1 then 2 and so forth.

Most of the time point3 and point4 are empty when there is a contour of 0 passing.

"""
struct GridCell
    has_line::Int

    lat1::Float64
    lon1::Float64
    point1::Int

    lat2::Float64
    lon2::Float64
    point2::Int

    lat3::Float64
    lon3::Float64
    point3::Int

    lat4::Float64
    lon4::Float64
    point4::Int
end

# This creates a default_cell to fill the matrix called the grid space initially which will be overwritten as the
# cells are starting to get filled. 
default_cell = GridCell(0, 0.0, 0.0, 0, 0.0, 0.0, 0, 0.0, 0.0, 0, 0.0, 0.0, 0)
grid_space = fill(default_cell, m, n)


"""
This function finds the bilinear interpolation for the given grid space defined by the lon and lat values. Based on this it 
    determines the 0 contour and creates the grid space with the points that are needed to be created. Based on the points on the 
    contour it creates adds points to the data structure surface_lines. This is used later to create the correctly directed lines 
    for createing the surface mesh. 

"""
real_bathymetry_x = Float64[]#range(argmin(abs.(lon .- lon_min)), argmin(abs.(lon .- lon_max)) ) |> collect # to store the real bathymetry values for plotting later
real_bathymetry_y = Float64[]#range(argmin(abs.(lat .- lat_min)), argmin(abs.(lat .- lat_max)) ) |> collect


# real_bathymetry_z = Matrix{Float64}(undef, argmin(abs.(lon .- lon_max)) - argmin(abs.(lon .- lon_min))+1, argmin(abs.(lat .- lat_max)) - argmin(abs.(lat .- lat_min))+1) # to store the real bathymetry values for plotting later
 real_bathymetry_z = Float64[]

function find_biliniear_interpolation()
    @showprogress 1 ("Processing longitude loop...") for i in (argmin(abs.(lon .- lon_min)): argmin(abs.(lon .- lon_max))) # longitude
        for j in (argmin(abs.(lat .- lat_min)):argmin(abs.(lat .- lat_max))) # latitude
            corner1 = z[i, j] # Q11

            corner2 = z[i+1, j] # Q21

            corner3 = z[i, j+1] # Q12

            corner4 = z[i+1, j+1] # Q22

            push!(real_bathymetry_x, lon[i])
            push!(real_bathymetry_y, lat[j])
            push!(real_bathymetry_z, corner1)

            # real_bathymetry_z[i-argmin(abs.(lon .- lon_min))+1, j-argmin(abs.(lat .- lat_min))+1] = corner1

            xs = [lon[i], lon[i]+displacement_lon] #the x axis is longitude
            ys = [lat[j], lat[j]+displacement_lat] #the y axis is latitude
            
            A = [corner1; corner3; corner2; corner4] # This is the matrix of corners needed to calculate the coefficients that represent the bilinear interpolation
            # for the specific grid space defined by the corners.

            coef = bilinear_coefficients(xs[1], xs[2], ys[1], ys[2], A)
            println("Coefficients: ", coef)
 
            # these are the possible points on the grid space. 
            point1 = (0,0,0,0) # lat, long, tag id, prev_existing
            point2 = (0,0,0,0)
            point3 = (0,0,0,0)
            point4 = (0,0,0,0)

            # if all of the points have the same sign, then there is no 0 contour passing through this grid space
            if !(all_same_sign(A))
                
                    # these are the possible points on the grid space that can be used to create the 0 contour
                    # possible_yn represent the latitude of the 0 contour that passes through the vertical strip defined by xs[n]
                    # that is also the case for possible_xn but defined by the horizontal strips ys[n]
                    possible_y1 = find_0_contour_latitude_on_grid(coef,  xs[1], ys[1], ys[2])
                    possible_y2 = find_0_contour_latitude_on_grid(coef,  xs[2], ys[1], ys[2])
                    possible_x1 = find_0_contour_longitude_on_grid(coef, ys[1], xs[1], xs[2])
                    possible_x2 = find_0_contour_longitude_on_grid(coef, ys[2], xs[1], xs[2])

                    # checks that if the possible points are not 0, then it maps them to the grid space and creates the point with the tag id,
                    # its lat and long values and whether the point already exists in the grid space.
                    if (possible_y1 != 0)
                        (prev_existing, number) = map_in_3D(possible_y1, xs[1], i ,j)
                        point1 = (possible_y1, xs[1], number, prev_existing)
                    end
                    if (possible_y2 != 0)
                        (prev_existing, number) = map_in_3D(possible_y2, xs[2], i ,j)
                        if point1 != (0,0,0,0)
                            point2 = (possible_y2, xs[2], number, prev_existing)
                        else
                            point1 = (possible_y2, xs[2], number, prev_existing)
                        end
                    end
                    if (possible_x1 != 0)
                        (prev_existing, number) = map_in_3D(ys[1], possible_x1, i ,j)
                        if point2 != (0,0,0,0)
                            point3 = (ys[1], possible_x1, number, prev_existing) 
                        elseif point1 != (0,0,0,0)
                            point2 = (ys[1], possible_x1, number, prev_existing)
                        else
                            point1 = (ys[1], possible_x1, number, prev_existing)
                        end
                    end
                    if (possible_x2 != 0)
                        (prev_existing, number) = map_in_3D(ys[2], possible_x2,  i ,j)
                        if point3 != (0,0,0,0)
                            point4 = (ys[2], possible_x2, number, prev_existing)
                        else
                            point2 = (ys[2], possible_x2, number, prev_existing)
                        end
                    end 


                # if points 3 and 4 are not empty it means that there is a 0 contour passing through every edge of the grid and 
                # the points should be connected correctly because there are 2 poissble ways to connect them. 

                #         3
                #      [     ]
                #   1  |     |  2
                #      [     ]
                #         4

                # so the lines that pass through can either be from 1-3 and 4-2 or 1-4 and 3-2
                if point3 != (0,0,0,0) && point4 != (0,0,0,0)
                    if sign(corner3) == sign(corner2)
                        
                        # sqrt(abs(point1[1] - point3[1])^2 + abs(point1[2] -point3[2])^2 determines whether the line is too short 
                        # and needs to be cut to prevent the mesh from having very small connections and thus very small triangles.
                        # Helps keep the mesh size fairly uniform throughout.

                        push!(surface_lines, (point1[3], point3[3], sqrt(abs(point1[1] - point3[1])^2 + abs(point1[2] -point3[2])^2 ) < min_dist))

                        
                        push!(surface_lines, (point2[3], point4[3], sqrt(abs(point2[1] - point4[1])^2 + abs(point2[2] -point4[2])^2 ) < min_dist))
                    else
                        
                        push!(surface_lines, (point1[3], point4[3], sqrt(abs(point1[1] - point4[1])^2 + abs(point1[2] -point4[2])^2 ) < min_dist))

                       
                        push!(surface_lines, (point2[3], point3[3], sqrt(abs(point2[1] - point3[1])^2 + abs(point2[2] -point3[2])^2 ) < min_dist))
                    end
                else
                
                    push!(surface_lines, (point1[3], point2[3], sqrt(abs(point1[1] - point2[1])^2 + abs(point1[2] -point2[2])^2 ) < min_dist))
                   
                end 
                grid_obj = GridCell(1, point1[1], point1[2], point1[3], point2[1], point2[2], point2[3], 
                                    point3[1], point3[2], point3[3], point4[1], point4[2], point4[3])
                grid_space[i, j] = grid_obj

            end
        end
    end
    gmsh.model.occ.synchronize()
end


find_biliniear_interpolation()

# This adds the surface line tags in order to create a directed closed loop
ordered_surface_lines = NTuple{3, Int}[]



"""
This funciton creates a closed loop based on the points and the line connectivities added to surface_lines

"""
function create_loop()
    initial = surface_lines[1][1]
    current = surface_lines[1][2]
    current_index = 1

    line_no = gmsh.model.occ.addLine(initial, current)
    push!(ordered_surface_lines, (line_no, initial, current))

    no_times = length(surface_lines)
    println("No of times: ", no_times)
    println("surface lines: ", surface_lines)
    
    short_edge_starting_point = 0
    prev_short = false
    while !isempty(surface_lines)
        
        deleteat!(surface_lines, current_index) # it deletes the chosen point so that the algorithm can march forward and find the 
        # next grid that shares the chosen current's tag. 
        
        for (idx,elem) in enumerate(surface_lines)
           
            if elem[1] == current
                if prev_short == true
                    line_no = gmsh.model.occ.addLine(short_edge_starting_point, elem[2])
                    push!(ordered_surface_lines, (line_no, short_edge_starting_point, elem[2]))
                    prev_short = false
                    short_edge_starting_point = 0
                
                else
                    if elem[3] == true
                        prev_short = true
                        short_edge_starting_point = elem[1]
                    else
                        line_no = gmsh.model.occ.addLine(elem[1], elem[2])
                        push!(ordered_surface_lines, (line_no, elem[1], elem[2]))
                    end
                end
                current = elem[2]
                current_index = idx

                # if the 
                if elem[3] == true
                    gmsh.model.occ.remove([0, elem[2]])
                end      
                break

            elseif elem[2] == current
                if prev_short == true
                    line_no = gmsh.model.occ.addLine(short_edge_starting_point, elem[1])
                    push!(ordered_surface_lines, (line_no, short_edge_starting_point, elem[1]))
                    prev_short = false
                    short_edge_starting_point = 0
            
                else
                    if elem[3] == true
                        prev_short = true
                        short_edge_starting_point = elem[2]
                    else
                        line_no = gmsh.model.occ.addLine(elem[2], elem[1])
                        push!(ordered_surface_lines, (line_no, elem[2], elem[1]))
                    end
                end

                current = elem[1]
                current_index = idx
                
                if elem[3] == true
                    gmsh.model.occ.remove([0, elem[1]])
                end

                break
    
            end
        end
    end
end

create_loop()

curve1 = gmsh.model.occ.addCurveLoop([line[1] for line in ordered_surface_lines])

surface1 = gmsh.model.occ.addPlaneSurface([curve1])

gmsh.model.occ.synchronize()
h = 30000 # mesh size

field_ids = Int[]
cell_coefficients = Matrix{Tuple{Float64, Float64, Float64, Float64}}(undef, m, n) #fill((0.0, 0.0, 0.0, 0.0), m, n)

function set_mesh_of_cells()
@showprogress 1 ("Processing each cell for mesh size...") for i in (argmin(abs.(lon .- lon_min)): argmin(abs.(lon .- lon_max))) # longitude
        for j in (argmin(abs.(lat .- lat_min)):argmin(abs.(lat .- lat_max))) # latitude

            x_min = deg2rad(lon[i]) * R
            x_max = deg2rad(lon[i+1]) * R
            y_min = deg2rad(lat[j]) * R
            y_max = deg2rad(lat[j+1]) * R 

            corner1 = z[i, j] # Q11

            corner2 = z[i+1, j] # Q21

            corner3 = z[i, j+1] # Q12

            corner4 = z[i+1, j+1] # Q22

            xs = [deg2rad(lon[i]) * R, deg2rad(lon[i]+displacement_lon) * R] #the x axis is longitude
            ys = [deg2rad(lat[j]) * R, deg2rad(lat[j]+displacement_lat) * R] #the y axis is latitude
            
            A = [corner1; corner3; corner2; corner4] # This is the matrix of corners needed to calculate the coefficients that represent the bilinear interpolation
            # for the specific grid space defined by the corners.

            coef = bilinear_coefficients(xs[1], xs[2], ys[1], ys[2], A)
            println("Coefficients: ", coef, " for cell (", i, ",", j, ")")
            cell_coefficients[i, j] = (coef[1], coef[2], coef[3], coef[4])

            # println(typeof(coef))
            # println(coef)
            # cell_coefficients[i, j] = coef

            x_middle = deg2rad(lon[i] + displacement_lon/2) * R
            y_middle = deg2rad(lat[j] + displacement_lat/2) * R

            (Hx, Hy) = polynomial_fit_derivative(coef, x_middle, y_middle)

            # println("Hx ", Hx, " Hy ", Hy)
            desired_mesh_size_in_cell = h/(sqrt(1  + Hx^2 + Hy^2))
            # println("Desired mesh size in cell: ", desired_mesh_size_in_cell)

            println(desired_mesh_size_in_cell)
            println("--------------")

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
end

# plt = Plots.contour(real_bathymetry_x, real_bathymetry_y, real_bathymetry_z', levels=10, fill=true, c=:viridis)
# display(plt)
# savefig(plt, "bathymetry_contour.png")

function cont()
    # plt = Plots.contour(real_bathymetry_x, real_bathymetry_y, real_bathymetry_z', levels=10, fill=true, c=:balance)

    plt = Plots.scatter(real_bathymetry_x, real_bathymetry_y, real_bathymetry_z)
    # Plots.contour!(real_bathymetry_x, real_bathymetry_y, real_bathymetry_z'; levels=[0], fill=false,
        #  linewidth=2, linecolor=:red, colorbar=false, color =:balance) 
    display(plt)
    # f(x, y) = (3x + y^2) * abs(sin(x) + cos(y))
    # x = range(0, 5, length=100)
    # y = range(0, 3, length=50)
    # z = @. f(x', y)
    # plt = Plots.contour(x, y, z)
    # display(plt)
end

# cont()
#@save "coefficients_of_cells" cell_coefficients

set_mesh_of_cells()


@save "boundary_line_orders" ordered_surface_lines
@save "cell_coefficients" cell_coefficients
@save "x_min" x_min
@save "y_min" y_min
@save "z_max" z_max
@save "L" L

min_field = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(min_field, "FieldsList", field_ids)

gmsh.model.mesh.field.setAsBackgroundMesh(min_field)

# find_next_point(0,0)

# gmsh.model.mesh.field.add("MathEval", 3)
# gmsh.model.mesh.field.setString(3, "F", "($h*3)") # for circular basin / ( ( 1 + ( 2*x )^2 + ( (2*y)^2 )  ) ^ (1/2) )

# gmsh.model.mesh.field.setAsBackgroundMesh(3)


# gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
# gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
# gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

# gmsh.option.setNumber("Mesh.Algorithm", 5)


gmsh.model.occ.removeAllDuplicates()
# gmsh.option.setNumber("Mesh.SaveAll", 1)

gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(3)
gmsh.write("Biliniear_Correlation.msh")
gmsh.finalize()