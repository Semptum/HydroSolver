function mesh_to_polygons(mesh::PlaneMesh{FloatT}) where FloatT
    N_cells = length(mesh.cCenters)
    N_nodes = length(mesh.nodes)
    faces = Matrix{Int64}(undef,(N_cells*2,3))
    vertices = Matrix{Float64}(undef,(N_cells*4,2))
    for i in 1:N_cells
        edges = mesh.fNodes[mesh.cells[i]]
        connec = edges_to_connectivity(edges)
        list = connectivity_to_list(connec)
        vertices[i,:] .= mesh.nodes[list[1]]
        vertices[i+N_cells,:] .= mesh.nodes[list[2]]
        vertices[i+2*N_cells,:] .= mesh.nodes[list[3]]
        vertices[i+3*N_cells,:] .= mesh.nodes[list[4]]
        faces[i,:] .= [i, i+N_cells, i+2*N_cells]
        faces[i+N_cells,:] .=  [i+2*N_cells,i+3*N_cells, i]
    end
    return faces, vertices
end

function col(val)
    return vcat(val,val,val,val)
end

function edges_to_connectivity(edges)
    C = Dict{Int,Vector{Int}}(j=>[]  for i in edges for j in i)
    for i in edges
        push!(C[i[1]],i[2])
        push!(C[i[2]],i[1])
    end
    return C
end

function connectivity_to_list(connectivity)
    list = Vector{Int}(undef,length(connectivity))
    list[1] = first(keys(connectivity))
    for i in 2:length(list)
        if connectivity[list[i-1]][1] in list
            list[i] = connectivity[list[i-1]][2]
        else
            list[i] = connectivity[list[i-1]][1]
        end
    end
    return list
end


function continuous_interpolation(mesh::PlaneMesh,value, point)
    i_closest = findmin(norm.([point] .- mesh.cCenters))[2]
    return value[i_closest]
end

function find_in_range(x_range,x)
    for i in 1:length(x_range)
        if x_range[i]>x
            return i-1
        end
    end
    return length(x_range)
end