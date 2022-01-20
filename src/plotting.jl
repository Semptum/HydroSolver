function mesh_to_polygons(mesh::CFD_Mesh{FloatT}) where FloatT
    N_cells = length(mesh.cCenters)
    list_polygons = Vector{Vector{SVector{2,FloatT}}}(undef,N_cells)
    for i in 1:N_cells
        edges = mesh.fNodes[mesh.cells[i]]
        connec = edges_to_connectivity(edges)
        list = connectivity_to_list(connec)
        list_polygons[i] = SVector{2,FloatT}.(mesh.nodes[list])
    end
    return list_polygons
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

function plot(mesh::CFD_Mesh)
    pol = mesh_to_polygons(mesh)
    f = Figure()
    ax = f[1, 1] = Axis(f)
    for i in pol
        poly!(ax,i)
    end
    return f
end