function rectangle_cartesian_mesher(x_ranges,y_ranges)
    PointT = SVector{2,Float64}
    FloatT = Float64
    Dim = 2
    R=[
         0 1
        -1 0
    ]
    existing_points = Dict{Tuple{Int,Int},Int}()
    face_index_by_points = Dict{Tuple{Int,Int},Int}()
    cell_index_by_points = Dict{Tuple{Int,Int},Int}()
    N_x, N_y = length(x_ranges), length(y_ranges)
    N_faces = N_x*(N_y-1)+N_y*(N_x-1)
    N_cells = (N_x-1)*(N_y-1)
    N_nodes = N_x*N_y

    cells = Vector{Vector{Int}}(undef,N_cells)
    cNeighbours = Vector{Vector{Int}}(undef, N_cells)
    cVols = Vector{FloatT}(undef, N_cells)
    cCenters = Vector{PointT}(undef, N_cells)
    cLabels = zeros(Int, N_cells)
    faces = Vector{Tuple{Int,Int}}(undef, N_faces)
    fNodes = Vector{NTuple{Dim,Int}}(undef, N_faces)
    fAVecs = Vector{PointT}(undef, N_faces)
    fCenters = Vector{PointT}(undef, N_faces)
    fLabels = zeros(Int,N_faces)
    nodes = Vector{PointT}(undef, N_nodes)
    nLabels = zeros(Int,N_nodes)

    index = 1
    for i_x in 1:(N_x), i_y in 1:(N_y)
        existing_points[(i_x,i_y)] = index
        nodes[index] = PointT([x_ranges[i_x],y_ranges[i_y]])
        index += 1
    end
    index = 1
    for i_x in 1:(N_x-1), i_y in 1:(N_y-1)
        i_left_up = existing_points[(i_x,i_y)]
        fNodes[index] = sort([i_left_up,existing_points[(i_x+1,i_y)]]) |> Tuple
        face_index_by_points[(fNodes[index][1],fNodes[index][2])]=index
        fAVecs[index] = R*PointT(nodes[fNodes[index][1]]-nodes[fNodes[index][2]])
        fCenters[index] = 1/2*(nodes[fNodes[index][1]]+nodes[fNodes[index][2]])
        index += 1
        fNodes[index] = sort([i_left_up,existing_points[(i_x,i_y+1)]]) |> Tuple
        face_index_by_points[(fNodes[index][1],fNodes[index][2])]=index
        fAVecs[index] = R*PointT(nodes[fNodes[index][1]]-nodes[fNodes[index][2]])
        fCenters[index] = 1/2*(nodes[fNodes[index][1]]+nodes[fNodes[index][2]])
        index += 1
    end
    for i in 1:(N_x-1)
        fNodes[index] = sort([existing_points[(i,N_y)],existing_points[(i+1,N_y)]]) |> Tuple
        face_index_by_points[(fNodes[index][1],fNodes[index][2])]=index
        fAVecs[index] = R*PointT(nodes[fNodes[index][1]]-nodes[fNodes[index][2]])
        fCenters[index] = 1/2*(nodes[fNodes[index][1]]+nodes[fNodes[index][2]])
        index += 1
    end
    for i in 1:(N_y-1)
        fNodes[index] = sort([existing_points[(N_x,i)],existing_points[(N_x,i+1)]]) |> Tuple
        face_index_by_points[(fNodes[index][1],fNodes[index][2])]=index
        fAVecs[index] = R*PointT(nodes[fNodes[index][1]]-nodes[fNodes[index][2]])
        fCenters[index] = 1/2*(nodes[fNodes[index][1]]+nodes[fNodes[index][2]])
        index += 1
    end
    index = 1
    for i_x in 1:(N_x-1), i_y in 1:(N_y-1)
        a,b,c,d =   existing_points[(i_x,i_y)],existing_points[(i_x+1,i_y)],
                    existing_points[(i_x+1,i_y+1)],existing_points[(i_x,i_y+1)]
        ab = face_index_by_points[(min(a,b),max(a,b))]
        bc = face_index_by_points[(min(c,b),max(c,b))]
        cd = face_index_by_points[(min(c,d),max(c,d))]
        da = face_index_by_points[(min(a,d),max(a,d))]
        cells[index] = [ab,bc,cd,da]
        cVols[index] = (x_ranges[i_x+1]-x_ranges[i_x])*(y_ranges[i_y+1]-y_ranges[i_y])
        cCenters[index] = PointT([1/2*(x_ranges[i_x+1]+x_ranges[i_x]),1/2*(y_ranges[i_y+1]+y_ranges[i_y])])
        cell_index_by_points[(i_x,i_y)] = index
        index+=1
    end
    for i_x in 1:(N_x-1), i_y in 1:(N_y-1)
        neighb = [(i_x+1,i_y),(i_x-1,i_y),(i_x,i_y+1),(i_x,i_y-1)]
        neighbours = neighb[[0<i[1]<N_x && 0<i[2]<N_y for i in neighb]]
        cNeighbours[cell_index_by_points[(i_x,i_y)]] = [cell_index_by_points[i] for i in neighbours]
    end
    for i_x in 1:(N_x-1), i_y in 1:(N_y-1)
        a,b,c,d =   existing_points[(i_x,i_y)],existing_points[(i_x+1,i_y)],
                    existing_points[(i_x+1,i_y+1)],existing_points[(i_x,i_y+1)]
        ab = face_index_by_points[(min(a,b),max(a,b))]
        bc = face_index_by_points[(min(c,b),max(c,b))]
        cd = face_index_by_points[(min(c,d),max(c,d))]
        da = face_index_by_points[(min(a,d),max(a,d))]
        cell = cell_index_by_points[(i_x,i_y)]
        up = i_y>1 ? cell_index_by_points[(i_x,i_y-1)] : -1
        down = i_y<N_y-1 ? cell_index_by_points[(i_x,i_y+1)] : -1
        left = i_x>1 ? cell_index_by_points[(i_x-1,i_y)] : -1
        right = i_x<N_x-1 ? cell_index_by_points[(i_x+1,i_y)] : -1
        faces[ab] = a==fNodes[ab][1] ? (up, cell) : (cell,up)
        faces[bc] = b==fNodes[bc][1] ? (right, cell) : (cell,right)
        faces[cd] = c==fNodes[cd][1] ? (down, cell) : (cell,down)
        faces[da] = d==fNodes[da][1] ? (left, cell) : (cell,left)
    end
    mesh = CFD_Mesh{Float64}(
        cells,
        cNeighbours,
        cVols,
        cCenters,
        cLabels,
        faces,
        fNodes,
        fAVecs,
        fCenters,
        fLabels,
        nodes,
        nLabels,
    )
    return mesh 
end 