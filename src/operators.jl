function gradient(mesh::CFD_Mesh{FloatT},::Type{Gradient{:LeastSquares}}) where {FloatT}
    N = length(mesh.cells)
    ∇ = [sparse(zeros(N,N)) for i in 1:2]
    for i in 1:N
        d = reshape(reinterpret(FloatT,hcat([mesh.cCenters[neighb]-mesh.cCenters[i] for neighb in mesh.cNeighbours[i]])), (2,length(mesh.cNeighbours[i])))'
        G = d'*d
        iG = inv(G)
        grad_form = iG*d'
        for dim in 1:2
            ∇[dim][i,mesh.cNeighbours[i]] = grad_form[dim,:]
            ∇[dim][i,i] = -sum(grad_form[dim,:])
        end
    end
    return ∇
end

function convection(mesh::CFD_Mesh{FloatT}, U, InterpMethod=Interpolation{:Upwind}) where {FloatT}
    N_cells = length(mesh.cells)
    N_faces = length(mesh.faces)

    interp=interpolate(mesh,U,InterpMethod) 
    # Matrix transforming values at centroids into values at faces

    summation=sparse(zeros(N_cells,N_faces))
    # The matrix contains on each line (cell) for each column (face): 
    #       1 if the face is oriented outwards and -1 if inwards
    for i in 1:N_cells
        summation[i,mesh.cells[i]] .= [f[1]==i ? 1 : -1  for f in  mesh.faces[mesh.cells[i]]] ./ mesh.cVols[i]
    end

    Uf=interp*U # The values of the velocity field on the faces: array of vectors
    Flux = [Uf[i]⋅mesh.fAVecs[i] for i in 1:N_faces] # The flux through each face : array of scalars

    M = summation*Diagonal(Flux)*interp
    # That way, M*U_new = summation * (Flux_old .* U_new) = 
    # = [sum(Flux_old(f)*U_new(f) for f in faces(i)) for i in cells]
    return M
end

function diffusion(mesh::CFD_Mesh{FloatT}, method::Type{Laplacian{:NoCorrection}}) where FloatT
    N_cells = length(mesh.cells)
    N_faces = length(mesh.faces)
    interp_matrix = interpolate(mesh,Interpolation{:Central})
    grad = gradient(mesh,Gradient{:LeastSquares})
    summation=sparse(zeros(N_cells,N_faces))
    # The matrix contains on each line (cell) for each column (face): 
    #       1 if the face is oriented outwards and -1 if inwards
    for i in 1:N_cells
        summation[i,mesh.cells[i]] .= [f[1]==i ? 1 : -1  for f in  mesh.faces[mesh.cells[i]]] ./ mesh.cVols[i]
    end

    diff_op = summation*(Diagonal([i[1] for i in mesh.fAVecs]) * interp_matrix * grad[1] + 
        Diagonal([i[2] for i in mesh.fAVecs]) * interp_matrix * grad[2])

    return diff_op
end