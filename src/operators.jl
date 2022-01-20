function gradient(mesh::NoDim_CFD_Mesh{FloatT, PointT},::Type{Gradient{:LeastSquares}}) where {PointT,FloatT}
    N = length(mesh.cells)
    Dim = length(mesh.cCenters[1])
    ∇ = [sparse(zeros(N,N)) for i in 1:Dim]
    for i in 1:N
        d = reshape(
                reinterpret(FloatT,
                    hcat([mesh.cCenters[neighb]-mesh.cCenters[i] for neighb in mesh.cNeighbours[i]])), 
            (Dim,length(mesh.cNeighbours[i])))'
        G = d'*d
        iG = inv(G)
        grad_form = iG*d'
        for dim in 1:Dim
            ∇[dim][i,mesh.cNeighbours[i]] = grad_form[dim,:]
            ∇[dim][i,i] = -sum(grad_form[dim,:])
        end
    end
    return ∇
end

function convection(mesh::NoDim_CFD_Mesh{FloatT, PointT}, U, BC, InterpMethod=Interpolation{:Upwind}) where {PointT,FloatT}
    N_cells = length(mesh.cells)
    N_faces = length(mesh.faces)

    interp = interpolate(mesh,U,InterpMethod)
    # Matrix transforming values at centroids into values at faces

    summation=sparse(zeros(N_cells,N_faces))
    # The matrix contains on each line (cell) for each column (face): 
    #       1 if the face is oriented outwards and -1 if inwards
    for i in 1:N_cells
        summation[i,mesh.cells[i]] .= [f[1]==i ? 1 : -1  for f in  mesh.faces[mesh.cells[i]]] ./ mesh.cVols[i]
    end
    interp_with_BC = apply_boundary_conditions(mesh,interp,BC)

    Uf=interp_with_BC*U # The values of the velocity field on the faces: array of vectors
    Flux = [Uf[i]⋅mesh.fAVecs[i] for i in 1:N_faces] # The flux through each face : array of scalars

    M = summation*Diagonal(Flux)*interp_with_BC
    # That way, M*U_new = summation * (Flux_old .* U_new) = 
    # = [sum(Flux_old(f)*U_new(f) for f in faces(i)) for i in cells]
    return M
end

function diffusion(mesh::NoDim_CFD_Mesh{FloatT, PointT}, method::Type{Laplacian{:NoCorrection}}) where {PointT,FloatT}
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

    diff_op = summation*sum(Diagonal([i[k] for i in mesh.fAVecs]) * interp_matrix * grad[k] for k in 1:length(mesh.cCenters[1]))
    return diff_op
end


function summation(mesh::NoDim_CFD_Mesh{FloatT, PointT}) where {PointT,FloatT}
    N_cells = length(mesh.cells)
    N_faces = length(mesh.faces)
    summation=sparse(zeros(N_cells,N_faces))
    # The matrix contains on each line (cell) for each column (face): 
    #       1 if the face is oriented outwards and -1 if inwards
    for i in 1:N_cells
        summation[i,mesh.cells[i]] .= [f[1]==i ? 1 : -1  for f in  mesh.faces[mesh.cells[i]]] ./ mesh.cVols[i]
    end
    summation
end


function apply_boundary_conditions(mesh,interp_matrix,BC::Vector{BoundaryConditions{T}}) where T
    interp_matrix_new = copy(interp_matrix)
    N_faces = length(BC)
    Y = zeros(T,N_faces)
    for i in 1:N_faces
        if typeof(BC[i])==Dirichlet{T}
            interp_matrix_new[i,:] .= 0
            Y[i] = BC[i].value
        end
        if typeof(BC[i])==Neumann{T}
            interp_matrix_new[i,:] .= 0
            P,N = mesh.faces[i]
            n = mesh.fAVecs[i]
            n /= norm(n)
            if P == -1
                P=N
                n = -n
            end
            dPf = mesh.fCenters[i]-mesh.cCenters[P]

            Y[i] = BC[i].gradient * (dPf⋅n)
            interp_matrix_new[i,P] = 1.0
        end
    end
    return AffineTransform(interp_matrix_new,Y)
end