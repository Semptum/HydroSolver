function advection(mesh::NoDim_CFD_Mesh{FloatT, PointT}, velocity;
        interp = interpolation(mesh; method=:LinearUpwind, U=velocity),
        summation = face_sum(mesh),
        BC_velocity,
        BC_quantity,
        kwargs...) where {PointT,FloatT}
    N_faces = length(mesh.faces)
    Dim = length(BC_velocity)
    interp_velocity = [apply_boundary_conditions(mesh,interp,BC_v) for BC_v in BC_velocity]
    interp_quantity = apply_boundary_conditions(mesh,interp,BC_quantity)

    v=split(velocity)
    Uf = reunite([interp_velocity[i] * v[i] for i in 1:Dim]) # The values of the velocity field on the faces: array of vectors
    Flux = [Uf[i]⋅mesh.fAVecs[i] for i in 1:N_faces] # The flux through each face : array of scalars

    M = summation * Diagonal(Flux) * interp_quantity
    return M
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


function diffusion(mesh::NoDim_CFD_Mesh{FloatT, PointT};
        interp = interpolation(mesh; method=:Central),
        grad = gradient(mesh; method=:LeastSquares),
        summation = face_sum(mesh),
        icV = inv_cell_vols(mesh),
        BC_quantity::Vector{BoundaryConditions{T}},
        kwargs...) where {T,PointT,FloatT}
    J,G,S=interp,grad,summation
    N_faces = length(mesh.faces)
    fV = split(face_vecs(mesh))
    Dim = length(mesh.cCenters[1])
    grad_on_faces = J * G

    Y = zeros(SVector{Dim,T},N_faces)
    for i in 1:N_faces
        if typeof(BC_quantity[i])==Neumann{T}
            grad_on_faces[i,:] .= [zero(PointT)]
            Y[i] = SVector{Dim,T}(
                [BC_quantity[i].gradient * mesh.fAVecs[i][k] * (mesh.faces[i][1]==-1 ? -1 : 1) 
            for k in 1:Dim])
        end
        if typeof(BC_quantity[i])==Dirichlet{T}
            grad_on_faces[i,:] .= [zero(PointT)]
            P,N = mesh.faces[i]
            if P == -1
                P=N
            end
            dPf = mesh.fCenters[i]-mesh.cCenters[P]
            grad_on_faces[i,P] = -dPf/norm(dPf)^2
            Y[i] = SVector{Dim,T}([BC_quantity[i].value * dPf[k]/norm(dPf)^2 for k in 1:Dim])
        end
    end
    grad_on_faces_tr = AffineTransform(grad_on_faces,Y)
    jg = split(grad_on_faces_tr)
    diff_op = sum(icV * S * fV[k] * jg[k] for k in 1:Dim)
    return diff_op
end