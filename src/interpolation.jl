function interpolate(mesh::NoDim_CFD_Mesh{FloatT, PointT},U,::Type{Interpolation{:Upwind}}) where {PointT,FloatT}
    N_cells = length(mesh.cells)
    N_faces = length(mesh.faces)
    interp_matrix = sparse(zeros(N_faces,N_cells))
    for i in 1:N_faces
        P,N = mesh.faces[i]
        if P==-1
            P=N
        elseif N==-1
            N=P
        end
        interp_matrix[i,U[P] ⋅ mesh.fAVecs[i] > 0  ? P : N] = 1
    end
    return interp_matrix
end

function interpolate(mesh::NoDim_CFD_Mesh{FloatT, PointT},U,::Type{Interpolation{:LinearUpwind}}) where {PointT,FloatT}
    N_cells = length(mesh.cells)
    N_faces = length(mesh.faces)
    interp_matrix = interpolate(mesh,U,Interpolation{:Upwind})
    grad = gradient(mesh,Gradient{:LeastSquares})
    face_cell = [interp_matrix[i,:].nzind[1] for i in 1:N_faces]
    d = mesh.fCenters-mesh.cCenters[face_cell]
    return sum(interp_matrix+Diagonal([i[k] for i in d])*grad[k][face_cell,:] for k in 1:length(grad))
end

function interpolate(mesh::NoDim_CFD_Mesh{FloatT, PointT},::Type{Interpolation{:Central}}) where {PointT,FloatT}
    N_cells = length(mesh.cells)
    N_faces = length(mesh.faces)
    interp_matrix = sparse(zeros(N_faces,N_cells))
    for i in 1:N_faces
        P,N = mesh.faces[i]
        if P==-1
            interp_matrix[i,N] = 1
        elseif N==-1
            interp_matrix[i,P] = 1
        else
            r_PN = mesh.cCenters[N] - mesh.cCenters[P]
            r_Pf = mesh.fCenters[i] - mesh.cCenters[P]
            α = r_Pf⋅r_PN/norm(r_PN)^2
            interp_matrix[i,[P,N]] .= [α,1-α]
        end
    end
    return interp_matrix
end