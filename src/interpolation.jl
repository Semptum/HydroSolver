function interpolate(mesh::CFD_Mesh{FloatT},U,::Type{Interpolation{:Upwind}}) where {FloatT}
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

function interpolate(mesh::CFD_Mesh{FloatT},U,::Type{Interpolation{:LinearUpwind}}) where {FloatT}
    N_cells = length(mesh.cells)
    N_faces = length(mesh.faces)
    interp_matrix = interpolate(mesh,U,Interpolation{:Upwind})
    grad = gradient(mesh,Gradient{:LeastSquares})
    face_cell = [interp_matrix[i,:].nzind[1] for i in 1:N_faces]
    d = mesh.fCenters-mesh.cCenters[face_cell]
    return interp_matrix+Diagonal([i[1] for i in d])*grad[1][face_cell,:] + Diagonal([i[2] for i in d])*grad[2][face_cell,:]
end

function interpolate(mesh::CFD_Mesh{FloatT},::Type{Interpolation{:Central}}) where {FloatT}
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
            α = r_Pf⋅r_PN/norm(r_PN)
            interp_matrix[i,[P,N]] .= [α,1-α]
        end
    end
    return interp_matrix
end