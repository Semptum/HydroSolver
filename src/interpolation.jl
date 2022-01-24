function interpolation(mesh::NoDim_CFD_Mesh{FloatT, PointT}, ::Type{Interpolation{:Upwind}}; U) where {PointT,FloatT}
    N_cells = length(mesh.cells)
    N_faces = length(mesh.faces)
    interp_matrix = spzeros(N_faces,N_cells)
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

function interpolation(mesh::NoDim_CFD_Mesh{FloatT, PointT}, ::Type{Interpolation{:LinearUpwind}};
        U,
        interp_matrix = interpolation(mesh,method=:Upwind, U=U),
        grad = gradient(mesh,method=:LeastSquares),
        kwargs...) where {PointT,FloatT}
    N_faces = length(mesh.faces)
    Dim = length(mesh.cCenters[1])
    face_cell = [interp_matrix[i,:].nzind[1] for i in 1:N_faces]
    d = mesh.fCenters-mesh.cCenters[face_cell]
    G = split(grad)
    return sum(
        interp_matrix + Diagonal([i[k] for i in d]) * G[k][face_cell,:]
    for k in 1:Dim)
end

function interpolation(mesh::NoDim_CFD_Mesh{FloatT, PointT}, ::Type{Interpolation{:Central}}) where {PointT,FloatT}
    N_cells = length(mesh.cells)
    N_faces = length(mesh.faces)
    interp_matrix = spzeros(N_faces,N_cells)
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

function rhiechow(mesh::NoDim_CFD_Mesh{FloatT, PointT};
        interp = interpolation(mesh,method=:Central),
        grad = gradient(mesh,method=:LeastSquares),
        kwargs...) where {PointT,FloatT}
    N_cells = length(mesh.cells)
    N_faces = length(mesh.faces)    
    basic_grad = spzeros(PointT,N_faces,N_cells)
    interp_grad = interp*grad
    for i in 1:N_faces
        P,N = mesh.faces[i]
        if P==-1 || N==-1
            basic_grad[i,:] .= interp_grad[i,:]
        else
            r_PN = mesh.cCenters[N] - mesh.cCenters[P]
            basic_grad[i,N] = r_PN/norm(r_PN)^2
            basic_grad[i,P] = -basic_grad[i,N]
        end
    end   
    return basic_grad-interp_grad
end

function interpolation(mesh; method=:LinearUpwind, kwargs...)
    interpolation(mesh, Interpolation{method}; kwargs...)
end