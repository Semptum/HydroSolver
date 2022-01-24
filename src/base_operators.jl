function _gradient(mesh::NoDim_CFD_Mesh{FloatT, PointT}, ::Type{Gradient{:LeastSquares}}; kwargs...) where {PointT,FloatT}
    N = length(mesh.cells)
    Dim = length(mesh.cCenters[1])
    ∇ = spzeros(PointT,N,N)
    for i in 1:N
        d = reshape(
                reinterpret(FloatT,
                    hcat([mesh.cCenters[neighb]-mesh.cCenters[i] for neighb in mesh.cNeighbours[i]])), 
            (Dim,length(mesh.cNeighbours[i])))'
        G = d'*d
        iG = inv(G)
        grad_form = iG*d'
        neighb = mesh.cNeighbours[i]
        for (j,n) in enumerate(neighb)
            if PointT <: Number
                ∇[i,n] = grad_form[1,:][j]
            else
                ∇[i,n] = PointT(grad_form[dim,:][j] for dim in 1:Dim)
            end
        end
        if PointT <: Number
            ∇[i,i] = -sum(grad_form[1,:])
        else
            ∇[i,i] = PointT(-sum(grad_form[dim,:]) for dim in 1:Dim)
        end
    end
    return ∇
end

function _laplacian(mesh::NoDim_CFD_Mesh{FloatT, PointT}, ::Type{Laplacian{:NoCorrection}};
        interp = interpolation(mesh; method=:Central),
        grad = gradient(mesh; method=:LeastSquares),
        summation = face_sum(mesh),
        kwargs...) where {PointT,FloatT}
    J,G,S=interp,grad,summation
    fV = split(face_vecs(mesh))
    g = split(G)
    Dim = length(mesh.cCenters[1])
    diff_op = sum(S * fV[k] * J * g[k] for k in 1:Dim)
    return diff_op
end

function _divergence(mesh::NoDim_CFD_Mesh{FloatT, PointT}, ::Type{Divergence{:FinVol}}; 
        interp = interpolation(mesh; method=:Central),
        summation = face_sum(mesh),
        fV = face_vecs(mesh),
        icV = inv_cell_vols(mesh),
        kwargs...) where {PointT,FloatT}
    fv = split(fV)
    Dim = length(fv)
    return reunite([icV*summation*fv[i]*interp for i in 1:Dim])
end


function face_sum(mesh::NoDim_CFD_Mesh{FloatT, PointT}; kwargs...) where {PointT,FloatT}
    N_cells = length(mesh.cells)
    N_faces = length(mesh.faces)
    summation=spzeros(N_cells,N_faces)
    # The matrix contains on each line (cell) for each column (face): 
    #       1 if the face is oriented outwards and -1 if inwards
    for i in 1:N_cells
        summation[i,mesh.cells[i]] .= [f[1]==i ? 1 : -1  for f in  mesh.faces[mesh.cells[i]]]
    end
    summation
end

function cell_vols(mesh)
    Diagonal(mesh.cVols)
end

function inv_cell_vols(mesh)
    Diagonal(1 ./ mesh.cVols)
end

function face_vecs(mesh)
    Diagonal(mesh.fAVecs)
end


function gradient(mesh::NoDim_CFD_Mesh{FloatT, PointT}; method=:LeastSquares, kwargs...) where {PointT,FloatT}
    _gradient(mesh, Gradient{method}; kwargs...)
end

function laplacian(mesh::NoDim_CFD_Mesh{FloatT, PointT}; method=:NoCorrection, kwargs...) where {PointT,FloatT}
    _laplacian(mesh, Laplacian{method}; kwargs...)
end

function divergence(mesh::NoDim_CFD_Mesh{FloatT, PointT}; method=:FinVol, kwargs...) where {PointT,FloatT}
    _divergence(mesh, Divergence{method}; kwargs...)
end
