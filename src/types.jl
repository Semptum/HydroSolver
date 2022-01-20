struct NoDim_CFD_Mesh{FloatT, PointT}
    cells::Vector{Vector{Int}}      # Each element is a cell: a list of face indices
    cNeighbours::Vector{Vector{Int}}# List of neighbouring cells' indices
    cVols::Vector{FloatT}           # Volumes of the corresponding cells
    cCenters::Vector{PointT}       # Centroids of the corresponding cells
    cLabels::Vector{Int}            # Labels of the corresponding cells
    faces::Vector{SVector{2,Int}}   # For each face the index of the owner and neighbout cells
    fNodes::Vector{Vector{Int}} # For each face the List of node indices
    fAVecs::Vector{PointT}         # Area vectors of the faces
    fCenters::Vector{PointT}       # Centroids of the faces
    fLabels::Vector{Int}            # Labels of the faces
    nodes::Vector{PointT}          # Coordinates of the nodes
    nLabels::Vector{Int}            # Labels of the nodes
end

CFD_Mesh = NoDim_CFD_Mesh{FloatT,SVector{2,Float64}} where FloatT
LinMesh = NoDim_CFD_Mesh{FloatT,FloatT} where FloatT

struct CFD_Solution{FloatT}
    mesh::CFD_Mesh{FloatT}
end


abstract type Gradient{Method} end
abstract type Laplacian{Method} end
abstract type Interpolation{Method} end

abstract type BoundaryConditions{T} end

struct Dirichlet{T} <: BoundaryConditions{T}
    value::T
end

struct Neumann{T} <: BoundaryConditions{T}
    gradient::T
end

struct None{T} <: BoundaryConditions{T} end

struct AffineTransform{MatrixT,VectorT}
    A::MatrixT
    B::VectorT
end

function +(trX::AffineTransform, trY::AffineTransform)
    A = trX.A + trY.A
    B = trX.B + trY.B
    MatrixT = typeof(A)
    VectorT = typeof(B)
    return AffineTransform{MatrixT,VectorT}(A,B)
end

function +(trX::AffineTransform, M::AbstractMatrix)
    A = trX.A + M
    B = trX.B
    MatrixT = typeof(A)
    VectorT = typeof(B)
    return AffineTransform{MatrixT,VectorT}(A,B)
end

function +(M::AbstractMatrix, trX::AffineTransform)
    trX+M
end

function -(M::AbstractMatrix, trX::AffineTransform)
    -trX+M
end

function -(trX::AffineTransform, M::AbstractMatrix)
    trX + -1*M
end

function -(trX::AffineTransform)
    trX *-1
end

function *(trX::AffineTransform, trY::AffineTransform)
    A = trX.A * trY.A
    B = trX.A * trY.B + trX.B
    MatrixT = typeof(A)
    VectorT = typeof(B)
    return AffineTransform{MatrixT,VectorT}(A,B)
end

function *(M::AbstractMatrix, trX::AffineTransform)
    A = M * trX.A
    B = M * trX.B
    MatrixT = typeof(A)
    VectorT = typeof(B)
    return AffineTransform{MatrixT,VectorT}(A,B)
end

function *(trX::AffineTransform, M::AbstractMatrix)
    A = trX.A * M
    B = trX.B
    MatrixT = typeof(A)
    VectorT = typeof(B)
    return AffineTransform{MatrixT,VectorT}(A,B)
end

function *(trX::AffineTransform, X::AbstractVector)
    return trX.A * X + trX.B
end


function *(trX::AffineTransform{MT,VT}, X::Number) where {VT,MT}
    return AffineTransform{MT,VT}(trX.A*X,trX.B*X)
end

function *(X::Number, trX::AffineTransform{MT,VT}) where {VT,MT}
    return AffineTransform{MT,VT}(trX.A*X,trX.B*X)
end


function *(x::SVector{2,T}, y::SVector{2,T}) where T
    x⋅y
end