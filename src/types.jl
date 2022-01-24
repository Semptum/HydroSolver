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

struct RefUnits
    ref_length::Unitful.Quantity
    ref_mass::Unitful.Quantity
    ref_time::Unitful.Quantity
    ref_current::Unitful.Quantity
    ref_temp::Unitful.Quantity
    ref_lum::Unitful.Quantity
    ref_mat::Unitful.Quantity
end

PlaneMesh = NoDim_CFD_Mesh{FloatT,SVector{2,Float64}} where FloatT
LinMesh = NoDim_CFD_Mesh{FloatT,FloatT} where FloatT


abstract type Gradient{Method} end
abstract type Laplacian{Method} end
abstract type Interpolation{Method} end
abstract type Divergence{Method} end

abstract type BoundaryConditions{T} end

struct Dirichlet{T} <: BoundaryConditions{T}
    value::T
end

struct DirichletNormal{T} <: BoundaryConditions{T}
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
