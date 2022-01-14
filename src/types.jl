struct CFD_Mesh{FloatT}
    cells::Vector{Vector{Int}}      # Each element is a cell: a list of face indices
    cNeighbours::Vector{Vector{Int}}# List of neighbouring cells' indices
    cVols::Vector{FloatT}           # Volumes of the corresponding cells
    cCenters::Vector{SVector{2,FloatT}}       # Centroids of the corresponding cells
    cLabels::Vector{Int}            # Labels of the corresponding cells
    faces::Vector{Tuple{Int,Int}}   # For each face the index of the owner and neighbout cells
    fNodes::Vector{NTuple{2,Int}} # For each face the List of node indices
    fAVecs::Vector{SVector{2,FloatT}}         # Area vectors of the faces
    fCenters::Vector{SVector{2,FloatT}}       # Centroids of the faces
    fLabels::Vector{Int}            # Labels of the faces
    nodes::Vector{SVector{2,FloatT}}          # Coordinates of the nodes
    nLabels::Vector{Int}            # Labels of the nodes
end


abstract type Gradient{Method} end
abstract type Laplacian{Method} end
abstract type Interpolation{Method} end
abstract type Laplacian{Method} end

