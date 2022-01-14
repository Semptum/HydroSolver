module HydroSolver

using StaticArrays
using LinearAlgebra
using SparseArrays

include("types.jl")
include("operators.jl")
include("interpolation.jl")
include("plotting.jl")
include("mesh_creation.jl")

end # module
