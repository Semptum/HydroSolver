module HydroSolver

using StaticArrays
using LinearAlgebra
using SparseArrays

using Makie

import Base.+
import Base.*
import Base.-

include("types.jl")
include("operators.jl")
include("interpolation.jl")
include("plotting.jl")
include("mesh_creation.jl")
include("utils.jl")

end # module
