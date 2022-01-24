module HydroSolver

using StaticArrays
using LinearAlgebra
using SparseArrays
using Statistics
using Makie
using Unitful

import Base.+
import Base.*
import Base.-

export inrefunit
export mmgs

include("types.jl")
include("operators.jl")
include("interpolation.jl")
include("plotting.jl")
include("mesh_creation.jl")
include("utils.jl")
include("units.jl")
include("base_operators.jl")

end # module
