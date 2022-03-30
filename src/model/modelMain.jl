@enum ObjectiveType minimization=0 maximization=1 allDiff=2
abstract type ObjectiveFunction end

struct Range{T<:Real}
    max     ::T
    min     ::T
end

include("./sequence.jl")
include("./geometry/geometryMain.jl")
include("./finalModels/finalModelsMain.jl")

abstract type Node end
include("./nodes/nodesMain.jl")
