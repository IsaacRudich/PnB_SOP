#module DDForDP
using Statistics
using UUIDs

include("./utilities/utilitiesMain.jl")
include("./intervalTrees/treesMain.jl")
include("./model/modelMain.jl")

@enum RRBSetting sop
include("./solvers/solversMain.jl")
#end #module

include("./io/ioMain.jl")
include("./examples/examples.jl")


function testRelaxedBound()
    model = testModelMain()

    (bound, path) = calculateRelaxedBound(model, 2)

    println(bound)
    for (key,value) in path
        println("   $key => $value")
    end
end
