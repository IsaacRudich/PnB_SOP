include("./model/modelMain.jl")
include("./heuristics/heuristicsMain.jl")

include("./relaxed/relaxedMain.jl")
include("./restricted/restrictedMain.jl")
include("./branchAndBound/branchAndBoundMain.jl")




"""
    checkForSolution(dd)

A debug method that checks a dd for a hard coded solution in the form of Vector{Int}

# Arguments
- `dd::Vector{Vector{RelaxedSequencingNode}}`: the dd to check
"""
function checkForSolution(dd::Vector{Vector{RelaxedSequencingNode}})
    solution = [1, 5, 2, 10, 3, 6, 4, 7, 11, 8, 9, 12, 13]
    solLength = length(dd[1][1])
    startingNode = nothing

    for node in dd[1]
        if node.state == solution[solLength]
            secondCheck = true
            for e in 1:solLength
                if !(solution[e] in node.someDown)
                    secondCheck = false
                    break
                end
            end
            if secondCheck
                startingNode = node
                break
            end
        end
    end
    if isnothing(startingNode)
        return false
    end

    for i in solLength+1:length(solution)
        nextNum = solution[i]

        for arc in startingNode.arcsOut
            if arc.assignment==nextNum
                startingNode = arc.destination
                break
            end
        end

        if length(startingNode) != i
            println("_______________check for solution_____________")
            println(startingNode)
            println("_______________checked for solution_____________")
            return false
        end
    end

    return true
end
