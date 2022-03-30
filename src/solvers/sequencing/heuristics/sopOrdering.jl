"""
     orderSOPJobsPrecedenceToEdge(model::SequencingModel)

Find a variable ordering for the model based on the edge weights and precedence graph
Return the variable ordering

# Arguments
- `model::SequencingModel`: the sequencing problem to be evaluated
"""
function orderSOPPrecedenceToEdge(model::SequencingModel)
    objType = valtype(model.objective.f.components)
    sortedKeys = orderEdgeWeight(model)
    #put jobs into the order by precedence and then average edge weight
    finalOrder = Vector{objType}()
    followerNumbers = getPrecedenceNumbers(model.precedenceP2F)
    maxLength = maximum(values(followerNumbers))
    while maxLength>0
        for i in 1:length(sortedKeys)
            if sortedKeys[i] in keys(followerNumbers) && followerNumbers[sortedKeys[i]]==maxLength
                push!(finalOrder, sortedKeys[i])
            end
        end
        maxLength = maxLength - 1
    end
    for i in 1:length(sortedKeys)
        if !(sortedKeys[i] in finalOrder)
            push!(finalOrder, sortedKeys[i])
        end
    end
    return finalOrder
end
