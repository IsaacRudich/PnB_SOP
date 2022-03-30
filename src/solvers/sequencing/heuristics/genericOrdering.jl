@enum OrderingOption lexigraphic edge soporder
"""
     findSplitOrdering(model::SequencingModel,type::OrderingOption=edge)

Find a split ordering for the nodes

# Arguments
- `model::SequencingModel`: the sequencing problem to be evaluated
- `type::SOPOrderingOption`: The desired method to use for finding a split ordering, default is edge
"""
function findSplitOrdering(model::SequencingModel,type::Union{OrderingOption,Nothing}=edge)
    if type == edge || isnothing(type)
        return orderEdgeWeight(model::SequencingModel)
    elseif type == lexigraphic
        return orderLexigraphic(model)
    elseif type == soporder
        return orderSOPPrecedenceToEdge(model)
    end
end

"""
     orderLexigraphic(model::SequencingModel)

Find a variable ordering for the model based on the order of the jobs in the problem definition
Return the variable ordering

# Arguments
- `model::SequencingModel`: the sequencing problem to be evaluated
"""
function orderLexigraphic(model::SequencingModel)
    objType = valtype(model.objective.f.components)
    assignmentOrdering = Vector{objType}()
    for i in 1:length(model)
        push!(assignmentOrdering, i)
    end
    return assignmentOrdering
end

"""
     orderEdgeWeight(model::SequencingModel)

Find a variable ordering for the model based on the distances between nodes
Return the variable ordering

# Arguments
- `model::SequencingModel`: the sequencing problem to be evaluated
"""
function orderEdgeWeight(model::SequencingModel)
    #get objective values of the arcs
    objType = valtype(model.objective.f.components)
    objectiveValues = Dict{objType, Vector{objType}}()
    for (key,value) in model.objective.f.components
        for i in 1:2
            if key[i] in keys(objectiveValues)
                push!(objectiveValues[key[i]], value)
            else
                objectiveValues[key[i]] = [value]
            end
        end
    end

    #get (and sort by) averages of objective values for each node
    averages = Dict{objType, Float64}()
    for (key, value) in objectiveValues
        averages[key] = mean(value)
    end
    return sort(collect(keys(averages)), by=x->averages[x], rev=true)
end
