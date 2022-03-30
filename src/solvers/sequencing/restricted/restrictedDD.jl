"""
    restrictedSequencingDD(model::SequencingModel, maxWidth::Int;rootNode::Union{RestrictedSequencingNode,RelaxedSequencingNode,RelaxedStateSequencingNode,Nothing}=nothing,bestKnownSolution::Union{Vector{Int},Nothing}=nothing,bestKnownValue::Union{T,Nothing}=nothing, rrbSetting::Union{RRBSetting,Nothing}=nothing, rrbPacket=nothing, rootDomain::BitVector=BitVector(),debugOn::Bool=false)where{T<:Real}

A quick restricted decision diagram for sequencing problems of fixed length
Returns: bestKnownSolution (Vector{Int}), bestKnownValue (Real), exact {Bool}
'exact' is true if the dd is exact 

# Arguments
- `model::SequencingModel`: The sequencing problem
- `maxWidth::Int`: The max width allowed in a decision diagram layer
- `rootNode::Union{RestrictedSequencingNode,RelaxedSequencingNode,RelaxedStateSequencingNode,Nothing}`: The starting node for the DD
- `bestKnownSolution::Union{Vector{Int},Nothing}`: An optional parameter that gives the DD a starting solution
- `bestKnownValue::Union{T,Nothing}`: An optional parameter that gives the DD a starting bound
- `rrbSetting::Union{RRBSetting,Nothing}`: An optional parameter that chooses an rough relaxed bounding method
- `rrbPacket`: Data to be passed to the rrbFunction
- `rootDomain::BitVector`: The domain of the root node if it is a RelaxedStateSequencingNode
- `debugOn::Bool`: An optional parameter that can be used to turn on debug statements in other functions
"""
function restrictedSequencingDD(model::SequencingModel, maxWidth::Int;rootNode::Union{RestrictedSequencingNode,RelaxedSequencingNode,RelaxedStateSequencingNode,Nothing}=nothing,bestKnownSolution::Union{Vector{Int},Nothing}=nothing,bestKnownValue::Union{T,Nothing}=nothing, rrbSetting::Union{RRBSetting,Nothing}=nothing, rrbPacket=nothing, rootDomain::BitVector=BitVector(),debugOn::Bool=false)where{T<:Real}
    doSecondLayerDomainCheck = false
    if isnothing(rootNode)
        rootNode = createRestrictedRootNode(model)
        queue = Vector{RestrictedSequencingNode}([rootNode])
        pathLength = length(rootNode.path)
    elseif isa(rootNode,RelaxedSequencingNode)
        pathLength = length(rootNode)
        rootDomain = falses(length(model))
        for arc in rootNode.arcsOut
            rootDomain[arc.assignment] = true
        end
        queue = Vector{RestrictedSequencingNode}([RestrictedSequencingNode(copy(rootNode.bestSeq),rootNode.allDown, rootNode.lengthToRoot,rootDomain)])
        if model.precedenceConstraints
            doSecondLayerDomainCheck = true
        end
    elseif isa(rootNode,RelaxedStateSequencingNode)
        pathLength = length(rootNode)
        queue = Vector{RestrictedSequencingNode}([RestrictedSequencingNode(copy(rootNode.bestSeq),rootNode.allDown, rootNode.lengthToRoot,rootDomain)])
        if model.precedenceConstraints
            doSecondLayerDomainCheck = true
        end
    else
        queue = Vector{RestrictedSequencingNode}([rootNode])
        pathLength = length(rootNode.path)
    end

    #track exactness
    exact = true


    #iterate down through the layers
    for layerIndex in 1:length(model)-pathLength
        if layerIndex == 2 && doSecondLayerDomainCheck
            secondLayerDomainCheck(model, queue)
        end
        #initialize array for new layer
        newLayer = Vector{RestrictedSequencingNode}()
        #for each node in the queue
        for node in queue
            if length(node.path)>0
                state = getState(node)
            else
                state = nothing
            end
            #for each unvisited location of the node
            for i in 1:length(model)
                #if i is in the domain add the new node
                if inDomain(node, i)
                    newNode = makeDecision(model, node, state, i)
                                            
                    if !isnothing(rrbSetting) && layerIndex+pathLength<length(model)-2
                        if  !isBetterSolutionValue(model, sequencingRRB(rrbSetting, newNode, rrbPacket), bestKnownValue)
                            continue
                        end
                    end

                    #if worth adding, add it
                    if isnothing(bestKnownValue) || model.objective.type==maximization
                        addToSortedLayer!(model, newLayer, newNode)
                    elseif isBetterSolutionValue(model, newNode.value, bestKnownValue)
                        addToSortedLayer!(model, newLayer, newNode)
                    end
                end
            end#end domain loop
            if length(newLayer)>maxWidth
                if exact
                    exact = false
                end
                newLayer = newLayer[1:maxWidth]
            end
        end#node in queue loop
        #trim queue
        if length(newLayer)>maxWidth
            if exact
                exact = false
            end
            queue = newLayer[1:maxWidth]
        elseif isempty(newLayer)
            return bestKnownSolution, bestKnownValue, exact
        else
            queue = newLayer
        end
    end#iterate over layers loop

    #update result
    if isnothing(bestKnownValue) || isBetterSolutionValue(model, queue[1].value, bestKnownValue)
        bestKnownSolution = queue[1].path
        bestKnownValue = queue[1].value
    end
    return bestKnownSolution, bestKnownValue, exact
end

"""
    makeDecision(model::SequencingModel, rootNode::RestrictedSequencingNode, state::Union{Nothing, Int}, decision::Int)

A method that generates the child node resulting from adding an arc to an existing node

# Arguments
- `model::SequencingModel`: The Sequencing problem
- `rootNode::RestrictedSequencingNode`: The starting node for the decision
- `state::Union{Nothing, Int}`: The last node visited
- `decision::Int`: The node being visited
"""
function makeDecision(model::SequencingModel, rootNode::RestrictedSequencingNode, state::Union{Nothing, Int}, decision::Int)
    #update path
    newPath = copy(rootNode.path)
    push!(newPath, decision)

    #update visited
    newVisited = copy(rootNode.visited)
    newVisited[decision] = true

    #update domain
    newDomain = copy(rootNode.domain)
    newDomain[decision] = false
    #handle precedence constraints
    if model.precedenceConstraints
        if haskey(model.precedenceP2F,decision)
            for i in model.precedenceP2F[decision]
                #check each of the relevant precedence constraints
                #if satisfied, put in the domain
                if haskey(model.precedenceF2P,i)
                    if checkPrecedence(newVisited,model.precedenceF2P[i])
                        newDomain[i] = true
                    end
                end
            end
        end
    end

    #update value
    if isnothing(state)
        return RestrictedSequencingNode(newPath, newVisited, rootNode.value, newDomain)
    else
        return RestrictedSequencingNode(newPath, newVisited, rootNode.value + evaluateDecision(model.objective, state, decision), newDomain)
    end
end

"""
    addToSortedLayer!(model::SequencingModel, layer::Vector{RestrictedSequencingNode}, newNode::RestrictedSequencingNode

A method that generates the child node resulting from adding an arc to an existing node

# Arguments
- `model::SequencingModel`: The Sequencing problem
- `layer::Vector{RestrictedSequencingNode}`: The list to add to
- `newNode::RestrictedSequencingNode`: the node to insert
"""
function addToSortedLayer!(model::SequencingModel, layer::Vector{RestrictedSequencingNode}, newNode::RestrictedSequencingNode)
    if isempty(layer)
        push!(layer, newNode)
    else
        insert!(
            layer,
            searchsorted(layer, newNode, by = x -> valueSequencingNode(model, x)).start,
            newNode
        )
    end
end

"""
    valueSequencingNode(model::CVRPModel, node::RestrictedCVRPNode)

Return a quick heuristic value of a node, lower is better

# Arguments
- `model::SequencingModel`: The sequencing problem
- `newNode::RestrictedSequencingNode`: the node to value
"""
function valueSequencingNode(model::SequencingModel, node::RestrictedSequencingNode)
    if model.objective.type == minimization
        return node.value
    else
        return (-1 * node.value)
    end
end

"""
    secondLayerDomainCheck(model::SequencingModel, queue::Vector{RestrictedSequencingNode})

A method that fixes the domain of the first nodes generated from a relaxed node

# Arguments
- `model::SequencingModel`: The Sequencing problem
- `queue::Vector{RestrictedSequencingNode}`: The list of nodes to check
"""
function secondLayerDomainCheck(model::SequencingModel, queue::Vector{RestrictedSequencingNode})
    for node in queue
        for i in 1:length(node.domain)
            if !node.visited[i] && !node.domain[i] && hasKey(model.objective, getState(node), i)
                #handle precedence constraints
                if haskey(model.precedenceF2P,i)
                    if checkPrecedence(node.visited,model.precedenceF2P[i])
                        node.domain[i] = true
                    end
                else
                    node.domain[i] = true
                end
            end
        end
    end
end