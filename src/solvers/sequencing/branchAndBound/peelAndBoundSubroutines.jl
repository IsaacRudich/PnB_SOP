"""
    peelNode!(model::SequencingModel, node::RelaxedSequencingNode, arcs::Vector{ExactArc};bestKnownValue::Union{T,Nothing}=nothing,rrbSetting::Union{RRBSetting,Nothing},rrbPacket)where{T<:Real}

Peels a node, like splitting a node, but the sorting of arcs is pre-determined
Returns peeledNode (RelaxedSequencingNode), the new node

# Arguments
- `model::SequencingModel`: The sequencing problem
- `node::RelaxedSequencingNode`: The node to peel
- `arcs::Vector{ExactArc}`: The set of 'in' arcs belonging to the new node
- `bestKnownValue::Union{T,Nothing}`: The value of the best known feasible solution
- `rrbSetting::Union{RRBSetting,Nothing}`: A rough relaxed bounding method, or lack thereof
- `rrbPacket`: Optional data to be passed to the rrbFunction

"""
function peelNode!(model::SequencingModel, node::RelaxedSequencingNode, arcs::Vector{ExactArc};bestKnownValue::Union{T,Nothing}=nothing,rrbSetting::Union{RRBSetting,Nothing},rrbPacket)where{T<:Real}
    peeledNode = RelaxedSequencingNode(
        node.length, #bestPath
        node.state, #current state
        copy(node.allUp), #allUp
        BitVector(), #allDown
        copy(node.someUp), #someUp
        BitVector(), #someDown
        node.lengthToRoot, #to root
        node.lengthToTerminal, #to terminal
        Vector{Int}(),
        node.exact, #exact
        Vector{ExactArc}(), #arcsin
        Vector{ExactArc}() #arcsout
    )

    #add new arcs out
    for arc in node.arcsOut
        newArc = ExactArc(peeledNode, arc.destination, arc.assignment, arc.value)
        push!(peeledNode.arcsOut, newArc)
        push!(arc.destination.arcsIn, newArc)
    end

    #fix the destination of the in arc for the new node
    for arc in arcs
        newArc = ExactArc(arc.origin, peeledNode, arc.assignment, arc.value)
        push!(peeledNode.arcsIn, newArc)
        push!(arc.origin.arcsOut, newArc)
    end

    #get rid of the extra incoming arcs
    deleteArcs!(arcs)

    updateNodeDownVariables!(peeledNode,model)
    filterArcs!(peeledNode.arcsOut,model,bestKnownValue=bestKnownValue, rrbSetting=rrbSetting, rrbPacket=rrbPacket)

    if isempty(node.arcsIn)
        deleteArcs!(node.arcsOut)
    elseif isempty(node.arcsOut) && length(node) != length(model)
        deleteArcs!(node.arcsIn)
    else
        updateNodeDownVariables!(node,model)
        filterArcs!(node.arcsOut,model,bestKnownValue=bestKnownValue, rrbSetting=rrbSetting, rrbPacket=rrbPacket)
    end

    return peeledNode
end


"""
peelDD!(
    model::SequencingModel,
    maxWidth::T,
    dd::Vector{Vector{RelaxedSequencingNode}},
    relaxedSolution::Vector{T},
    frontierNode::RelaxedSequencingNode,
    fNodeIndex::U,
    splitOrder::Vector{T},
    bestKnownValue::Union{U,Nothing},
    rrbSetting::Union{RRBSetting,Nothing},
    rrbPacket,
    loggingOn::Bool,
    debugOn::Bool
)where{T<:Int,U<:Real}

Peel a diagram given the node to start from, and the diagram
Returns (dd, bound, relaxedSol, ddExact), (newDD, peeledBound, peeledRelaxedSol, peeledExact), bestKnownSolution, bestKnownValue
(dd, bound, relaxedSol, ddExact): is a Tuple{Vector{Vector{RelaxedSequencingNode}}, Real, Vector{Int}, Bool} storing all the updated info for 'dd'
(newDD, peeledBound, peeledRelaxedSol, peeledExact): is also a Tuple{Vector{Vector{RelaxedSequencingNode}}, Real, Vector{Int}, Bool} storing all the info for the peeled dd
bestKnownSolution, bestKnownValue: are Vector{Int}, Real

# Arguments
- `model::SequencingModel`: The sequencing problem
- `maxWidth::Int`: The max width allowed in a decision diagram layer
- `dd::Vector{Vector{RelaxedSequencingNode}}`: The starting diagram
- `relaxedSolution::::Vector{T}`: The best path through 'dd'
- `frontierNode::RelaxedSequencingNode`: The node start from
- `fNodeIndex::U`: The index of the layer of 'dd' that 'frontierNode' is on
- `splitOrder::Vector{T}`: A heuristic ordering of nodes
- `bestKnownValue::Union{T,Nothing}`: The value of the best known feasible solution
- `rrbSetting::Union{RRBSetting,Nothing}`: A rough relaxed bounding method, or lack thereof
- `rrbPacket`: Optional data to be passed to the rrbFunction
- `loggingOn::Bool`: An optional parameter that can be used to turn off logging to the console and local file
- `debugOn::Bool`: An optional parameter that can be used to turn on debug statements in other functions

"""
function peelDD!(
    model::SequencingModel,
    maxWidth::T,
    dd::Vector{Vector{RelaxedSequencingNode}},
    relaxedSolution::Vector{T},
    frontierNode::RelaxedSequencingNode,
    fNodeIndex::U,
    splitOrder::Vector{T},
    bestKnownValue::Union{U,Nothing},
    rrbSetting::Union{RRBSetting,Nothing},
    rrbPacket,
    loggingOn::Bool,
    debugOn::Bool
)where{T<:Int,U<:Real}
    cleanRelaxedDD(dd) = relaxedSequencingDD(model, maxWidth,splitOrder, dd = dd, debugOn=debugOn,bestKnownSolution=bestKnownSolution, bestKnownValue=bestKnownValue, rrbSetting=rrbSetting, rrbPacket=rrbPacket)
    bestKnownSolution = Vector{Int}()
    filter!(x->x!=frontierNode, dd[fNodeIndex])
    deleteArcs!(frontierNode.arcsIn)
    newDD = Vector{Vector{RelaxedSequencingNode}}()
    push!(newDD, [frontierNode])
    for i in 1:length(relaxedSolution)-length(frontierNode)
        arcMap = Dict{RelaxedSequencingNode, Vector{ExactArc}}()
        #create a map of child nodes to the arcs that should go with them
        for node in newDD[i]
            for arc in node.arcsOut
                if haskey(arcMap, arc.destination)
                    push!(arcMap[arc.destination], arc)
                else
                    arcMap[arc.destination] = [arc]
                end
            end
        end
        #iterate over each child node of the current layer and peel it
        push!(newDD, Vector{RelaxedSequencingNode}())
        for key in keys(arcMap)
            newNode = peelNode!(model, key, arcMap[key],bestKnownValue=bestKnownValue, rrbSetting=rrbSetting, rrbPacket=rrbPacket)
            #check the peeled node to make sure it is worth adding
            if !isempty(newNode.arcsOut) || length(model)==length(newNode)
                push!(newDD[i+1], newNode)
            else #delete in arcs so parent nodes may be removed in subsequent passes
                deleteArcs!(newNode.arcsIn)
            end
            #check the old node to make sure it doesnt need to be removed
            if isempty(key.arcsOut) && (i+fNodeIndex)<length(dd)
                #delete in arcs so parent nodes may be removed in subsequent passes
                deleteArcs!(key.arcsIn)
                filter!(x->x!=key, dd[i+fNodeIndex])
            elseif isempty(key.arcsIn)
                #delete in arcs so parent nodes may be removed in subsequent passes
                deleteArcs!(key.arcsOut)
                filter!(x->x!=key, dd[i+fNodeIndex])
            end
        end

        for node in dd[i+fNodeIndex]
            if isempty(node.arcsIn)
                deleteArcs!(node.arcsOut)
            end
        end
        filter!(x->length(x.arcsIn)!=0, dd[i+fNodeIndex])
        for node in dd[i+fNodeIndex]
            updateNodeDownVariables!(node,model)
        end
    end#iterate through the layers
    #check and upddate the two DDs
    if !isempty(last(dd))
        bottomUpUpdate(model, dd,bestKnownValue)
        bound = last(dd)[1].lengthToRoot
        relaxedSol = last(dd)[1].bestSeq
        bestKnownSolution, bestKnownValue = checkFeasibleTerminalSolutions!(model, last(dd)[1],bestKnownValue=bestKnownValue,bestKnownSolution=bestKnownSolution)
        if !isnothing(bound) && isBetterSolutionValue(model, bound, bestKnownValue)
            ddExact = false
        else
            dd = nothing
            bound = bestKnownValue
            relaxedSol = nothing
            ddExact = true
        end
    else
        dd = nothing
        bound = bestKnownValue
        relaxedSol = nothing
        ddExact = true
    end
    if !isempty(last(newDD))
        bottomUpUpdate(model, newDD, bestKnownValue)
        peeledBound = last(newDD)[1].lengthToRoot
        peeledRelaxedSol = last(newDD)[1].bestSeq
        bestKnownSolution, bestKnownValue = checkFeasibleTerminalSolutions!(model, last(newDD)[1],bestKnownValue=bestKnownValue,bestKnownSolution=bestKnownSolution)
        if !isnothing(peeledBound) && isBetterSolutionValue(model, peeledBound, bestKnownValue)
            newDD, peeledBound, peeledRelaxedSol, bestKnownSolution, bestKnownValue, peeledExact = cleanRelaxedDD(newDD)
            if peeledExact
                newDD = Vector{Vector{RelaxedSequencingNode}}()
            end
        else
            newDD = nothing
            peeledBound = bestKnownValue
            peeledRelaxedSol = nothing
            peeledExact = true
        end
    else
        newDD = nothing
        peeledBound = bestKnownValue
        peeledRelaxedSol = nothing
        peeledExact = true
    end

    return (dd, bound, relaxedSol, ddExact), (newDD, peeledBound, peeledRelaxedSol, peeledExact), bestKnownSolution, bestKnownValue
end

"""
peelAndRemove!(
    model::SequencingModel,
    dd::Vector{Vector{RelaxedSequencingNode}},
    frontierNode::RelaxedSequencingNode,
    fNodeIndex::U,
    bestKnownValue::Union{U,Nothing},
    rrbSetting::Union{RRBSetting,Nothing},
    rrbPacket,
    loggingOn::Bool,
    debugOn::Bool
)where{T<:Int,U<:Real}

A version of the peel process that deletes the peeled dd as it is created
Returns (dd, bound, relaxedSol, ddExact), bestKnownSolution, bestKnownValue
(dd, bound, relaxedSol, ddExact): is a Tuple{Vector{Vector{RelaxedSequencingNode}}, Real, Vector{Int}, Bool} storing all the updated info for 'dd'
bestKnownSolution, bestKnownValue: are Vector{Int}, Real

# Arguments
- `model::SequencingModel`: The sequencing problem
- `dd::Vector{Vector{RelaxedSequencingNode}}`: The starting diagram
- `frontierNode::RelaxedSequencingNode`: The node start from
- `fNodeIndex::U`: The index of the layer of 'dd' that 'frontierNode' is on
- `bestKnownValue::Union{T,Nothing}`: The value of the best known feasible solution
- `rrbSetting::Union{RRBSetting,Nothing}`: A rough relaxed bounding method, or lack thereof
- `rrbPacket`: Optional data to be passed to the rrbFunction
- `loggingOn::Bool`: An optional parameter that can be used to turn off logging to the console and local file
- `debugOn::Bool`: An optional parameter that can be used to turn on debug statements in other functions

"""
function peelAndRemove!(
    model::SequencingModel,
    dd::Vector{Vector{RelaxedSequencingNode}},
    frontierNode::RelaxedSequencingNode,
    fNodeIndex::U,
    bestKnownValue::Union{U,Nothing},
    rrbSetting::Union{RRBSetting,Nothing},
    rrbPacket,
    loggingOn::Bool,
    debugOn::Bool
)where{T<:Int,U<:Real}
    bestKnownSolution = Vector{Int}()

    #remove the node from the dd
    deleteArcs!(frontierNode.arcsIn)
    deleteArcs!(frontierNode.arcsOut)

    filter!(x->x!=frontierNode, dd[fNodeIndex])

    #iterate downwards
    for i in (fNodeIndex+1):length(dd)
        for node in dd[i]
            if isempty(node.arcsIn)
                deleteArcs!(node.arcsOut)
            end
        end
        filter!(x->length(x.arcsIn)!=0, dd[i])
        for node in dd[i]
            updateNodeDownVariables!(node,model)
            filterArcs!(node.arcsOut,model;bestKnownValue=bestKnownValue, rrbSetting=rrbSetting, rrbPacket=rrbPacket)
        end
    end

    #update upwards values
    if !isempty(last(dd))
        bottomUpUpdate(model, dd,bestKnownValue)
        
        bound = last(dd)[1].lengthToRoot
        relaxedSol = last(dd)[1].bestSeq
        bestKnownSolution, bestKnownValue = checkFeasibleTerminalSolutions!(model, last(dd)[1],bestKnownValue=bestKnownValue,bestKnownSolution=bestKnownSolution)
        
        if !isnothing(bound) && isBetterSolutionValue(model, bound, bestKnownValue)
            ddExact = false
        else
            dd = nothing
            bound = bestKnownValue
            relaxedSol = nothing
            ddExact = true
        end
    else
        dd = nothing
        bound = bestKnownValue
        relaxedSol = nothing
        ddExact = true
    end

    return (dd, bound, relaxedSol, ddExact), bestKnownSolution, bestKnownValue
end

"""
    addDDToQueue!(queue::Vector{Tuple{T, U, Vector{V}, Bool}}, block::Tuple{T, U, Vector{V}, Bool}, model::SequencingModel) where{T<:Vector{Vector{RelaxedSequencingNode}}, U<:Real, V<:Int}

Adds a tuple containing a relaxed dd and related useful information to a sorted processing queue (at the correct index)

# Arguments
- `queue::Vector{Tuple{T, U, Vector{V}, Bool}}`: The queue
- `block::Tuple{T, U, Vector{V}, Bool}`: The block to add to the queue
- `model::SequencingModel`: The sequencing problem

"""
function addDDToQueue!(queue::Vector{Tuple{T, U, Vector{V}, Bool}}, block::Tuple{T, U, Vector{V}, Bool}, model::SequencingModel) where{T<:Vector{Vector{RelaxedSequencingNode}}, U<:Real, V<:Int}
    #put worst bound at back of queue
    if isempty(queue)
        push!(queue, block)
    else
        insert!(
            queue,
            searchsorted(queue, block, by = x -> x[2],rev= (model.objective!=minimization)).start,
            block
        )
    end
end

"""
    finishPeelAndBound(bestKnownSolution::Vector{Int}, bestKnownValue::T;startTime::Union{U,Nothing}=nothing, timeElapsed::V=0,fileName::Union{Nothing,String}=nothing)where{T<:Real, U<:Real, V<:Real}

Print the optimal solution

# Arguments
- `bestKnownSolution::Vector{Int}`: The solution as an ordered list of nodes
- `bestKnownValue::T`: The cost of the solution
- `startTime::U`: the time the program started (in seconds since the epoch)
- `additionalTimeElapsed::V`: additional time used (in seconds)
- `fileName::Union{Nothing,String}`: An optional file to write the log statement to
"""
function finishPeelAndBound(bestKnownSolution::Vector{Int}, bestKnownValue::T;startTime::Union{U,Nothing}=nothing, timeElapsed::V=0,fileName::Union{Nothing,String}=nothing)where{T<:Real, U<:Real, V<:Real}
    logSequencingSolution(bestKnownSolution, bestKnownValue,firstLine="Optimal Solution Found:", fileName=fileName)
    if !isnothing(startTime)
        logRuntime(startTime, timeElapsed, fileName=fileName)
    end
    println()
end