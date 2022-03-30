"""
relaxedSequencingDD(
    model::SequencingModel,
    maxWidth::U,
    splitOrder::Vector{U};

    dd::Union{Vector{Vector{RelaxedSequencingNode}},Vector{Vector{RelaxedStateSequencingNode}},Nothing}=nothing,
    startingRoot::Union{RelaxedSequencingNode,RelaxedStateSequencingNode,Nothing}=nothing,
    bestKnownSolution::Union{Vector{U},Nothing}=nothing,
    bestKnownValue::Union{T,Nothing}=nothing,
    rrbSetting::Union{RRBSetting,Nothing}=nothing,
    rrbPacket=nothing,
    timeLimit::Union{Nothing, V}=nothing,
    timeElapsed::W=0.0,
    startTime::X=time(),
    widthOne::Bool=false,
    debugOn::Bool=false
)where{T<:Real,U<:Int}

Generate a relaxed bound for a sequencing problem using a relaxed decision diagram
Returns dd (Vector{Vector{RelaxedSequencingNode}}), bound (Real), relaxedSol (Vector{Int}), bestKnownSolution (Vector{Int}), bestKnownValue (Real), exact (Bool)
'exact' is true if the diagram exactly represents the problem (or subproblem) it was given

# Arguments
- `model::SequencingModel`: The sequencing problem
- `maxWidth::Int`: The max width allowed in a decision diagram layer
- `splitOrder::Vector{U}`: A heuristic ordering of nodes by importance
- `dd::Union{Vector{Vector{RelaxedSequencingNode}},Vector{Vector{RelaxedStateSequencingNode}},Nothing}`: An optional diagram to start from
- `startingRoot::Union{RelaxedSequencingNode,RelaxedStateSequencingNode,Nothing}`: An optional root node to start from
- `bestKnownSolution::Union{Vector{Int},Nothing}`: An optional parameter that gives the DD a starting solution
- `bestKnownValue::Union{T,Nothing}`: An optional parameter, The value of the best known feasible solution
- `rrbSetting::Union{RRBSetting,Nothing}`: An optional parameter that chooses an rough relaxed bounding method
- `rrbPacket`: Optional data to be passed to the rrbFunction
- `timeLimit::Union{Nothing, V}`: An optional parameter that adds a time limit
- `timeElapsed::W`: An optional parameter denoting how much time passed already
- `startTime::X`: An optional parameter denoting when the countdown started if there is a timer
- `widthOne::Bool`: An optional parameter to indicate that a widdth one diagram is the starting diagram
- `rootDomain::BitVector`: the domain of the rootNode if using widthOne
- `debugOn::Bool`: An optional parameter that can be used to turn on debug statements in other functions
"""
function relaxedSequencingDD(
    model::SequencingModel,
    maxWidth::U,
    splitOrder::Vector{U};
    dd::Union{Vector{Vector{RelaxedSequencingNode}},Vector{Vector{RelaxedStateSequencingNode}},Nothing}=nothing,
    startingRoot::Union{RelaxedSequencingNode,RelaxedStateSequencingNode,Nothing}=nothing,
    bestKnownSolution::Union{Vector{U},Nothing}=nothing,
    bestKnownValue::Union{T,Nothing}=nothing,
    rrbSetting::Union{RRBSetting,Nothing}=nothing,
    rrbPacket=nothing,
    timeLimit::Union{Nothing, V}=nothing,
    timeElapsed::W=0.0,
    startTime::X=time(),
    widthOne::Bool=false,
    rootDomain::BitVector=BitVector(),
    debugOn::Bool=false
)where{T<:Real,U<:Int, V<:Real, W<:Real, X<:Real}
    timedOut() = isTimedOut(timeLimit,startTime, timeElapsed)
    #create initial DD from initial constraints and domains
    if isnothing(dd)
        if isnothing(startingRoot)
            if widthOne
                dd = createWidthOneDDFromNode(model::SequencingModel)
            else
                dd = createSquareStartDD(model)
            end
        else
            if widthOne
                dd = createWidthOneDDFromNode(model,rootNode=startingRoot,masterDomain=rootDomain)
            else
                dd = createSquareStartDD(model, rootNode=startingRoot)
            end
        end
    end
    #for each layer
    for i in 1:length(dd)-1
        #only run this loop if the dd isnt too wide (or its the exact layer)
        if length(dd[i])<maxWidth || i==2
            #for each job in the ordering
            for assignment in splitOrder
                queue = selectNodes(dd[i],assignment)
                for node in queue
                    #split node into two new nodes
                    newNode = splitNode!(node, assignment, model)
            
                    #trim useless arcs
                    filterArcs!(node.arcsOut,model;bestKnownValue=bestKnownValue, rrbSetting=rrbSetting, rrbPacket=rrbPacket)
                    filterArcs!(newNode.arcsOut,model;bestKnownValue=bestKnownValue, rrbSetting=rrbSetting, rrbPacket=rrbPacket)
        
                    #if newNode has remaining out arcs, add it to the dd
                    if !isempty(newNode.arcsOut)
                        push!(dd[i], newNode)
                    else #delete in arcs so parent nodes may be removed in subsequent passes
                        deleteArcs!(newNode.arcsIn)
                    end
                    #if old node has no out arcs, remove it from the dd
                    if isempty(node.arcsOut)
                        deleteArcs!(node.arcsIn)
                        filter!(x->x!=node, dd[i])
                    end
                    #continue until max width of layer is reached except on the first layer
                    if length(dd[i])>=maxWidth && i!= 2
                        break
                    end
                    #timeout check
                    if timedOut()
                        return Vector{Vector{RelaxedSequencingNode}}(), 0, bestKnownSolution,bestKnownSolution,bestKnownValue, false
                    end
                end#queue loop
                #continue until max width of layer is reached except on the first layer
                if length(dd[i])>=maxWidth && i!= 2
                    break
                end
                #timeout check
                if timedOut()
                    return Vector{Vector{RelaxedSequencingNode}}(), 0, bestKnownSolution,bestKnownSolution,bestKnownValue, false
                end
            end#end assignment iteration
        end#if statement

        #timeout check
        if timedOut()
            return Vector{Vector{RelaxedSequencingNode}}(), 0, bestKnownSolution,bestKnownSolution,bestKnownValue, false
        end
        if i<length(dd)
            for node in dd[i+1]
                if isempty(node.arcsIn)
                    deleteArcs!(node.arcsOut)
                end
            end
            filter!(x->length(x.arcsIn)!=0, dd[i+1])
        end
        if dd isa Vector{Vector{RelaxedStateSequencingNode}}
            for node in dd[i]
                updateOutgoingArcValues!(model, node)
            end
        end
        #prep the next layer
        for node in dd[i+1]
            updateNodeDownVariables!(node,model)
            if node isa RelaxedStateSequencingNode
                updateOutgoingArcValues!(model, node)
            end
            filterArcs!(node.arcsOut,model;bestKnownValue=bestKnownValue, rrbSetting=rrbSetting, rrbPacket=rrbPacket)
            
            if isempty(node.arcsOut) && i+1 != length(dd)
                deleteArcs!(node.arcsIn)
            end
        end
        filter!(x->length(x.arcsIn)!=0, dd[i+1])
        #timeout check
        if timedOut()
            return Vector{Vector{RelaxedSequencingNode}}(), 0, bestKnownSolution,bestKnownSolution,bestKnownValue, false
        end
    end#end layer iteration

    #check if empty
    if isempty(last(dd))
        return Vector{Vector{RelaxedSequencingNode}}(), bestKnownValue, bestKnownSolution,bestKnownSolution,bestKnownValue, true
    end

    #timeout check
    if timedOut()
        return Vector{Vector{RelaxedSequencingNode}}(), 0, bestKnownSolution,bestKnownSolution,bestKnownValue, false
    end

    #iteratively delete nodes going up
    bottomUpUpdate(model, dd, bestKnownValue)

    if widthOne
        return dd, last(dd)[1].lengthToRoot, [], bestKnownSolution, bestKnownValue, last(dd)[1].exact
    else
        #update best value
        bestKnownSolution, bestKnownValue = checkFeasibleTerminalSolutions!(model, last(dd)[1],bestKnownValue=bestKnownValue,bestKnownSolution=bestKnownSolution)
        return dd, last(dd)[1].lengthToRoot, last(dd)[1].bestSeq, bestKnownSolution, bestKnownValue, last(dd)[1].exact
    end 
end

"""
    getFrontierNode(dd::Vector{Vector{RelaxedSequencingNode}}, relaxedSolution::Vector{T})where{T<:Int}

Get the deepest layer indexed exact node from the best path through the diagram

# Arguments
- `dd::Vector{Vector{RelaxedSequencingNode}}`: The dd to search
- `relaxedSolution::Vector{T}`: The problem solution
"""
function getFrontierNode(dd::Vector{Vector{RelaxedSequencingNode}}, relaxedSolution::Vector{T})where{T<:Int}
    currentNode = dd[1][1]
    for i in length(currentNode):length(relaxedSolution)
        for arc in currentNode.arcsOut
            if relaxedSolution[i+1] == arc.assignment
                if arc.destination.exact
                    currentNode = arc.destination
                    break
                else
                    return currentNode, (i - length(dd[1][1]) + 1)
                end
            end
        end
    end
end

"""
    getLastExactNode(dd::Vector{Vector{RelaxedSequencingNode}}, relaxedSolution::Vector{T})where{T<:Int}

Get the first node from the best path through the diagram that is a parent to a non-exact node

# Arguments
- `dd::Vector{Vector{RelaxedSequencingNode}}`: The dd to search
- `relaxedSolution::Vector{T}`: The problem solution
"""
function getLastExactNode(dd::Vector{Vector{RelaxedSequencingNode}}, relaxedSolution::Vector{T})where{T<:Int}
    currentNode = dd[1][1]
    for i in length(currentNode):length(relaxedSolution)
        nextNode = nothing
        for arc in currentNode.arcsOut
            if relaxedSolution[i+1] == arc.assignment
                if arc.destination.exact
                    nextNode = arc.destination
                else
                    return currentNode, (i - length(dd[1][1]) + 1)
                end
            elseif !arc.destination.exact
                return currentNode, (i - length(dd[1][1]) + 1)
            end
        end
        currentNode = nextNode
    end
end

"""
    getLastExactLayer(dd::Vector{Vector{RelaxedStateSequencingNode}})

Get the last exact layer of a dd

# Arguments
- `dd::Vector{Vector{RelaxedStateSequencingNode}}`: The dd to search
"""
function getLastExactLayer(dd::Vector{Vector{RelaxedStateSequencingNode}})
    currentLayer = dd[1]
    for layer in dd
        for node in layer
            if !node.exact
                return currentLayer
            end
        end
        currentLayer = layer
    end
end