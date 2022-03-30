"""
branchAndBound(
    model::SequencingModel,
    maxWidth::U;

    orderingType::OrderingOption=edge,
    bestKnownSolution::Union{Array{U},Nothing}=nothing,
    bestKnownValue::Union{T,Nothing}=nothing,
    rrbSetting::Union{RRBSetting,Nothing}=nothing,
    rrbPacket=nothing,
    timeLimit::Union{Nothing, W}=nothing,
    timeElapsed::X=0.0,
    loggingOn::Bool=true,
    debugOn::Bool=false,
    fileName::Union{Nothing,String}=nothing
)where{T<:Real,U<:Int, W<:Real, X<:Real}

Solve a sequencing problem to optimality using a branch and bound scheme
Returns bestKnownSolution (Vector{Int}), bestKnownValue (Real), bound (Real)

# Arguments
- `model::SequencingModel`: The sequencing problem
- `maxWidth::Int`: The max width allowed in a decision diagram layer
- `orderingType::OrderingOption`:Which heuristic to use to sort the edges, defaults to 'edge'
- `bestKnownSolution::Union{Vector{Int},Nothing}`: An optional parameter that gives the DD a starting solution
- `bestKnownValue::Union{T,Nothing}`: An optional parameter, The value of the best known feasible solution
- `rrbSetting::Union{RRBSetting,Nothing}`: An optional parameter that chooses an rough relaxed bounding method
- `rrbPacket`: Optional data to be passed to the rrbFunction
- `timeLimit::Union{Nothing, W}`: An optional parameter that adds a time limit
- `timeElapsed::X`: An optional parameter denoting how much time passed already
- `loggingOn::Bool`: An optional parameter that can be used to turn off logging to the console and local file
- `debugOn::Bool`: An optional parameter that can be used to turn on debug statements in other functions
- `fileName::Union{Nothing,String}`: The file to write results to
"""
function branchAndBound(
    model::SequencingModel,
    maxWidth::U;
    orderingType::OrderingOption=edge,
    bestKnownSolution::Union{Array{U},Nothing}=nothing,
    bestKnownValue::Union{T,Nothing}=nothing,
    rrbSetting::Union{RRBSetting,Nothing}=nothing,
    rrbPacket=nothing,
    timeLimit::Union{Nothing, W}=nothing,
    timeElapsed::X=0.0,
    loggingOn::Bool=true,
    debugOn::Bool=false,
    fileName::Union{Nothing,String}=nothing
)where{T<:Real,U<:Int, W<:Real, X<:Real}
    startTime = time()
    bestBound = nothing
    timedOut() = isTimedOut(timeLimit,startTime, timeElapsed)
    cleanRestrictedDD(node) = restrictedSequencingDD(model, maxWidth,rootNode = node, debugOn=debugOn,bestKnownSolution=bestKnownSolution, bestKnownValue=bestKnownValue, rrbSetting=rrbSetting, rrbPacket=rrbPacket)
    cleanRestrictedDD(node, domain) = restrictedSequencingDD(model, maxWidth,rootNode = node, debugOn=debugOn,bestKnownSolution=bestKnownSolution, bestKnownValue=bestKnownValue, rrbSetting=rrbSetting, rrbPacket=rrbPacket, rootDomain=domain)
    cleanRelaxedDD(dd) = relaxedSequencingDD(model, maxWidth,splitOrder, dd = dd, debugOn=debugOn,bestKnownSolution=bestKnownSolution, bestKnownValue=bestKnownValue, rrbSetting=rrbSetting, rrbPacket=rrbPacket, timeLimit=timeLimit, timeElapsed=timeElapsed, startTime=startTime, widthOne=true)
    splitOrder = findSplitOrdering(model,orderingType)

     #get initial bounds
     if isnothing(bestKnownSolution)
        bestKnownSolution, bestKnownValue, exact = restrictedSequencingDD(model,maxWidth,debugOn=debugOn)
        if exact
            if loggingOn
                finishBranchAndBound(bestKnownSolution, bestKnownValue,startTime=startTime, timeElapsed=timeElapsed, fileName=fileName)
            end
            return bestKnownSolution, bestKnownValue, bestKnownValue
        end
    end
    if timedOut()
        return bestKnownSolution, bestKnownValue, 0
    end
    #improve initial bounds using rrb
    if !isnothing(rrbSetting)
        bestKnownSolution, bestKnownValue, exact = cleanRestrictedDD(nothing)
        #if solution is found finish
        if exact
            if loggingOn
                finishBranchAndBound(bestKnownSolution, bestKnownValue,startTime=startTime, timeElapsed=timeElapsed, fileName=fileName)
            end
            return bestKnownSolution, bestKnownValue, bestKnownValue
        end
        if timedOut()
            return bestKnownSolution, bestKnownValue, 0
        end
    end
    if loggingOn
        logSequencingSolution(bestKnownSolution, bestKnownValue,firstLine="Initial Solution:",fileName=fileName)
        logRuntime(startTime, timeElapsed, fileName=fileName)
    end
    #generate initial relaxed DD
    relaxedDD, relaxedBound, relaxedSolution, bestKnownSolution, bestKnownValue, exact = cleanRelaxedDD(nothing)

    if exact || relaxedBound == bestKnownValue
        if loggingOn
            finishBranchAndBound(bestKnownSolution, bestKnownValue,startTime=startTime, timeElapsed=timeElapsed, fileName=fileName)
        end
        return bestKnownSolution, bestKnownValue, bestKnownValue
    elseif timedOut()
        return bestKnownSolution, bestKnownValue, 0
    else
        bestBound = relaxedBound
        if loggingOn
            logOptimalityGap(bestBound, bestKnownValue, fileName=fileName)
            logRuntime(startTime, timeElapsed, fileName=fileName)
        end
    end

    #initialize queue
    queue = Vector{Tuple{RelaxedStateSequencingNode, BitVector}}()
    for node in getLastExactLayer(relaxedDD)
        filterArcs!(node.arcsOut,model,bestKnownValue = bestKnownValue,rrbSetting=rrbSetting, rrbPacket=rrbPacket)
        if !isempty(node.arcsOut)
            addNodeToBnBQueue!(queue, node,model)
        end
    end
    relaxedDD = nothing
    if loggingOn
        logQueueLength(length(queue),fileName=fileName)
    end
    #process queue
    while !isempty(queue)
        currentValue = bestKnownValue
        currentBlock = pop!(queue)
        if debugOn
            if isempty(queue)
                println("In Queue: ",length(queue)+1, "  From: ",currentBlock[1].lengthToTerminal+currentBlock[1].lengthToRoot, " - ",currentBlock[1].lengthToTerminal+currentBlock[1].lengthToRoot, "\n")
            else
                println("In Queue: ",length(queue)+1, "  From: ",currentBlock[1].lengthToTerminal+currentBlock[1].lengthToRoot, " - ",queue[1][1].lengthToTerminal+queue[1][1].lengthToRoot, "\n")
            end
        end

        bestKnownSolution, bestKnownValue, restrictedExact = cleanRestrictedDD(currentBlock[1], currentBlock[2])
        #time limit block
        if timedOut()
            return bestKnownSolution, bestKnownValue, bestBound
        end

        if !restrictedExact
            #relax the node
            relaxedDD, localBound, emptyPath, bestKnownSolution, bestKnownValue, localExactness = cleanRelaxedDD(createWidthOneDDFromNode(model,rootNode=currentBlock[1],masterDomain=currentBlock[2]))
            #time limit block
            if timedOut()
                return bestKnownSolution, bestKnownValue, bestBound
            end
            if !localExactness
                #add the DDs to the queue
                for node in getLastExactLayer(relaxedDD)
                    if isBetterSolutionValue(model, node.lengthToRoot+node.lengthToTerminal, bestKnownValue)
                        filterArcs!(node.arcsOut,model,bestKnownValue = bestKnownValue,rrbSetting=rrbSetting, rrbPacket=rrbPacket)
                        if !isempty(node.arcsOut)
                            addNodeToBnBQueue!(queue, node,model)
                        end
                    end
                end
            end
        end

        while !isempty(queue) && isBetterSolutionValue(model, bestKnownValue, queue[1][1].lengthToRoot + queue[1][1].lengthToTerminal)
            deleteat!(queue, 1)
        end

        if !isempty(queue)
            currentBound = last(queue)[1].lengthToRoot + last(queue)[1].lengthToTerminal
            if loggingOn && isBetterSolutionValue(model, bestKnownValue, currentValue)
                logSequencingSolution(bestKnownSolution, bestKnownValue,firstLine="Improved Solution:", fileName=fileName)
                logOptimalityGap(currentBound, bestKnownValue, fileName=fileName)
                logRuntime(startTime, timeElapsed, fileName=fileName)
                logQueueLength(length(queue),fileName=fileName)
            elseif loggingOn && isBetterSolutionValue(model, bestBound, currentBound)
                logOptimalityGap(currentBound, bestKnownValue, fileName=fileName)
                logRuntime(startTime, timeElapsed, fileName=fileName)
                logQueueLength(length(queue),fileName=fileName)
            end
            bestBound = currentBound
        end

        #time limit block
        if timedOut()
            return bestKnownSolution, bestKnownValue, bestBound
        end
    end

    if loggingOn
        finishBranchAndBound(bestKnownSolution, bestKnownValue, startTime=startTime, timeElapsed=timeElapsed, fileName=fileName)
    end
    return bestKnownSolution, bestKnownValue, bestKnownValue
end


"""
    addNodeToBnBQueue!(queue::Vector{Tuple{RelaxedStateSequencingNode, BitVector}}, node::RelaxedStateSequencingNode,model::SequencingModel)

Adds a tuple containing a realxed node and its domain

# Arguments
- `queue::Vector{Tuple{RelaxedStateSequencingNode, BitVector}}`: The queue
- `node::RelaxedStateSequencingNode`: The block to add to the queue
- `model::SequencingModel`: The sequencing problem

"""
function addNodeToBnBQueue!(queue::Vector{Tuple{RelaxedStateSequencingNode, BitVector}}, node::RelaxedStateSequencingNode,model::SequencingModel)
    #get domain then delete ouut arcs
    domain = zeros(length(node.state))
    for arc in node.arcsOut
        domain[arc.assignment] = true
    end
    deleteArcs!(node.arcsOut)
    #delete in arcs
    deleteArcs!(node.arcsIn)

    block = (node, domain)
    #put worst bound at back of queue
    if isempty(queue)
        push!(queue, block)
    else
        insert!(
            queue,
            searchsorted(queue, block, by = x -> x[1].lengthToTerminal+x[1].lengthToRoot,rev= (model.objective!=minimization)).start,
            block
        )
    end
end

"""
    finishBranchAndBound(bestKnownSolution::Vector{Int}, bestKnownValue::T;startTime::Union{U,Nothing}=nothing, timeElapsed::V=0,fileName::Union{Nothing,String}=nothing)where{T<:Real, U<:Real, V<:Real}

Print the optimal solution

# Arguments
- `bestKnownSolution::Vector{Int}`: The solution as an ordered list of nodes
- `bestKnownValue::T`: The cost of the solution
- `startTime::U`: the time the program started (in seconds since the epoch)
- `additionalTimeElapsed::V`: additional time used (in seconds)
- `fileName::Union{Nothing,String}`: An optional file to write the log statement to
"""
function finishBranchAndBound(bestKnownSolution::Vector{Int}, bestKnownValue::T;startTime::Union{U,Nothing}=nothing, timeElapsed::V=0,fileName::Union{Nothing,String}=nothing)where{T<:Real, U<:Real, V<:Real}
    logSequencingSolution(bestKnownSolution, bestKnownValue,firstLine="Optimal Solution Found:",fileName=fileName)
    if !isnothing(startTime)
        logRuntime(startTime, timeElapsed,fileName=fileName)
    end
    println()
end