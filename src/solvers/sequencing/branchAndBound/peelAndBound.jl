"""
peelAndBound(
    model::SequencingModel,
    maxWidth::U;

    orderingType::OrderingOption=edge,
    bestKnownSolution::Union{Array{U},Nothing}=nothing,
    bestKnownValue::Union{T,Nothing}=nothing,
    peelSetting::PeelSetting=peelf,
    rrbSetting::Union{RRBSetting,Nothing}=nothing,
    rrbPacket=nothing,
    timeLimit::Union{Nothing, W}=nothing,
    timeElapsed::X=0.0,
    loggingOn::Bool=true,
    debugOn::Bool=false,
    fileName::Union{Nothing,String}=nothing
)where{T<:Real,U<:Int, W<:Real, X<:Real}

Solve a sequencing problem to optimality using a peel and bound scheme
Returns bestKnownSolution (Vector{Int}), bestKnownValue (Real), bound (Real)

# Arguments
- `model::SequencingModel`: The sequencing problem
- `maxWidth::Int`: The max width allowed in a decision diagram layer
- `orderingType::OrderingOption`:Which heuristic to use to sort the edges, defaults to 'edge'
- `bestKnownSolution::Union{Vector{Int},Nothing}`: An optional parameter that gives the DD a starting solution
- `bestKnownValue::Union{T,Nothing}`: An optional parameter, The value of the best known feasible solution
- `peelSetting::PeelSetting`:Which heuristic to use to decide where to start the peel from, defaults to 'peelf'
- `rrbSetting::Union{RRBSetting,Nothing}`: An optional parameter that chooses an rough relaxed bounding method
- `rrbPacket`: Optional data to be passed to the rrbFunction
- `timeLimit::Union{Nothing, W}`: An optional parameter that adds a time limit
- `timeElapsed::X`: An optional parameter denoting how much time passed already
- `loggingOn::Bool`: An optional parameter that can be used to turn off logging to the console and local file
- `debugOn::Bool`: An optional parameter that can be used to turn on debug statements in other functions
- `fileName::Union{Nothing,String}`: The file to write results to
"""
function peelAndBound(
    model::SequencingModel,
    maxWidth::U;
    orderingType::OrderingOption=edge,
    bestKnownSolution::Union{Array{U},Nothing}=nothing,
    bestKnownValue::Union{T,Nothing}=nothing,
    peelSetting::PeelSetting=peelf,
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
    cleanRelaxedDD(dd) = relaxedSequencingDD(model, maxWidth,splitOrder, dd = dd, debugOn=debugOn,bestKnownSolution=bestKnownSolution, bestKnownValue=bestKnownValue, rrbSetting=rrbSetting, rrbPacket=rrbPacket, timeLimit=timeLimit, timeElapsed=timeElapsed, startTime=startTime)
    splitOrder = findSplitOrdering(model,orderingType)
    cleanPeel!(relaxedDD, relaxedSolution, frontierNode, fNodeIndex) = peelDD!(model, maxWidth, relaxedDD, relaxedSolution, frontierNode, fNodeIndex, splitOrder,bestKnownValue,rrbSetting,rrbPacket,loggingOn,debugOn)
    cleanPeelAndRemove!(dd, frontierNode, fNodeIndex, bestKnownValue) = peelAndRemove!(model,dd,frontierNode,fNodeIndex,bestKnownValue,rrbSetting,rrbPacket,loggingOn,debugOn)
    
    #get initial bounds
    if isnothing(bestKnownSolution)
        bestKnownSolution, bestKnownValue, exact = restrictedSequencingDD(model,maxWidth,debugOn=debugOn)
        if exact
            if loggingOn
                finishPeelAndBound(bestKnownSolution, bestKnownValue,startTime=startTime, timeElapsed=timeElapsed, fileName=fileName)
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
                finishPeelAndBound(bestKnownSolution, bestKnownValue,startTime=startTime, timeElapsed=timeElapsed, fileName=fileName)
            end
            return bestKnownSolution, bestKnownValue, bestKnownValue
        end
        if timedOut()
            return bestKnownSolution, bestKnownValue, 0
        end
    end
    if loggingOn
        logSequencingSolution(bestKnownSolution, bestKnownValue,firstLine="Initial Solution:", fileName=fileName)
        logRuntime(startTime, timeElapsed, fileName=fileName)
    end
    #generate initial relaxed DD
    relaxedDD, relaxedBound, relaxedSolution, bestKnownSolution, bestKnownValue, exact = cleanRelaxedDD(nothing)

    if exact || relaxedBound == bestKnownValue
        if loggingOn
            finishPeelAndBound(bestKnownSolution, bestKnownValue,startTime=startTime, timeElapsed=timeElapsed, fileName=fileName)
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
    ddType = Vector{Vector{RelaxedSequencingNode}}
    valueType = typeof(bestKnownValue)
    queue = Vector{Tuple{ddType, valueType, Vector{Int}, Bool}}()
    push!(queue, (relaxedDD, bestBound, relaxedSolution, false))

    if loggingOn
        logQueueLength(length(queue),fileName=fileName)
    end
    #process queue
    while !isempty(queue)
        currentValue = bestKnownValue
        currentBlock = pop!(queue)
        if debugOn
            if isempty(queue)
                println("In Queue: ",length(queue)+1, "  From: ",currentBlock[2], " - ",currentBlock[2], "\n")
            else
                println("In Queue: ",length(queue)+1, "  From: ",currentBlock[2], " - ",queue[1][2], "\n")
            end
        end

        if peelSetting == peelf
            frontierNode, fNodeIndex = getFrontierNode(currentBlock[1], currentBlock[3])
        else
            frontierNode, fNodeIndex = getLastExactNode(currentBlock[1], currentBlock[3])
        end
        filterArcs!(frontierNode.arcsOut,model;bestKnownValue=bestKnownValue, rrbSetting=rrbSetting, rrbPacket=rrbPacket)

        bestKnownSolution, bestKnownValue, restrictedExact = cleanRestrictedDD(frontierNode)
        #time limit block
        if timedOut()
            return bestKnownSolution, bestKnownValue, bestBound
        end

        if restrictedExact
            #delete the node
            ddInfo, bestPeeledSolution, bestPeeledValue = cleanPeelAndRemove!(currentBlock[1], frontierNode, fNodeIndex, bestKnownValue)

            if !ddInfo[4] && isBetterSolutionValue(model,ddInfo[2], bestKnownValue)
                addDDToQueue!(queue, ddInfo, model)
            elseif isBetterSolutionValue(model, bestPeeledValue, bestKnownValue)
                bestKnownValue = bestPeeledValue
                bestKnownSolution = bestPeeledSolution
            end
        else
            #peel the dd
            ddInfo, newDDInfo, bestPeeledSolution, bestPeeledValue = cleanPeel!(currentBlock[1], currentBlock[3], frontierNode, fNodeIndex)

            if isBetterSolutionValue(model, bestPeeledValue, bestKnownValue)
                bestKnownValue = bestPeeledValue
                bestKnownSolution = bestPeeledSolution
            end

            #add the DDs to the queue
            if !ddInfo[4] && isBetterSolutionValue(model,ddInfo[2], bestKnownValue)
                addDDToQueue!(queue, ddInfo, model)
            end
            if !newDDInfo[4] && isBetterSolutionValue(model,newDDInfo[2], bestKnownValue)
                addDDToQueue!(queue, newDDInfo, model)
            end
        end

        while !isempty(queue) && isBetterSolutionValue(model, bestKnownValue, queue[1][2])
            deleteat!(queue, 1)
        end
        if !isempty(queue)
            if loggingOn && isBetterSolutionValue(model, bestKnownValue, currentValue)
                logSequencingSolution(bestKnownSolution, bestKnownValue,firstLine="Improved Solution:", fileName=fileName)
                logOptimalityGap(last(queue)[2], bestKnownValue, fileName=fileName)
                logRuntime(startTime, timeElapsed, fileName=fileName)
                logQueueLength(length(queue),fileName=fileName)
            elseif loggingOn && isBetterSolutionValue(model, bestBound, last(queue)[2])
                logOptimalityGap(last(queue)[2], bestKnownValue, fileName=fileName)
                logRuntime(startTime, timeElapsed, fileName=fileName)
                logQueueLength(length(queue),fileName=fileName)
            end
            bestBound = last(queue)[2]
        end

        #time limit block
        if timedOut()
            return bestKnownSolution, bestKnownValue, bestBound
        end
    end

    if loggingOn
        finishPeelAndBound(bestKnownSolution, bestKnownValue, startTime=startTime, timeElapsed=timeElapsed, fileName=fileName)
    end
    return bestKnownSolution, bestKnownValue, bestKnownValue
end