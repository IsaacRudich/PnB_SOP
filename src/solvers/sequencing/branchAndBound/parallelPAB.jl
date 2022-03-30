#export JULIA_NUM_THREADS=4 (before starting julia)

"""
parallelPeelAndBound(
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
- `debugOn::Bool`: An optional parameter that can be used to turn on debug statements in other functions,
- `fileName::Union{Nothing,String}`: The file to write results to
"""
function parallelPeelAndBound(
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

    #get initial bounds
    if isnothing(bestKnownSolution)
        bestKnownSolution, bestKnownValue, exact = restrictedSequencingDD(model,2048,debugOn=debugOn)
        if exact
            if loggingOn
                finishPeelAndBound(bestKnownSolution, bestKnownValue,startTime=startTime, timeElapsed=timeElapsed)
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
                finishPeelAndBound(bestKnownSolution, bestKnownValue,startTime=startTime, timeElapsed=timeElapsed)
            end
            return bestKnownSolution, bestKnownValue, bestKnownValue
        end
        if timedOut()
            return bestKnownSolution, bestKnownValue, 0
        end
    end
    if loggingOn
        logSequencingSolution(bestKnownSolution, bestKnownValue,firstLine="Initial Solution:")
        logRuntime(startTime, timeElapsed)
        println()
    end
    #generate initial relaxed DD
    relaxedDD, relaxedBound, relaxedSolution, bestKnownSolution, bestKnownValue, exact = cleanRelaxedDD(nothing)
    if exact || relaxedBound == bestKnownValue
        if loggingOn
            finishPeelAndBound(bestKnownSolution, bestKnownValue,startTime=startTime, timeElapsed=timeElapsed)
        end
        return bestKnownSolution, bestKnownValue, bestKnownValue
    else
        bestBound = relaxedBound
        if loggingOn
            logOptimalityGap(bestBound, bestKnownValue)
            logRuntime(startTime, timeElapsed)
            println()
        end
    end

    #time limit block
    if timedOut()
        return bestKnownSolution, bestKnownValue, bestBound
    end

    #initialize queue
    ddType = Vector{Vector{RelaxedSequencingNode}}
    queue = Vector{Tuple{ddType, typeof(bestKnownValue), Vector{Int}, Bool}}()
    push!(queue, (relaxedDD, bestBound, relaxedSolution, false))

    queueLock = ReentrantLock()
    countLock = ReentrantLock()
    solutionLock = ReentrantLock()
    runningJobs = Ref(0)
    bestKnownValue = Ref(bestKnownValue)
    #process queue 
    while !isempty(queue) || runningJobs[] > 0 
        while runningJobs[] >= Threads.nthreads() 
            sleep(.001)
        end
        currentBlock = nothing
        begin
            lock(queueLock)
            lock(countLock)
            try
                if !isempty(queue)
                    currentBlock = pop!(queue)
                    runningJobs[] += 1
                    if loggingOn
                        println("Running Jobs: ", runningJobs[], " | Queue: ", length(queue), " | Local Bound: ", currentBlock[2])
                        logRuntime(startTime, timeElapsed)
                    end
                end
            finally
                unlock(queueLock)
                unlock(countLock)
            end
        end
        if !isnothing(currentBlock)
            Threads.@spawn asyncProcessNode!(model, maxWidth, queue, queueLock, solutionLock, countLock, runningJobs,currentBlock, copy(bestKnownSolution), bestKnownValue[], bestKnownSolution, bestKnownValue, splitOrder, peelSetting, rrbSetting, rrbPacket, timeLimit, timeElapsed, startTime, debugOn, loggingOn)
        end #process node
    end

    #log solution
    if loggingOn
        finishPeelAndBound(bestKnownSolution, bestKnownValue[], startTime=startTime, timeElapsed=timeElapsed)
    end
    return bestKnownSolution, bestKnownValue[], bestKnownValue[]
end















"""
    asyncProcessNode!(model::SequencingModel, maxWidth::Int, queue::Vector{Tuple{Vector{Vector{RelaxedSequencingNode}}, U, Vector{Int}, Bool}}, queueLock::ReentrantLock, solutionLock::ReentrantLock,countLock::ReentrantLock, jobCount::Ref{Int}, currentBlock::Tuple{Vector{Vector{RelaxedSequencingNode}}, U, Vector{Int}, Bool}, bestKnownSolution::Union{Vector{Int},Nothing}, bestKnownValue::Union{T,Nothing}, solutionRef::Union{Vector{Int},Nothing}, valueRef::Ref{T}, splitOrder::Vector{Int}, peelSetting::PeelSetting, rrbSetting::Union{RRBSetting,Nothing}, rrbPacket, timeLimit::Union{Nothing, W}, timeElapsed::X, startTime::Y, debugOn::Bool, loggingOn::Bool)where{T<:Real, U<:Real, W<:Real, X<:Real, Y<:Real}

An asynchronous call to the process the next node
Returns bestKnownSolution (Vector{Int}), bestKnownValue (Real), jobs added (Int)

# Arguments
- `model::SequencingModel`: The sequencing problem
- `maxWidth::Int`: The max width allowed in a decision diagram layer
- `queue`::Vector{Tuple{Vector{Vector{RelaxedSequencingNode}}, U, Vector{Int}, Bool}}: The node queue
- `queueLock::ReentrantLock`: The queue lock for parallel processing
- `solutionLock::ReentrantLock`:
- `countLock::::ReentrantLock`:
- `jobCount::Ref{Int}`:
- `currentBlock::Tuple{Vector{Vector{RelaxedSequencingNode}}, U, Vector{Int}, Bool}`: the dd to be processed
- `bestKnownSolution::Union{Vector{Int},Nothing}`: An optional parameter that gives the DD a starting solution
- `bestKnownValue::Union{T,Nothing}`: An optional parameter, The value of the best known feasible solution
- `solutionRef::Union{Vector{Int},Nothing}`:
- `valueRef::Ref{T}`:
- `splitOrder::Vector{Int}`: A heuristic ordering of nodes by importance
- `peelSetting::PeelSetting`:Which heuristic to use to decide where to start the peel from, defaults to 'peelf'
- `rrbSetting::Union{RRBSetting,Nothing}`: An optional parameter that chooses an rough relaxed bounding method
- `rrbPacket`: Optional data to be passed to the rrbFunction
- `timeLimit::Union{Nothing, W}`: An optional parameter that adds a time limit
- `timeElapsed::X`: An optional parameter denoting how much time passed already
- `startTime::Y`: An optional parameter denoting when the countdown started if there is a timer
- `loggingOn::Bool`: An optional parameter that can be used to turn off logging to the console and local file
- `debugOn::Bool`: An optional parameter that can be used to turn on debug statements in other functions
"""
function asyncProcessNode!(model::SequencingModel, maxWidth::Int, queue::Vector{Tuple{Vector{Vector{RelaxedSequencingNode}}, U, Vector{Int}, Bool}}, queueLock::ReentrantLock, solutionLock::ReentrantLock,countLock::ReentrantLock, jobCount::Ref{Int}, currentBlock::Tuple{Vector{Vector{RelaxedSequencingNode}}, U, Vector{Int}, Bool}, bestKnownSolution::Union{Vector{Int},Nothing}, bestKnownValue::Union{T,Nothing}, solutionRef::Union{Vector{Int},Nothing}, valueRef::Ref{T}, splitOrder::Vector{Int}, peelSetting::PeelSetting, rrbSetting::Union{RRBSetting,Nothing}, rrbPacket, timeLimit::Union{Nothing, W}, timeElapsed::X, startTime::Y, debugOn::Bool, loggingOn::Bool)where{T<:Real, U<:Real, W<:Real, X<:Real, Y<:Real}
    if peelSetting == peelf
        frontierNode, fNodeIndex = getFrontierNode(currentBlock[1], currentBlock[3])
    else
        frontierNode, fNodeIndex = getLastExactNode(currentBlock[1], currentBlock[3])
    end
    filterArcs!(frontierNode.arcsOut,model;bestKnownValue=bestKnownValue, rrbSetting=rrbSetting, rrbPacket=rrbPacket)
    bestKnownSolution, bestKnownValue, restrictedExact = restrictedSequencingDD(model, maxWidth,rootNode = frontierNode, debugOn=debugOn,bestKnownSolution=bestKnownSolution, bestKnownValue=bestKnownValue, rrbSetting=rrbSetting, rrbPacket=rrbPacket)

    #time limit block
    if isTimedOut(timeLimit,startTime, timeElapsed)
        return bestKnownSolution, bestKnownValue
    end

    #if node is closed
    if restrictedExact
        #delete the node
        ddInfo, bestPeeledSolution, bestPeeledValue = peelAndRemove!(model,currentBlock[1],frontierNode,fNodeIndex,bestKnownValue,rrbSetting,rrbPacket,loggingOn,debugOn)
        
        #if local bound
        if !ddInfo[4] && isBetterSolutionValue(model,ddInfo[2], bestKnownValue)
            begin
                lock(queueLock)
                try
                    addDDToQueue!(queue, ddInfo, model)
                finally
                    unlock(queueLock)
                end
            end
        elseif isBetterSolutionValue(model, bestPeeledValue, bestKnownValue)
            bestKnownValue = bestPeeledValue
            bestKnownSolution = bestPeeledSolution
        end
    else #if node is not closed 
        #peel the dd
        ddInfo, newDDInfo, bestPeeledSolution, bestPeeledValue = peelDD!(model, maxWidth, currentBlock[1], currentBlock[3], frontierNode, fNodeIndex, splitOrder,bestKnownValue,rrbSetting,rrbPacket,loggingOn,debugOn)
        if isBetterSolutionValue(model, bestPeeledValue, bestKnownValue)
            bestKnownValue = bestPeeledValue
            bestKnownSolution = bestPeeledSolution
        end
        #add the DDs to the queue
        #if both nodes need to be added
        if !ddInfo[4] && isBetterSolutionValue(model,ddInfo[2], bestKnownValue) && !newDDInfo[4] && isBetterSolutionValue(model,newDDInfo[2], bestKnownValue)
            begin
                lock(queueLock)
                try
                    addDDToQueue!(queue, ddInfo, model)
                    addDDToQueue!(queue, newDDInfo, model)
                finally
                    unlock(queueLock)
                end
            end
        elseif !ddInfo[4] && isBetterSolutionValue(model,ddInfo[2], bestKnownValue)
            begin
                lock(queueLock)
                try
                    addDDToQueue!(queue, ddInfo, model)
                finally
                    unlock(queueLock)
                end
            end
        elseif !newDDInfo[4] && isBetterSolutionValue(model,newDDInfo[2], bestKnownValue)
            begin
                lock(queueLock)
                try
                    addDDToQueue!(queue, newDDInfo, model)
                finally
                    unlock(queueLock)
                end
            end
        end
    end
    #update solution value
    if isBetterSolutionValue(model, bestKnownValue, valueRef[])
        begin
            lock(solutionLock)
            try
                valueRef[] = bestKnownValue
                solutionRef = bestKnownSolution
                if loggingOn
                    logSequencingSolution(bestKnownSolution, bestKnownValue,firstLine="Improved Solution:")
                    logRuntime(startTime, timeElapsed)
                    println()
                end
            finally
                unlock(solutionLock)
            end
        end
    end
    #update job count
    begin
        lock(countLock)
        try
            jobCount[] -=1
        finally
            unlock(countLock)
        end
    end
    return bestKnownSolution, bestKnownValue
    #=if !isempty(queue)
        if loggingOn && isBetterSolutionValue(model, bestKnownValue, currentValue)
            logSequencingSolution(bestKnownSolution, bestKnownValue,firstLine="Improved Solution:")
            logOptimalityGap(last(queue)[2], bestKnownValue)
            logRuntime(startTime, timeElapsed)
        elseif loggingOn && isBetterSolutionValue(model, bestBound, last(queue)[2])
            println("Running Jobs: $running_jobs")
            logOptimalityGap(last(queue)[2], bestKnownValue)
            logRuntime(startTime, timeElapsed)
        end
        begin
            lock(boundLock)
            try
                bestBound = last(queue)[2]
            finally
                unlock(boundLock)
            end
        end
    end

    #time limit block
    isTimedOut(timeLimit,startTime, timeElapsed)
        return bestKnownSolution, bestKnownValue, bestBound
    end=#
end