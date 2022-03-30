"""
    logSequencingSolution(bestKnownSolution::Vector{Int}, bestKnownValue::T;firstLine::Union{String, Nothing}=nothing,fileName::Union{Nothing,String}=nothing) where T<:Real

Print a solution to the terminal

# Arguments
- `bestKnownSolution::Vector{Int}`: The solution as an ordered list of nodes
- `bestKnownValue::T`: The cost of the solution
- `firstLine::Union{String, Nothing}`:An optional preface to the solution
- `fileName::Union{Nothing,String}`: An optional file to write the log statement to
"""
function logSequencingSolution(bestKnownSolution::Vector{Int}, bestKnownValue::T;firstLine::Union{String, Nothing}=nothing,fileName::Union{Nothing,String}=nothing) where T<:Real
    if !isnothing(firstLine)
        println(firstLine)
    end
    print("Sequence: ",)
    for i in 1:length(bestKnownSolution)
        print(bestKnownSolution[i], " ")
    end
    if isinteger(bestKnownValue)
        println("\nCost: ", Int(bestKnownValue))
    else
        println("\nCost: ", bestKnownValue)
    end
    if !isnothing(fileName)
        f = open(fileName,"a")
            if !isnothing(firstLine)
                write(f, string(firstLine, "\n"))
            end
            write(f, "Sequence: ")
            for i in 1:length(bestKnownSolution)
                write(f, string(bestKnownSolution[i], " "))
            end
            if isinteger(bestKnownValue)
                write(f, string("\nCost: ", Int(bestKnownValue),"\n"))
            else
                write(f, string("\nCost: ", bestKnownValue,"\n"))
            end
        close(f)
    end
end

"""
    logOptimalityGap(lower::T, upper::U;fileName::Union{Nothing,String}=nothing)where{T<:Real,U<:Real}

Print an optimality gap to the terminal

# Arguments
- `lower::T`: The smaller number
- `upper::U`: The larger number
- `fileName::Union{Nothing,String}`: An optional file to write the log statement to
"""
function logOptimalityGap(lower::T, upper::U;fileName::Union{Nothing,String}=nothing)where{T<:Real,U<:Real}
    println("Bounds: [",lower,",",upper,"]")
    println("Optimality Gap: ",round((upper-lower)/upper,digits=4)*100,"%")
    if !isnothing(fileName)
        f = open(fileName,"a")
            write(f, string("Bounds: [",lower,",",upper,"]", "\n"))
            write(f, string("Optimality Gap: ",round((upper-lower)/upper,digits=4)*100,"%", "\n"))
            write(f, string("Data: ",lower,",",upper,",",round((upper-lower)/upper,digits=4)*100,"\n"))
        close(f)
    end
end


"""
    logDDWidths(dd::Vector{Vector{RelaxedSequencingNode}})

Print the width of each layer of a dd to the terminal

# Arguments
- `dd::Vector{Vector{RelaxedSequencingNode}}`: The dd
"""
function logDDWidths(dd::Vector{Vector{RelaxedSequencingNode}})
    print("DD Layer Widths: [")
    if !isempty(dd)
        print(length(dd[1]))
        for i in 2:length(dd)
            print(",",length(dd[i]))
        end
    end
    println("]")
end


"""
    logRuntime(startTime::U, additionalTimeElapsed::V;fileName::Union{Nothing,String}=nothing)where{U<:Real, V<:Real}

logs the current runtime

# Arguments
- `startTime::U`: the time the program started (in seconds since the epoch)
- `additionalTimeElapsed::V`: additional time used (in seconds)
- `fileName::Union{Nothing,String}`: An optional file to write the log statement to
"""
function logRuntime(startTime::U, additionalTimeElapsed::V;fileName::Union{Nothing,String}=nothing)where{U<:Real, V<:Real}
    println(round(time() - startTime + additionalTimeElapsed, digits=2), " seconds")
    if !isnothing(fileName)
        f = open(fileName,"a")
            write(f, string(round(time() - startTime + additionalTimeElapsed, digits=2), " seconds", "\n"))
        close(f)
    end
end

"""
    logQueueLength(qLength::U;fileName::Union{Nothing,String}=nothing)where{U<:Int} 

Log the length of the queue

# Arguments
- `qLength::U`: the queue length
- `fileName::Union{Nothing,String}`: An optional file to write the log statement to
"""
function logQueueLength(qLength::U;fileName::Union{Nothing,String}=nothing)where{U<:Int}
    println("Nodes In Queue: ", qLength)
    if !isnothing(fileName)
        f = open(fileName,"a")
            write(f, string("Nodes In Queue: ",qLength, "\n"))
        close(f)
    end
end