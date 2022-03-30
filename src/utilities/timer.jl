"""
    isTimedOut(timeLimit::Union{T,Nothing},startTime::U, additionalTimeElapsed::V)where{T<:Real, U<:Real, V<:Real}

Check if a timer has been exceeded, always true if timeLimit is nothing

# Arguments
- `timeLimit::Union{T,Nothing}`: the time limit in seconds (or nothing)
- `startTime::U`: the time the program started (in seconds since the epoch)
- `additionalTimeElapsed::V`: additional time used (in seconds)
"""
function isTimedOut(timeLimit::Union{T,Nothing},startTime::U, additionalTimeElapsed::V)where{T<:Real, U<:Real, V<:Real}
    if isnothing(timeLimit)
        return false
    end
    if (time() - startTime + additionalTimeElapsed) >= timeLimit
        return true
    else
        return false
    end
end