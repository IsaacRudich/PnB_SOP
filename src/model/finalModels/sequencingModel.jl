struct SequencingModel{U<:ObjectiveFunction, V<:Int}
    objective               ::U
    length                  ::V
    seqStart                ::Union{Nothing, V}
    seqEnd                  ::Union{Nothing, V}
    precedenceConstraints   ::Bool #true if there are precedence constraints
    precedenceP2F           ::Union{Nothing, Dict{Int, Vector{V}}}#map of priors to followers
    precedenceF2P           ::Union{Nothing, Dict{Int, Vector{V}}} #map of followers to priors
    precedenceNumbers       ::Union{Nothing, Dict{V, V}}#how many things must come before each node
end

Base.length(model::SequencingModel) = model.length

"""
    hasStart(model::SequencingModel)

Check if the sequence has a predefined start

# Arguments
- `model::SequencingModel`: The problem model
"""
function hasStart(model::SequencingModel)
    return !isnothing(model.seqStart)
end

"""
    hasEnd(model::SequencingModel)

Check if the sequence has a predefined end

# Arguments
- `model::SequencingModel`: The problem model
"""
function hasEnd(model::SequencingModel)
    return !isnothing(model.seqEnd)
end


"""
    isBetterSolutionValue(model::SequencingModel, newValue::T, bestKnownValue::T)where{T<:Real}

Check if a solution value is better than the best known value

# Arguments
- `model::SequencingModel`: The problem model
"""
function isBetterSolutionValue(model::SequencingModel, newValue::T, bestKnownValue::T)where{T<:Real}
    if model.objective.type == minimization
        return newValue<bestKnownValue
    else
        return newValue>bestKnownValue
    end
end

"""
    checkValidity(model::SequencingModel, seq::Vector{T})where{T<:Int}

Check if a solution is feasible

# Arguments
- `model::SequencingModel`: The problem model
- `seq::Vector{T}`: The sequence to check
"""
function checkValidity(model::SequencingModel, seq::Vector{T})where{T<:Int}
    #precedence
    if checkPrecedence(model, seq)
        #alldiff
        used = falses(length(model))
        for e in seq
            used[e] = true
        end
        if used == trues(length(model))
            return true
        else
            return false
        end
    else
        return false
    end
end

"""
    checkPrecedence(model::SequencingModel, seq::Vector{T})where{T<:Int}

Check if a solution satisfies all the precedence constraints

# Arguments
- `model::SequencingModel`: The problem model
- `seq::Vector{T}`: The sequence to check
"""
function checkPrecedence(model::SequencingModel, seq::Vector{T})where{T<:Int}
    used = falses(length(seq))
    for e in seq
        used[e] = true
        if haskey(model.precedenceF2P,e)
            for val in values(model.precedenceF2P[e])
                if !used[val]
                    return false
                end
            end
        end
    end
    return true
end
