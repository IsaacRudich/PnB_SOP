struct ExtensionalArcFunction{T<:Real}
    components      ::Dict{Tuple{Int, Int}, T}
end

struct ExtensionalArcObjective<:ObjectiveFunction
    f           ::ExtensionalArcFunction
    type        ::ObjectiveType
end

"""
    evaluateDecision(obj::ExtensionalArcObjective, i1::Int, i2::Int)

Takes in a decision and returns the impact on the function

# Arguments
- `obj::ExtensionalArcObjective`: the function to be evaluated
- `i1::Int`: the index of the first point
- `i2::Int`: the index of the second point
"""
function evaluateDecision(obj::ExtensionalArcObjective, i1::Int, i2::Int)
    return obj.f.components[(i1, i2)]
end

"""
    evaluateDecisions(obj::ExtensionalArcObjective,decisions::Vector{Int})

Takes in a list of decisions and returns the impact on the function

# Arguments
- `obj::ExtensionalArcObjective`: the function to be evaluated
- `decisions::Vector{Int}`: the list
"""
function evaluateDecisions(obj::ExtensionalArcObjective,decisions::Vector{Int})
    sum = 0
    for i in 1:length(decisions)-1
        sum += evaluateDecision(obj, decisions[i], decisions[i+1])
    end
    return sum
end
"""
    hasKey(obj::ExtensionalArcObjective,i1::Int, i2::Int)

Takes in an arc and return its validity

# Arguments
- `obj::ExtensionalArcObjective`: the function to be evaluated
- `i1::Int`: the index of the first point
- `i2::Int`: the index of the second point
"""
function hasKey(obj::ExtensionalArcObjective,i1::Int, i2::Int)
    return haskey(obj.f.components,(i1, i2))
end
