"""
    checkPrecedence(visited::BitVector,priors::Vector{Int})

Return true if the precedence constraint is satisfied

# Arguments
- `visited::BitVector`: A list of visited nodes
- `priors::Vector{Int}`: A list of nodes that must have been visited
"""
function checkPrecedence(visited::BitVector,priors::Vector{Int})
    for i in priors
        if !visited[i]
            return false
        end
    end
    return true
end

"""
    getPrecedenceNumbers(precedenceF2P::Dict{T, Vector{T}})where{T<:Int}

Return a list of pre-calculated precedence numbers

# Arguments
- `precedenceF2P::Dict{T, Vector{T}}`: The map of followers to their priors
"""
function getPrecedenceNumbers(precedenceF2P::Dict{T, Vector{T}})where{T<:Int}
    precedenceNums = Dict{Int, Int}()
    for num in keys(precedenceF2P)
        if !haskey(precedenceNums,num)
            getPrecedenceNumber!(num, precedenceF2P,precedenceNums)
        end
    end
    return precedenceNums
end

"""
    getPrecedenceNumber!(checkNum::T, precedenceF2P::Dict{T, Vector{T}},precedenceNums::Dict{Int, Int})where{T<:Int}

Return a  precedence number, but also updates the values of any missing numbers it has to find along the way

# Arguments
- `checkNum::T`: The number to calculate
- `precedenceF2P::Dict{T, Vector{T}}`: The map of followers to their priors
- `precedenceNums::Dict{Int, Int}`: The list of calculated numbers
"""
function getPrecedenceNumber!(checkNum::T, precedenceF2P::Dict{T, Vector{T}},precedenceNums::Dict{Int, Int})where{T<:Int}
    if haskey(precedenceF2P,checkNum)
        maxVal = 0
        for element in precedenceF2P[checkNum]
            if !haskey(precedenceNums,element)
                precedenceNums[element] = getPrecedenceNumber!(element, precedenceF2P,precedenceNums)
            end
            if precedenceNums[element]>maxVal
                maxVal = precedenceNums[element]
            end
        end
        precedenceNums[checkNum] = max(maxVal+1,length(precedenceF2P[checkNum]))
    else
        precedenceNums[checkNum] = 0
    end
end
