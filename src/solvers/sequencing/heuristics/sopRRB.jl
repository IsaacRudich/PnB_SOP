"""
    precalculateSOPRRB(model::SequencingModel)

A method that preloads and returns the calculations required to find a rough relaxed bound for SOPs

# Arguments
- `model::SequencingModel`: The SOP
"""
function precalculateSOPRRB(model::SequencingModel)
    paramType = valtype(model.objective.f.components)
    values = Dict{Int, Vector{Tuple{Int, paramType}}}()
    lastArc = Vector{Tuple{Int,paramType}}()
    for (states,objValue) in model.objective.f.components
        startState = states[1]
        endState = states[2]
        if haskey(values,startState)
            insert!(
                values[startState],
                searchsorted(values[startState], (endState, objValue), by = x -> x[2]).start,
                (endState, objValue)
            )
        else
            values[startState] = [(endState, objValue)]
        end
        if endState == model.seqEnd
            insert!(
                lastArc,
                searchsorted(lastArc, (startState, objValue), by = x -> x[2]).start,
                (startState, objValue)
            )
        end
    end
    return values, lastArc
end

"""
    sopRRB(visited::BitVector, state::U, startingValue::T, values::Dict{Int, Vector{Tuple{U, T}}},lastArc::Vector{Tuple{U, T}})where{T<:Real,U<:Int}

A method that takes in the current state of the sop and the precalculated values and returns a rough relaxe bound

# Arguments
- `visited::BitVector`: The visited nodes
- `state::U`: The current node
- `values::Dict{U, Vector{Tuple{U, T}}}`: The presorted out arc values
- `lastArc::Vector{Tuple{Int, T}}`: The presorted in arc values for the last edge
"""
function sopRRB(visited::BitVector, state::U, startingValue::T, values::Dict{Int, Vector{Tuple{U, T}}},lastArc::Vector{Tuple{U, T}})where{T<:Real,U<:Int}
    rrb = startingValue
    #add the first arc
    for pair in values[state]
        if !visited[pair[1]] && pair[1]!=length(visited)
            rrb += pair[2]
            break
        end
    end
    #add the last arc
    for pair in lastArc
        if !visited[pair[1]]
            rrb += pair[2]
            break
        end
    end

    worstValue = 0
    for i in 1:length(visited)-1
        if !visited[i]
            for pair in values[i]
                if !visited[pair[1]] && pair[1]!=length(visited)
                    rrb += pair[2]
                    if pair[2] > worstValue
                        worstValue = pair[2]
                    end
                    break
                end
            end
        end
    end
    rrb -= worstValue
    return rrb
end

"""
    sopRRAB(arc::Union{ExactArc,SimpleArc},values::Dict{Int, Vector{Tuple{U, T}}},lastArc::Vector{Tuple{U, T}})where{T<:Real,U<:Int}

A method that takes in an arc from a relaxed DD of the sop, and the precalculated values, then returns a rough relaxed bound

# Arguments
- `arc::Union{ExactArc,SimpleArc}`: The arc to value
- `values::Dict{U, Vector{Tuple{U, T}}}`: The presorted out arc values
- `lastArc::Vector{Tuple{Int, T}}`: The presorted in arc values for the last edge
"""
function sopRRAB(arc::Union{ExactArc,SimpleArc},values::Dict{Int, Vector{Tuple{U, T}}},lastArc::Vector{Tuple{U, T}})where{T<:Real,U<:Int}
    rrb = arc.value
    visited = copy(arc.origin.allDown)
    visited[arc.assignment] = true
    #add the first arc
    for pair in values[arc.assignment]
        if !visited[pair[1]] && pair[1]!=length(visited)
            rrb += pair[2]
            break
        end
    end
    #add the last arc
    for pair in lastArc
        if !visited[pair[1]]
            rrb += pair[2]
            break
        end
    end

    worstValues = Vector{Int}()
    for i in 1:length(visited)-1
        if !visited[i]
            valueFound = false
            for pair in values[i]
                if !visited[pair[1]] && pair[1]!=length(visited)
                    if isempty(worstValues)
                        push!(worstValues, pair[2])
                    else
                        insert!(
                            worstValues,
                            searchsorted(worstValues, pair[2]).start,
                            pair[2]
                        )
                    end
                    valueFound = true
                    break
                end
            end
            if !valueFound
                if isempty(worstValues)
                    push!(worstValues, 0)
                else
                    insert!(
                        worstValues,
                        searchsorted(worstValues, 0).start,
                        0
                    )
                end
            end
        end
    end
    worstValues = worstValues[1:(length(visited)-length(arc.destination)-2)]
    rrb += sum(worstValues)
    return rrb
end
