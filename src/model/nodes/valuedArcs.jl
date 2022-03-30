mutable struct ValuedArc{T<:Node,U<:Real,V<:Int}
    uuid                    ::UUID
    origin                  ::T
    destination             ::T
    assignment              ::V
    value                   ::Union{Nothing, U}
    seq                     ::Union{Nothing, Vector{V}}
    seqValue                ::Union{Nothing, U}
end

Base.show(io::IO,arc::ValuedArc) = Base.print(
    io,
    "\n   Assignment: ", arc.assignment, "\n",
    "   Arc Value: ", arc.value, "\n",
    "   Sequence: ", arc.seq, "\n",
    "   Sequence Value: ", arc.seqValue, "\n"
)

#change hashing of nodes to use uuids
Base.hash(arc::ValuedArc) = hash(arc.uuid)
Base.:(==)(x::ValuedArc,y::ValuedArc) = x.uuid == y.uuid
#do UUIDs automatically
function ValuedArc(origin::T, destination::T, assignment::U, value::V,seq::Vector{W}, seqValue::V)where{T<:Node, U<:Real, V<:Real,W<:Int}
    return ValuedArc(uuid4(), origin, destination, assignment, value,seq,seqValue)
end
function ValuedArc(origin::T, destination::T, assignment::U)where{T<:Node, U<:Real}
    return ValuedArc{T,U,Int}(uuid4(), origin, destination, assignment, nothing,nothing,nothing)
end

"""
    deepCopy(arc::ValuedArc)

Deep copy an arc

# Arguments
- `arcs::ValuedArc`: the arc
"""
function deepCopy(arc::ValuedArc)
    return ValuedArc(arc.uuid, arc.origin,arc.destination,arc.assignment,arc.value,arc.seq,arc.seqValue)
end

"""
    deepCopy(arcs::Vector{ValuedArc})

Deep copy a list of arcs

# Arguments
- `arcs::Array{ValuedArc}`: a list of arcs
"""
function deepCopy(arcs::Vector{ValuedArc})
    newArcs = Vector{ValuedArc}()
    for arc in arcs
        push!(newArcs, deepCopy(arc))
    end
    return newArcs
end

"""
    findBestArcPath!(model::SequencingModel, arc::ValuedArc)

Find the best value and path an arc can take by iterating over the possible origin sequences
Return the value of the arc, or nothing if no feasible sequence includes this arc

# Arguments
- `model::SequencingModel`: The Sequencing problem
- `arc::ValuedArc`: The arc
"""
function findBestArcPath!(model::SequencingModel, arc::ValuedArc)
    bestVal = nothing
    bestSeqVal = nothing
    bestPath = nothing
    for parentArc in arc.origin.arcsIn
        if hasKey(model.objective, parentArc.assignment, arc.assignment)
            if isnothing(bestSeqVal) || isBetterSolutionValue(model, evaluateDecision(model.objective, parentArc.assignment, arc.assignment)+parentArc.seqValue, bestSeqVal)
                bestVal = evaluateDecision(model.objective, parentArc.assignment, arc.assignment)
                bestSeqVal = bestVal + parentArc.seqValue
                bestPath = copy(parentArc.seq)
            end
        end
    end
    if !isnothing(bestSeqVal)
        push!(bestPath,arc.assignment)
        arc.seq = bestPath
        arc.value = bestVal
        arc.seqValue = bestSeqVal
        return arc.value
    else
        return nothing
    end
end

"""
    findBestArc(model::SequencingModel, arcs::Vector{ValuedArc})

Find the best arc in a list

# Arguments
- `model::SequencingModel`: The Sequencing problem
- `arcs::Vector{ValuedArc}`: The list of arcs
"""
function findBestArc(model::SequencingModel, arcs::Vector{ValuedArc})
    bestArc = arcs[1]
    for arc in arcs
        if isBetterSolutionValue(model, arc.seqValue, bestArc.seqValue)
            bestArc = arc
        end
    end
    return bestArc
end
