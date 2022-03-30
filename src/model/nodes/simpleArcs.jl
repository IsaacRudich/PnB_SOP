mutable struct SimpleArc{T<:Node,U<:Real,V<:Int}
    uuid                    ::UUID
    origin                  ::T
    destination             ::T
    assignment              ::V
    value                   ::U        
end

Base.show(io::IO,arc::SimpleArc) = Base.print(
    io,
    "\n   Assignment: ", arc.assignment, "\n",
    "   Arc Value: ", arc.value, "\n"
)

#change hashing of nodes to use uuids
Base.hash(arc::SimpleArc) = hash(arc.uuid)
Base.:(==)(x::SimpleArc,y::SimpleArc) = x.uuid == y.uuid

#do UUIDs automatically
function SimpleArc(origin::T, destination::T, assignment::U, value::V)where{T<:Node, U<:Int, V<:Real}
    return SimpleArc(uuid4(), origin, destination, assignment, value)
end

"""
    deepCopy(arc::SimpleArc)

Deep copy an arc

# Arguments
- `arcs::SimpleArc`: the arc
"""
function deepCopy(arc::SimpleArc)
    return SimpleArc(arc.uuid, arc.origin,arc.destination,arc.assignment,arc.value)
end

"""
    deepCopy(arcs::Vector{SimpleArc})

Deep copy a list of arcs

# Arguments
- `arcs::Array{SimpleArc}`: a list of arcs
"""
function deepCopy(arcs::Vector{SimpleArc})
    newArcs = Vector{SimpleArc}()
    for arc in arcs
        push!(newArcs, deepCopy(arc))
    end
    return newArcs
end

"""
    shallowCopy(arc::SimpleArc)

Shallow copy an arc

# Arguments
- `arcs::SimpleArc`: the arc
"""
function shallowCopy(arc::SimpleArc)
    return SimpleArc(arc.origin,arc.destination,arc.assignment,arc.value)
end

"""
    shallowCopy(arcs::Vector{SimpleArc})

functionally the same as deepCopy a list of arcs, but with different uuids

# Arguments
- `arcs::Array{SimpleArc}`: a list of arcs
"""
function shallowCopy(arcs::Vector{SimpleArc})
    newArcs = Vector{SimpleArc}()
    for arc in arcs
        push!(newArcs, shallowCopy(arc))
    end
    return newArcs
end

"""
    findBestArc(model::SequencingModel, arcs::Vector{SimpleArc})

Find the best arc in a list, where best means containing the most optimal path

# Arguments
- `model::SequencingModel`: The Sequencing problem
- `arcs::Vector{SimpleArc}`: The list of arcs
"""
function findBestArc(model::SequencingModel, arcs::Vector{SimpleArc})
    bestArc = arcs[1]
    bestArcValue = bestArc.origin.lengthToRoot + bestArc.value
    for arc in arcs
        if isBetterSolutionValue(model, arc.origin.lengthToRoot + arc.value, bestArcValue)
            bestArc = arc
            bestArcValue = arc.origin.lengthToRoot + arc.value
        end
    end
    return bestArc, bestArcValue
end

"""
    getArcValue(model::SequencingModel, arc::SimpleArc)

Get the best arc value by checking the states it can be coming from, arc.origin.state must be a BitVector

# Arguments
- `model::SequencingModel`: The Sequencing problem
- `arc::SimpleArc`: The arc
"""
function getArcValue(model::SequencingModel, arc::SimpleArc)
    arcValue = nothing
    for (i, bit) in enumerate(arc.origin.state)
        if bit
            if hasKey(model.objective, i, arc.assignment)
                if isnothing(arcValue) || isBetterSolutionValue(model, evaluateDecision(model.objective, i, arc.assignment),arcValue)
                    arcValue = evaluateDecision(model.objective, i, arc.assignment)
                end
            end
        end
    end
    return arcValue
end

"""
    deleteArc!(arc::SimpleArc)
Remove an arc from both relevant nodes

# Arguments
- `arc::SimpleArc`: the arc to remove
"""
function deleteArc!(arc::SimpleArc)
    filter!(x->x!=arc, arc.origin.arcsOut)
    filter!(x->x!=arc, arc.destination.arcsIn)
end

"""
    deleteArcs!(arcsToRemove::Vector{SimpleArc})
Delete a list of arcs

# Arguments
- `arcsToRemove::Vector{SimpleArc}`: the arcs to remove
"""
function deleteArcs!(arcsToRemove::Vector{SimpleArc})
    toRemove = Vector{SimpleArc}()
    for arc in arcsToRemove
        push!(toRemove, arc)
    end
    for arc in toRemove
        deleteArc!(arc)
    end
end