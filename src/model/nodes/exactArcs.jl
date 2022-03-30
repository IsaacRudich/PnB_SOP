struct ExactArc{T<:Node,U<:Real,V<:Int}
    uuid                    ::UUID
    origin                  ::T
    destination             ::T
    assignment              ::V
    value                   ::U        
end

Base.show(io::IO,arc::ExactArc) = Base.print(
    io,
    "\n   Assignment: ", arc.assignment, "\n",
    "   Arc Value: ", arc.value, "\n"
)

#change hashing of nodes to use uuids
Base.hash(arc::ExactArc) = hash(arc.uuid)
Base.:(==)(x::ExactArc,y::ExactArc) = x.uuid == y.uuid

#do UUIDs automatically
function ExactArc(origin::T, destination::T, assignment::U, value::V)where{T<:Node, U<:Int, V<:Real}
    return ExactArc(uuid4(), origin, destination, assignment, value)
end

"""
    deepCopy(arc::ExactArc)

Deep copy an arc

# Arguments
- `arcs::ExactArc`: the arc
"""
function deepCopy(arc::ExactArc)
    return ExactArc(arc.uuid, arc.origin,arc.destination,arc.assignment,arc.value)
end

"""
    deepCopy(arcs::Vector{ExactArc})

Deep copy a list of arcs

# Arguments
- `arcs::Array{ExactArc}`: a list of arcs
"""
function deepCopy(arcs::Vector{ExactArc})
    newArcs = Vector{ExactArc}()
    for arc in arcs
        push!(newArcs, deepCopy(arc))
    end
    return newArcs
end

"""
    shallowCopy(arc::ExactArc)

Shallow copy an arc

# Arguments
- `arcs::ExactArc`: the arc
"""
function shallowCopy(arc::ExactArc)
    return ExactArc(arc.origin,arc.destination,arc.assignment,arc.value)
end

"""
    shallowCopy(arcs::Vector{ExactArc})

functionally the same as deepCopy a list of arcs, but with different uuids

# Arguments
- `arcs::Array{ExactArc}`: a list of arcs
"""
function shallowCopy(arcs::Vector{ExactArc})
    newArcs = Vector{ExactArc}()
    for arc in arcs
        push!(newArcs, shallowCopy(arc))
    end
    return newArcs
end

"""
    findBestArc(model::SequencingModel, arcs::Vector{ExactArc})

Find the best arc in a list, where best means containing the most optimal path

# Arguments
- `model::SequencingModel`: The Sequencing problem
- `arcs::Vector{ExactArc}`: The list of arcs
"""
function findBestArc(model::SequencingModel, arcs::Vector{ExactArc})
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
    deleteArc!(arc::ExactArc)
Remove an arc from both relevant nodes

# Arguments
- `arc::ExactArc`: the arc to remove
"""
function deleteArc!(arc::ExactArc)
    filter!(x->x!=arc, arc.origin.arcsOut)
    filter!(x->x!=arc, arc.destination.arcsIn)
end

"""
    deleteArcs!(arcsToRemove::Vector{ExactArc})
Delete a list of arcs

# Arguments
- `arcsToRemove::Vector{ExactArc}`: the arcs to remove
"""
function deleteArcs!(arcsToRemove::Vector{ExactArc})
    toRemove = Vector{ExactArc}()
    for arc in arcsToRemove
        push!(toRemove, arc)
    end
    for arc in toRemove
        deleteArc!(arc)
    end
end