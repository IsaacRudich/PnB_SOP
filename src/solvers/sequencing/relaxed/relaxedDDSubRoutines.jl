"""
    selectNodes(layer::Union{Vector{RelaxedSequencingNode},Vector{RelaxedStateSequencingNode}}, assignment::T)where{T<:Real}

Find all the nodes in a layer where some value is in someDown, but not in allDown
Return the subset of nodes

# Arguments
- `layer::Union{Vector{RelaxedSequencingNode},Vector{RelaxedStateSequencingNode}}`: the layer to search
- `assignment::T`:the value to check
"""
function selectNodes(layer::Union{Vector{RelaxedSequencingNode},Vector{RelaxedStateSequencingNode}}, assignment::T)where{T<:Real}
    nodes = typeof(layer)()
    if nodes isa Vector{RelaxedSequencingNode}
        for node in layer
            if node.state!=assignment && node.someDown[assignment] && !node.allDown[assignment]
                for parent in getParents(node)
                    if parent.allDown[assignment]
                        push!(nodes, node)
                        break
                    end
                end
            end
        end
    else
        for node in layer
            if node.someDown[assignment] && !node.allDown[assignment]
                if !node.state[assignment]
                    for parent in getParents(node)
                        if parent.allDown[assignment]
                            push!(nodes, node)
                            break
                        end
                    end
                else
                    push!(nodes, node)
                end
            end
        end
    end
    return nodes
end

"""
    splitNode!(node::Union{RelaxedSequencingNode, RelaxedStateSequencingNode}, assignment::U, model::SequencingModel) where{U<:Real}
Split a node during relaxation; updates all the down related variables, but not the up ones

# Arguments
- `node::Union{RelaxedSequencingNode, RelaxedStateSequencingNode}`: the node to split
- `assignment::U`: the equivalency class being created
- `model::SequencingModel`: the sequencing problem being evaluated
"""
function splitNode!(node::Union{RelaxedSequencingNode, RelaxedStateSequencingNode}, assignment::U, model::SequencingModel) where{U<:Real}
    if node isa RelaxedSequencingNode
        newNode = RelaxedSequencingNode(
            node.length, #bestPath
            node.state, #current state
            copy(node.allUp), #allUp
            BitVector(), #allDown
            copy(node.someUp), #someUp
            BitVector(), #someDown
            node.lengthToRoot, #to root
            node.lengthToTerminal, #to terminal
            Vector{Int}(),
            node.exact, #exact
            Vector{ExactArc}(), #arcsin
            Vector{ExactArc}() #arcsout
        )
    else
        newNode = RelaxedStateSequencingNode(
            node.length, #bestPath
            copy(node.state), #current state
            copy(node.allUp), #allUp
            BitVector(), #allDown
            copy(node.someUp), #someUp
            BitVector(), #someDown
            node.lengthToRoot, #to root
            node.lengthToTerminal, #to terminal
            Vector{Int}(),
            node.exact, #exact
            Vector{SimpleArc}(), #arcsin
            Vector{SimpleArc}() #arcsout
        )
    end

    #add new arcs out
    for arc in node.arcsOut
        if node isa RelaxedSequencingNode
            newArc = ExactArc(newNode, arc.destination, arc.assignment, arc.value)
        else
            newArc = SimpleArc(newNode, arc.destination, arc.assignment, arc.value)
        end
        push!(newNode.arcsOut, newArc)
        push!(arc.destination.arcsIn, newArc)
    end

    #redirect each incoming arc where it should go
    if node isa RelaxedSequencingNode
        toRemove = Vector{ExactArc}()
    else
        toRemove = Vector{SimpleArc}()
    end
    for arc in node.arcsIn
        if node isa RelaxedSequencingNode
            if arc.origin.allDown[assignment]
                arcIn = ExactArc(arc.origin, newNode, arc.assignment, arc.value)
                push!(newNode.arcsIn, arcIn)
                push!(arc.origin.arcsOut, arcIn)
                push!(toRemove, arc)
            end
        else
            if arc.origin.allDown[assignment] || arc.assignment == assignment
                arcIn = SimpleArc(arc.origin, newNode, arc.assignment, arc.value)
                push!(newNode.arcsIn, arcIn)
                push!(arc.origin.arcsOut, arcIn)
                push!(toRemove, arc)
            end
        end
    end
    deleteArcs!(toRemove)

    updateNodeDownVariables!(node,model)
    updateNodeDownVariables!(newNode,model)

    return newNode
end

"""
    filterArcs!(arcs::Union{Vector{ExactArc},Vector{SimpleArc}},model::SequencingModel;bestKnownValue::Union{Nothing, U}=nothing,rrbSetting::Union{RRBSetting,Nothing}=nothing, rrbPacket=nothing)where{U<:Real}
Remove arcs that cannot contain the optimal solution

# Arguments
- `arcs::Union{Vector{ExactArc},Vector{SimpleArc}}`: the list of arcs to check
- `model::SequencingModel`: the sequencing problem being evaluated
- `bestKnownValue::Union{Nothing, U}`: the best known solution value
- `rrbSetting::Union{RRBSetting,Nothing}`:The rrb setting to use
- `rrbPacket`:The rrb pre calculated packet
"""
function filterArcs!(arcs::Union{Vector{ExactArc},Vector{SimpleArc}},model::SequencingModel;bestKnownValue::Union{Nothing, U}=nothing,rrbSetting::Union{RRBSetting,Nothing}=nothing, rrbPacket=nothing)where{U<:Real}
    toFilter = typeof(arcs)()
    for i in 1:length(arcs)
        if filterArcRules(arcs[i],model,bestKnownValue=bestKnownValue, rrbSetting=rrbSetting, rrbPacket=rrbPacket)
            push!(toFilter, arcs[i])
        end
    end
    for arc in toFilter
        deleteArc!(arc)
    end
end

"""
    filterArcRules(arc::Union{ExactArc,SimpleArc},model::SequencingModel;bestKnownValue::Union{Nothing, U}=nothing,rrbSetting::Union{RRBSetting,Nothing}=nothing, rrbPacket=nothing)where{U<:Real}
Returns true if an arc cannot contain the optinal solution, returns false otherwise

# Arguments
- `arcs::Union{ExactArc,SimpleArc}`: the  arc to check
- `model::SequencingModel`: the sequencing problem being evaluated
- `bestKnownValue::Union{Nothing, U}`: the best known solution value
- `rrbSetting::Union{RRBSetting,Nothing}`:The rrb setting to use
- `rrbPacket`:The rrb pre calculated packet
"""
function filterArcRules(arc::Union{ExactArc,SimpleArc},model::SequencingModel;bestKnownValue::Union{Nothing, U}=nothing,rrbSetting::Union{RRBSetting,Nothing}=nothing, rrbPacket=nothing)where{U<:Real}
    #the label  is used in every path to the root
    if arc.origin.allDown[arc.assignment]
        return true
    end
    #the label is used in every path to the terminal
    if arc.destination.allUp[arc.assignment]
        return true
    end
    #the number of remaining decisions matches the number of remaining labels, and this one is in use
    if sum(arc.destination.someUp) == (length(model)-length(arc.destination)) && arc.destination.someUp[arc.assignment]
        return true
    end
    #total available paths passing through the arc is less than total number of decisions
    pathsThrough = arc.origin.someDown .| arc.destination.someUp
    pathsThrough[arc.assignment] = true
    if sum(pathsThrough) < length(model)
        return true
    end
    #if the best feasible solution is better than the best optimal path through the arc
    if !isnothing(bestKnownValue)
        if !isBetterSolutionValue(model, arc.origin.lengthToRoot + arc.destination.lengthToTerminal + arc.value, bestKnownValue)
            return true
        end
    end

    #handle precedence constraints in both directions
    if model.precedenceConstraints
        if haskey(model.precedenceF2P,arc.assignment)
            for val in model.precedenceF2P[arc.assignment]
                if !arc.origin.someDown[val]
                    return true
                end
            end
        end
        if haskey(model.precedenceP2F,arc.assignment)
            for val in model.precedenceP2F[arc.assignment]
                if !arc.destination.someUp[val]
                    return true
                end
            end
        end
    end

    #handle rrb
    if !isnothing(rrbSetting) && !isnothing(bestKnownValue) && length(arc.destination)!=length(model)
        if isBetterSolutionValue(model, bestKnownValue, sequencingRRAB(rrbSetting, arc, rrbPacket))
            return true
        end
    end
    return false
end

"""
    bottomUpUpdate(model::SequencingModel, dd::Union{Vector{Vector{RelaxedSequencingNode}},Vector{Vector{RelaxedStateSequencingNode}}}, bestKnownValue::T)where(T<:Real)
Update the upwards variables and delete useless arcs

# Arguments
- `model::SequencingModel`: the sequencing problem being evaluated
- `dd::Union{Vector{Vector{RelaxedSequencingNode}},Vector{Vector{RelaxedStateSequencingNode}}}`: the dd to update
- `bestKnownValue::T`: the best known value, used as a label for infeasible arcs
"""
function bottomUpUpdate(model::SequencingModel, dd::Union{Vector{Vector{RelaxedSequencingNode}},Vector{Vector{RelaxedStateSequencingNode}}}, bestKnownValue::T)where(T<:Real)
    #remove garbage nodes
    for i in Iterators.reverse(2:length(dd)-1)
        for node in dd[i]
            if isempty(node.arcsOut)
                deleteArcs!(node.arcsIn)
            end
        end
        filter!(x->length(x.arcsIn)!=0, dd[i])
    end
    filter!(x->length(x.arcsOut)!=0, dd[1])
    filter!(x->length(x.arcsIn)!=0, dd[length(dd)])

    #update values for remaining nodes
    for i in Iterators.reverse(1:length(dd)-1)
        for node in dd[i]
            updateNodeUpVariables!(node,model)
        end
    end
end

"""
    checkFeasibleTerminalSolutions!(model::SequencingModel, finalNode::RelaxedSequencingNode;bestKnownValue::Union{T,Nothing}=nothing,bestKnownSolution::Union{Vector{U},Nothing}=nothing)where{T<:Real,U<:Int}

Get the best feasible soluution into the terminal node, if there is one better than the input one

# Arguments
- `model::SequencingModel`: the sequencing problem being evaluated
- `finalNode::RelaxedSequencingNode`: the terminal node
- `bestKnownValue::Union{T,Nothing}`: the value of the best known solution
- `bestKnownSolution::Union{Vector{U},Nothing}`: the best known solution
"""
function checkFeasibleTerminalSolutions!(model::SequencingModel, finalNode::RelaxedSequencingNode;bestKnownValue::Union{T,Nothing}=nothing,bestKnownSolution::Union{Vector{U},Nothing}=nothing)where{T<:Real,U<:Int}
    for arc in finalNode.arcsIn
        #if its better and valid, store it
        if isBetterSolutionValue(model, arc.value+arc.origin.lengthToRoot, bestKnownValue) && checkValidity(model, push!(copy(arc.origin.bestSeq),arc.assignment))
            bestKnownValue = arc.value+arc.origin.lengthToRoot
            bestKnownSolution = push!(copy(arc.origin.bestSeq),arc.assignment)
        end
    end
    return bestKnownSolution, bestKnownValue
end