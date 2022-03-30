mutable struct RelaxedStateSequencingNode{T<:Real,U<:Int,V<:Int}<:Node
    uuid                    ::UUID
    length                  ::U
    state                   ::BitVector
    allUp                   ::BitVector
    allDown                 ::BitVector
    someUp                  ::BitVector
    someDown                ::BitVector
    lengthToRoot            ::T
    lengthToTerminal        ::T
    bestSeq                 ::Vector{V}
    exact                   ::Bool
    arcsIn                  ::Vector{SimpleArc}
    arcsOut                 ::Vector{SimpleArc}
end

#change hashing of nodes to use uuids
Base.hash(node::RelaxedStateSequencingNode) = hash(node.uuid)
Base.:(==)(node1::RelaxedStateSequencingNode, node2::RelaxedStateSequencingNode) = node1.uuid == node2.uuid
Base.length(node::RelaxedStateSequencingNode) = node.length
#do UUIDs automatically
function RelaxedStateSequencingNode(pathLength::U,state::BitVector,allUp::BitVector,allDown::BitVector,someUp::BitVector,someDown::BitVector,lengthToRoot::T,lengthToTerminal::T, bestSeq::Vector{V}, exact::Bool,arcsIn::Vector{SimpleArc},arcsOut::Vector{SimpleArc})where{T<:Real,U<:Int, V<:Int}
    return RelaxedStateSequencingNode(uuid4(), pathLength,state,allUp,allDown,someUp,someDown,lengthToRoot,lengthToTerminal, bestSeq, exact,arcsIn,arcsOut)
end

"""
    deepCopy(node::RelaxedStateSequencingNode)

Deep copy a node

# Arguments
- `node::RelaxedStateSequencingNode`: the node
"""
function deepCopy(node::RelaxedStateSequencingNode)
    return RelaxedStateSequencingNode(
        node.uuid, #uuid
        node.length, #path length
        copy(node.state), #current state
        copy(node.allUp), #allUp
        copy(node.allDown), #allDown
        copy(node.someUp), #someUp
        copy(node.someDown), #someDown
        node.lengthToRoot, #to root
        node.lengthToTerminal, #to terminal
        copy(node.bestSeq),
        node.exact, #exact
        deepCopy(node.arcsIn), #arcsin
        deepCopy(node.arcsOut) #arcsout
    )
end

Base.show(io::IO,node::RelaxedStateSequencingNode) = Base.print(
    io,
    "\n   pathLength: ", node.length, "\n",
    "   state: ", showBitVector(node.state), "\n",
    "   allUp: ", node.allUp, "\n",
    "   allDown: ", node.allDown, "\n",
    "   someUp: ", node.someUp, "\n",
    "   someDown: ", node.someDown, "\n",
    "   lengthToRoot: ", node.lengthToRoot, "\n",
    "   lengthToTerminal: ", node.lengthToTerminal, "\n",
    "   pathToNode: ", node.bestSeq, "\n",
    "   isExact: ", node.exact, "\n",
    "   arcsIn: ", showArcAssignments(node.arcsIn), "\n",
    "   arcsOut: ", showArcAssignments(node.arcsOut), "\n",
)

function showBitVector(b::BitVector)
    str = ""
    for (i,bit) in enumerate(b)
        if bit
            str = string(str, "$i ")
        end
    end
    return str
end

function showArcAssignments(arcs::Vector{SimpleArc})
    str = ""
    for arc in arcs
        str = string(str, arc.assignment, " ")
    end
    return str
end

"""
    createRootRelaxedStateSequencingNodeFromModel(model::SequencingModel)

Create a root node for a relaxed DD from the sequencing model

# Arguments
- `model::SequencingModel`: the sequencing problem to be evaluated
"""
function createRootRelaxedStateSequencingNodeFromModel(model::SequencingModel)
    restrictedNode = createRestrictedRootNode(model)

    state = falses(length(model))
    state[last(restrictedNode.path)] = true
    relaxedNode = RelaxedStateSequencingNode(
        length(restrictedNode.path), #length
        state, #current state
        falses(length(model)), #allUp
        copy(restrictedNode.visited), #allDown
        falses(length(model)), #someUp
        copy(restrictedNode.visited), #someDown
        0, #to root
        0, #to terminal
        copy(restrictedNode.path),
        true, #exact
        Vector{SimpleArc}(), #arcsin
        Vector{SimpleArc}() #arcsout
    )

    return relaxedNode, restrictedNode.domain
end

"""
createPartialRelaxedStateNode(;length::T=0, state::BitVector=falses(1), allUp::BitVector=BitVector(), allDown::BitVector=BitVector(), someUp::BitVector=BitVector(), someDown::BitVector=BitVector(), bestSeq::Vector{U}=Vector{Int}(),lengthToRoot::V=0, lengthToTerminal::V=0, exact::Bool=false, arcsIn::Vector{SimpleArc}= Vector{SimpleArc}(), arcsOut::Vector{SimpleArc}= Vector{SimpleArc}())where{T<:Int, U<:Int, V<:Real}

Create a node using default values for any combination of fields

# Arguments
- `length::T`: the length field
- `state::BitVector`: the state field
- `allUp::BitVector`: the allUp field
- `allDown::BitVector`: the allDown field
- `someUp::BitVector`: the someUp field
- `someDown::BitVector`: the someDown field
- `bestSeq::Vector{Int}`: the bestSeq field
- `lengthToRoot::V`: the lengthToRoot field
- `lengthToTerminal::V`: the lengthToTerminal field
- `exact::Bool`: the exact field
- `arcsIn::Vector{SimpleArc}`: the arcsIn field
- `arcsOut::Vector{SimpleArc}`: the arcsOut field
"""
function createPartialRelaxedStateNode(;length::T=0, state::BitVector=falses(1), allUp::BitVector=BitVector(), allDown::BitVector=BitVector(), someUp::BitVector=BitVector(), someDown::BitVector=BitVector(), bestSeq::Vector{U}=Vector{Int}(),lengthToRoot::V=0, lengthToTerminal::V=0, exact::Bool=false, arcsIn::Vector{SimpleArc}= Vector{SimpleArc}(), arcsOut::Vector{SimpleArc}= Vector{SimpleArc}())where{T<:Int, U<:Int, V<:Real}
    return RelaxedStateSequencingNode(
        length, #length
        state, #current state
        allUp, #allUp
        allDown, #allDown
        someUp, #someUp
        someDown, #someDown
        lengthToRoot, #to root
        lengthToTerminal, #to terminal
        bestSeq,
        exact, #exact
        arcsIn, #arcsin
        arcsOut #arcsout
    )
end

"""
    getParents(node::RelaxedStateSequencingNode)

Find and return the parents of a node

# Arguments
- `node::RelaxedStateSequencingNode`: the node whose parents to find
"""
function getParents(node::RelaxedStateSequencingNode)
    parents = Set{RelaxedStateSequencingNode}()
    for arc in node.arcsIn
        push!(parents, arc.origin)
    end
    return parents
end

"""
    getChildren(node::RelaxedStateSequencingNode)

Find and return the children
# Arguments
- `node::RelaxedSequencingNode`: the node whose children to find
"""
function getChildren(node::RelaxedStateSequencingNode)
    children = Set{RelaxedStateSequencingNode}()
    for arc in node.arcsOut
        push!(children, arc.destination)
    end
    return children
end

"""
    updateNodeDownVariables!(node::RelaxedStateSequencingNode,model::SequencingModel)

Update the state variable for allDown and someDown and lengthToRoot and exact

# Arguments
- `node::RelaxedStateSequencingNode`: the node to update
- `model::SequencingModel`: the sequencing problem being evaluated
"""
function updateNodeDownVariables!(node::RelaxedStateSequencingNode,model::SequencingModel)
    node.allDown = falses(length(model))
    node.someDown = falses(length(model))
    node.state = falses(length(model))
    node.exact = true

    for arc in node.arcsIn
        node.state[arc.assignment] = true
        node.someDown[arc.assignment] = true
        node.allDown[arc.assignment] = true
    end
    #update the best path variables
    bestArc, node.lengthToRoot = findBestArc(model, node.arcsIn)
    node.bestSeq = push!(copy(bestArc.origin.bestSeq), bestArc.assignment)
    #union domains with parents and check for exactness
    for parent in getParents(node)
        node.someDown = node.someDown .| parent.someDown
        node.allDown = node.allDown .| parent.allDown
        if !parent.exact
            node.exact = false
        end
    end
    #check each all down to make sure its valid
    for i in 1:length(node.allDown)
        if node.allDown[i]
            for arc in node.arcsIn
                if arc.assignment != i
                    if !arc.origin.allDown[i]
                        node.allDown[i]=false
                        break
                    end
                end
            end
        end
    end
    #final check for exactness (first part is just a check to see if the second part should run)
    if node.exact && (sum(node.state)!=1 || node.someDown != node.allDown || node.length!=sum(node.allDown))
        node.exact = false
    end
end

"""
    updateOutgoingArcValues!(model::SequencingModel, node::RelaxedStateSequencingNode)

Update the arcs leaving a node
# Arguments
- `model::SequencingModel`: the sequencing problem being evaluated
- `node::RelaxedStateSequencingNode`: the node to update
"""
function updateOutgoingArcValues!(model::SequencingModel, node::RelaxedStateSequencingNode)
    toFilter = Vector{SimpleArc}()
    for arc in node.arcsOut
        tempVal = getArcValue(model, arc)
        if isnothing(tempVal)
            push!(toFilter, arc)
        else
            arc.value = tempVal
        end
    end
    deleteArcs!(toFilter)
end

"""
    updateNodeUpVariables!(node::RelaxedStateSequencingNode,model::SequencingModel)

Update the state variable for someUp, allUp, and lengthToTerminal

# Arguments
- `node::RelaxedStateSequencingNode`: the node to update
- `model::SequencingModel`: the sequencing problem being evaluated
"""
function updateNodeUpVariables!(node::RelaxedStateSequencingNode,model::SequencingModel)
    node.allUp = falses(length(model))
    node.someUp = falses(length(model))
    if isempty(node.arcsOut)
        node.lengthToTerminal = 0
    else
        node.lengthToTerminal = node.arcsOut[1].value + node.arcsOut[1].destination.lengthToTerminal
    end
    for arc in node.arcsOut
        node.someUp[arc.assignment] = true
        node.allUp[arc.assignment] = true
        #update the length to terminal
        if isBetterSolutionValue(model, arc.value + arc.destination.lengthToTerminal, node.lengthToTerminal)
            node.lengthToTerminal = arc.value + arc.destination.lengthToTerminal
        end
    end
    
    #union domains with children
    for child in getChildren(node)
        node.someUp = node.someUp .| child.someUp
        node.allUp = node.allUp .| child.allUp
    end
    #check each all up to make sure its valid
    for i in 1:length(node.allUp)
        if node.allUp[i]
            for arc in node.arcsOut
                if arc.assignment != i
                    if !arc.destination.allUp[i]
                        node.allUp[i]=false
                        break
                    end
                end
            end
        end
    end

    
end


"""
    createWidthOneDDFromNode(model::SequencingModel; rootNode::Union{RelaxedStateSequencingNode,Nothing}=nothing)

Create a width one decision diagram for a relaxed DD from the sequencing model

# Arguments
- `model::SequencingModel`: the sequencing problem to be evaluated
- `rootNode::Union{RelaxedStateSequencingNode,Nothing}`: optional root node, default is creating one from the model
"""
function createWidthOneDDFromNode(model::SequencingModel; rootNode::Union{RelaxedStateSequencingNode,Nothing}=nothing,masterDomain::BitVector=BitVector())
    if isnothing(rootNode)
        rootNode, masterDomain = createRootRelaxedStateSequencingNodeFromModel(model::SequencingModel)
    end

    relaxedDD = Vector{Vector{RelaxedStateSequencingNode}}()
    push!(relaxedDD, [rootNode])

    #for each layer going down
    for layer in 1:length(model)-length(rootNode)
        #new node starts as a copy of the old node
        currentNode = relaxedDD[layer][1]
        newNode = deepCopy(currentNode)
        push!(relaxedDD,[newNode])
        #update length
        newNode.length = currentNode.length+1
        #reset arcs
        newNode.arcsIn = Vector{SimpleArc}()
        newNode.arcsOut = Vector{SimpleArc}()
        #iterate over domain
        for i in 1:length(masterDomain)
            if masterDomain[i]
                #add arcs connecting the two
                newArc = SimpleArc(currentNode, newNode, i, 0)
                tempVal = getArcValue(model, newArc)
                if !isnothing(tempVal)
                    newArc.value = tempVal
                    push!(currentNode.arcsOut, newArc)
                    push!(newNode.arcsIn, newArc)
                end
            end
        end
        updateNodeDownVariables!(newNode,model)
        updateOutgoingArcValues!(model, newNode)

        if model.precedenceConstraints
            for i in 1:length(masterDomain)
                if !masterDomain[i] && !newNode.allDown[i]
                    if haskey(model.precedenceF2P,i)
                        addToDomain = true
                        for j in model.precedenceF2P[i]
                            if !newNode.someDown[j]
                                addToDomain = false
                                break
                            end
                        end
                        if addToDomain
                            masterDomain[i] = true
                        end
                    else
                        masterDomain[i] = true
                    end
                end
            end
        end

        #check for penultimate node
        if hasEnd(model) && layer == length(model)-length(rootNode)-1
            masterDomain = falses(length(masterDomain))
            masterDomain[model.seqEnd] = true
        end
    end

    if model.precedenceConstraints
        for layer in relaxedDD
            for arc in layer[1].arcsOut
                #handle precedence constraints
                if haskey(model.precedenceF2P,arc.assignment)
                    for i in model.precedenceF2P[arc.assignment]
                        #check each of the relevant precedence constraints
                        #if satisfied, put in the domain
                        if !layer[1].someDown[i]
                            deleteArc!(arc)
                            break
                        end
                    end
                end
            end
        end
    end

    #for each layer going up
    for layer in Iterators.reverse(1:length(relaxedDD))
        updateNodeUpVariables!(relaxedDD[layer][1],model)
    end

    return relaxedDD
end
