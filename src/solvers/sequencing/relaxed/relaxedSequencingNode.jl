mutable struct RelaxedSequencingNode{T<:Real,U<:Int,V<:Int}<:Node
    uuid                    ::UUID
    length                  ::U
    state                   ::V
    allUp                   ::BitVector
    allDown                 ::BitVector
    someUp                  ::BitVector
    someDown                ::BitVector
    lengthToRoot            ::T
    lengthToTerminal        ::T
    bestSeq                 ::Vector{V}
    exact                   ::Bool
    arcsIn                  ::Vector{ExactArc}
    arcsOut                 ::Vector{ExactArc}
end

#change hashing of nodes to use uuids
Base.hash(node::RelaxedSequencingNode) = hash(node.uuid)
Base.:(==)(node1::RelaxedSequencingNode, node2::RelaxedSequencingNode) = node1.uuid == node2.uuid
Base.length(node::RelaxedSequencingNode) = node.length
#do UUIDs automatically
function RelaxedSequencingNode(pathLength::U,state::V,allUp::BitVector,allDown::BitVector,someUp::BitVector,someDown::BitVector,lengthToRoot::T,lengthToTerminal::T,bestSeq::Vector{V}, exact::Bool,arcsIn::Vector{ExactArc},arcsOut::Vector{ExactArc})where{T<:Real,U<:Int,V<:Int}
    return RelaxedSequencingNode(uuid4(), pathLength,state,allUp,allDown,someUp,someDown,lengthToRoot,lengthToTerminal,bestSeq, exact,arcsIn,arcsOut)
end

"""
    deepCopy(node::RelaxedSequencingNode)

Deep copy a node

# Arguments
- `node::RelaxedSequencingNode`: the node
"""
function deepCopy(node::RelaxedSequencingNode)
    return RelaxedSequencingNode(
        node.uuid, #uuid
        node.length, #path length
        node.state, #current state
        copy(node.allUp), #allUp
        copy(node.allDown), #allDown
        copy(node.someUp), #someUp
        copy(node.someDown), #someDown
        node.lengthToRoot, #to root
        node.lengthToTerminal, #to terminal
        copy(node.bestSeq), #best Seq
        node.exact, #exact
        deepCopy(node.arcsIn), #arcsin
        deepCopy(node.arcsOut) #arcsout
    )
end

Base.show(io::IO,node::RelaxedSequencingNode) = Base.print(
    io,
    "\n   pathLength: ", node.length, "\n",
    "   state: ", node.state, "\n",
    "   allUp: ", node.allUp, "\n",
    "   allDown: ", node.allDown, "\n",
    "   someUp: ", node.someUp, "\n",
    "   someDown: ", node.someDown, "\n",
    "   lengthToRoot: ", node.lengthToRoot, "\n",
    "   lengthToTerminal: ", node.lengthToTerminal, "\n",
    "   best Sequence: ", node.bestSeq, "\n",
    "   isExact: ", node.exact, "\n",
    "   arcsIn: ", showInArcOrigins(node.arcsIn), "\n",
    "   arcsOut: ", showOutArcAssignments(node.arcsOut), "\n",
)

function showOutArcAssignments(arcs::Vector{ExactArc})
    str = ""
    for arc in arcs
        str = string(str, arc.assignment, " ")
    end
    return str
end

function showInArcOrigins(arcs::Vector{ExactArc})
    str = ""
    for arc in arcs
        str = string(str, arc.origin.state, " ")
    end
    return str
end

"""
    createRootRelaxedSequencingNodeFromModel(model::SequencingModel)

Create a root node for a relaxed DD from the sequencing model

# Arguments
- `model::SequencingModel`: the sequencing problem to be evaluated
"""
function createRootRelaxedSequencingNodeFromModel(model::SequencingModel)
    restrictedNode = createRestrictedRootNode(model)
    relaxedNode = RelaxedSequencingNode(
        length(restrictedNode.path), #length
        last(restrictedNode.path), #current state
        falses(length(model)), #allUp
        copy(restrictedNode.visited), #allDown
        falses(length(model)), #someUp
        copy(restrictedNode.visited), #someDown
        0, #to root
        0, #to terminal
        copy(restrictedNode.path),
        true, #exact
        Vector{ExactArc}(), #arcsin
        Vector{ExactArc}() #arcsout
    )
    return relaxedNode, restrictedNode.domain
end

"""
    createPartialRelaxedNode(;length::T=0, state::U=-1, allUp::BitVector=BitVector(), allDown::BitVector=BitVector(), someUp::BitVector=BitVector(), someDown::BitVector=BitVector(), bestSeq::Vector{U}=Vector{Int}(), lengthToRoot::V=0, lengthToTerminal::V=0, exact::Bool=false, arcsIn::Vector{ExactArc}= Vector{ExactArc}(), arcsOut::Vector{ExactArc}= Vector{ExactArc}())where{T<:Int, U<:Int, V<:Real}

Create a node using default values for any combination of fields

# Arguments
- `length::T`: the length field
- `state::U`: the state field
- `allUp::BitVector`: the allUp field
- `allDown::BitVector`: the allDown field
- `someUp::BitVector`: the someUp field
- `someDown::BitVector`: the someDown field
- `bestSeq::Vector{Int}`: the bestSeq field
- `lengthToRoot::V`: the lengthToRoot field
- `lengthToTerminal::V`: the lengthToTerminal field
- `exact::Bool`: the exact field
- `arcsIn::Vector{ExactArc}`: the arcsIn field
- `arcsOut::Vector{ExactArc}`: the arcsOut field
"""
function createPartialRelaxedNode(;length::T=0, state::U=-1, allUp::BitVector=BitVector(), allDown::BitVector=BitVector(), someUp::BitVector=BitVector(), someDown::BitVector=BitVector(), bestSeq::Vector{U}=Vector{Int}(), lengthToRoot::V=0, lengthToTerminal::V=0, exact::Bool=false, arcsIn::Vector{ExactArc}= Vector{ExactArc}(), arcsOut::Vector{ExactArc}= Vector{ExactArc}())where{T<:Int, U<:Int, V<:Real}
    return RelaxedSequencingNode(
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
    getParents(node::RelaxedSequencingNode)

Find and return the parents of a node

# Arguments
- `node::RelaxedSequencingNode`: the node whose parents to find
"""
function getParents(node::RelaxedSequencingNode)
    parents = Set{RelaxedSequencingNode}()
    for arc in node.arcsIn
        push!(parents, arc.origin)
    end
    return parents
end

"""
    getChildren(node::RelaxedSequencingNode)

Find and return the children
# Arguments
- `node::RelaxedSequencingNode`: the node whose children to find
"""
function getChildren(node::RelaxedSequencingNode)
    children = Set{RelaxedSequencingNode}()
    for arc in node.arcsOut
        push!(children, arc.destination)
    end
    return children
end

"""
    updateNodeDownVariables!(node::RelaxedSequencingNode,model::SequencingModel)

Update the state variable for allDown and someDown and lengthToRoot and exact

# Arguments
- `node::RelaxedSequencingNode`: the node to update
- `model::SequencingModel`: the sequencing problem being evaluated
"""
function updateNodeDownVariables!(node::RelaxedSequencingNode,model::SequencingModel)
    node.allDown = falses(length(model))
    node.someDown = falses(length(model))
    node.exact = true

    node.someDown[node.state] = true
    node.allDown[node.state] = true
    #update the best path variables
    bestArc, node.lengthToRoot = findBestArc(model, node.arcsIn)
    node.bestSeq = push!(copy(bestArc.origin.bestSeq), node.state)
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
    if node.exact && (node.someDown != node.allDown || node.length!=sum(node.allDown))
        node.exact = false
    end
end

"""
    updateNodeUpVariables!(node::RelaxedSequencingNode,model::SequencingModel)

Update the state variable for someUp, allUp, and lengthToTerminal

# Arguments
- `node::RelaxedSequencingNode`: the node to update
- `model::SequencingModel`: the sequencing problem being evaluated
"""
function updateNodeUpVariables!(node::RelaxedSequencingNode,model::SequencingModel)
    node.allUp = falses(length(model))
    node.someUp = falses(length(model))
    node.lengthToTerminal = node.arcsOut[1].value + node.arcsOut[1].destination.lengthToTerminal
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
    createSquareStartDD(model::SequencingModel;rootNode::Union{Nothing, RelaxedSequencingNode}=nothing, masterDomain::BitVector=BitVector())

Create an initial dd from the model where all the arcs have exact values

# Arguments
- `model::SequencingModel`: the sequencing problem being evaluated
- `rootNode::Union{Nothing, RelaxedSequencingNode}`: a root node if not starting from the model
- `masterDomain::BitVector`: the domain of the root node if there is a root node, must be included if root node is included
"""
function createSquareStartDD(model::SequencingModel;rootNode::Union{Nothing, RelaxedSequencingNode}=nothing, masterDomain::BitVector=BitVector())
    if isnothing(rootNode)
        rootNode, masterDomain = createRootRelaxedSequencingNodeFromModel(model)
    end
    relaxedDD = Vector{Vector{RelaxedSequencingNode}}()
    push!(relaxedDD, [rootNode])

    #for each layer going down
    terminalNode = createPartialRelaxedNode(
        length = length(model),
        state=0,
        allUp=falses(length(masterDomain)),
        allDown=falses(length(masterDomain)),
        someUp=falses(length(masterDomain)),
        someDown=falses(length(masterDomain))
    )
    if hasEnd(model)
        lastLayerLoopIndex = length(model)-length(rootNode)
        terminalNode.state = model.seqEnd
    else
        lastLayerLoopIndex = length(model)-length(rootNode)+1
    end
    for layer in 2:lastLayerLoopIndex
        push!(relaxedDD,Vector{RelaxedSequencingNode}())
        #new node starts as a state and length
        for (index, inDomain) in enumerate(masterDomain)
            if inDomain
                newNode = createPartialRelaxedNode(length = layer, state=index)
                push!(relaxedDD[layer],newNode)
            end
        end
        #update domain based on precedence
        if model.precedenceConstraints
            for (i, inDomain) in enumerate(masterDomain)
                if !inDomain && (!hasStart(model) || i!=model.seqStart)
                    #add satisfied precedence constraints to the domain
                    if layer >= get(model.precedenceNumbers,i,0)
                        masterDomain[i] = true
                    end
                end
            end
        end
    end#layer iteration
    push!(relaxedDD, [terminalNode])

    #add all arcs except terminal arcs and update down values
    for layer in 1:length(relaxedDD)-2
        for parent in relaxedDD[layer]
            if layer != 1
                updateNodeDownVariables!(parent,model)
            end
            for child in relaxedDD[layer+1]
                if hasKey(model.objective, parent.state, child.state) && !parent.allDown[child.state]
                    newArc = ExactArc(parent, child, child.state, evaluateDecision(model.objective, parent.state, child.state))
                    push!(parent.arcsOut, newArc)
                    push!(child.arcsIn, newArc)
                end
            end
        end
        toRemove = Vector{RelaxedSequencingNode}()
        for child in relaxedDD[layer+1]
            if isempty(child.arcsIn)
                push!(toRemove, child)
            end
        end
        filter!(x->!(x in toRemove), relaxedDD[layer+1])
    end
    #add arcs to terminal
    if hasEnd(model)
        for parent in relaxedDD[length(relaxedDD)-1]
            if hasKey(model.objective, parent.state, model.seqEnd)
                newArc = ExactArc(parent, terminalNode, model.seqEnd, evaluateDecision(model.objective, parent.state, model.seqEnd))
                push!(parent.arcsOut, newArc)
                push!(terminalNode.arcsIn, newArc)
            end
            updateNodeDownVariables!(parent,model)
        end
    else
        for parent in relaxedDD[length(relaxedDD)-1]
            newArc = ExactArc(parent, terminalNode, -1, 0)
            push!(parent.arcsOut, newArc)
            push!(terminalNode.arcsIn, newArc)
            updateNodeDownVariables!(parent,model)
        end
    end
    updateNodeDownVariables!(terminalNode,model)

    #update all the up variables
    for layer in Iterators.reverse(1:length(relaxedDD)-1)
        for node in relaxedDD[layer]
            updateNodeUpVariables!(node,model)
        end
    end

    return relaxedDD
end