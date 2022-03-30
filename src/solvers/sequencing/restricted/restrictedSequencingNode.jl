struct RestrictedSequencingNode{T<:Real}<:Node
    path                    ::Vector{Int}
    visited                 ::BitVector
    value                   ::T
    domain                  ::BitVector
end

Base.show(io::IO,node::RestrictedSequencingNode) = Base.print(
    io,
    "\n   path: ", node.path, "\n",
    "   value: ", node.value, "\n",
    "   visited: ", node.visited, "\n",
    "   domain: ", node.domain, "\n"
)

"""
    deepCopy(node::RestrictedSequencingNode)

Make a deep copy of a RestrictedSequencingNode

# Arguments
- `node::RestrictedSequencingNode`: The node to copy
"""
function deepCopy(node::RestrictedSequencingNode)
    return RestrictedSequencingNode(
        copy(node.path),
        copy(node.visited),
        node.value,
        copy(node.domain)
    )
end

"""
    getState(node::RestrictedSequencingNode)

Get the last location visited in the path leading to the given node

# Arguments
- `node::RestrictedSequencingNode`: The node 
"""
function getState(node::RestrictedSequencingNode)
    return last(node.path)
end


"""
    inDomain(node::RestrictedSequencingNode, i::Int)

Check if a value is in the domain of a node (check if a location has been visited or not)

# Arguments
- `node::RestrictedSequencingNode`: The node whose domain is being checked
- `i::Int`: The location/value to check for in the domain
"""
function inDomain(node::RestrictedSequencingNode, i::Int)
    return node.domain[i]
end

"""
    remove!(node::RestrictedSequencingNode, i::Int)

Remove a value from the domain of a node

# Arguments
- `node::RestrictedSequencingNode`: the node
- `i::Int`: the value to remove
"""
function remove!(node::RestrictedSequencingNode, i::Int)
    node.domain[i] = false
end

"""
    add!(node::RestrictedSequencingNode, i::Int)

Add a a value to the domain of a node

# Arguments
- `node::RestrictedSequencingNode`: the domain
- `i::Int`: the value to add
"""
function add!(node::RestrictedSequencingNode, i::Int)
    node.domain[i] = true
end

"""
    initializeDomain(model::SequencingModel)

Initalize the domain for a given sequencing problem

# Arguments
- `model::SequencingModel`: the problem instance
"""
function initializeDomain(model::SequencingModel)
    d = trues(model.length)

    #handle precedence constraints
    for follower in keys(model.precedenceF2P)
        d[follower]=false
    end

    return d
end

"""
    createRestrictedRootNode(model::SequencingModel)

Create the root node for a problem using the given model

# Arguments
- `model::SequencingModel`: The problem instance
"""
function createRestrictedRootNode(model::SequencingModel)
    if hasStart(model)
        v = falses(model.length)
        v[model.seqStart] = true
        d = initializeDomain(model)
        d[model.seqStart] = false
        return RestrictedSequencingNode(Int[model.seqStart], v,0, d)
    else
        return RestrictedSequencingNode(Int[], falses(model.length),0, initializeDomain(model))
    end
end
