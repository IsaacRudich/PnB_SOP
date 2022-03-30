"""
    addNodeToIntervalTree!(parentNode::IntervalTreeNode, newNode::IntervalTreeNode)

Add a node to the tree recursively

# Arguments
- `parentNode::IntervalTreeNode`: the root node
- `newNode::IntervalTreeNode`: the node to add to the tree
"""
function addNodeToIntervalTree!(parentNode::IntervalTreeNode, newNode::IntervalTreeNode)
        #don't do anything if the added node is invalid
        if newNode.min > newNode.max || newNode.stepSize <=0
                return
        end

        #fix the max of the node before adding it
        if !isADivisibleByB(newNode.max-newNode.min, newNode.stepSize)
                newNode = IntervalTreeNode(newNode.min, newNode.max - roundedMod(newNode.max-newNode.min,newNode.stepSize), newNode.stepSize)
        end
        #update subtreeMaximum as you insert
        if newNode.max > parentNode.subtreeMax
                parentNode.subtreeMax = newNode.subtreeMax
        end

        #add node to the tree recursively
        if parentNode < newNode
                if parentNode.rightNode === nothing
                        parentNode.rightNode = newNode
                        newNode.parentNode = parentNode
                else
                        addNodeToIntervalTree!(parentNode.rightNode, newNode)
                end
        elseif parentNode == newNode
                return
        else
                if parentNode.leftNode === nothing
                        parentNode.leftNode = newNode
                        newNode.parentNode = parentNode
                else
                        addNodeToIntervalTree!(parentNode.leftNode, newNode)
                end
        end

        return
end

"""
    findIntersectingIntervals!(query::Interval,currentNode::Union{Nothing, IntervalTreeNode},intersections::Array{IntervalTreeNode})

Recursively find all intersecting intervals

# Arguments
- `query::Interval`: the interval to check against
- `currentNode::Union{Nothing, IntervalTreeNode}`: the node the recursive function is calling (start with root)
- `intersections::Array{IntervalTreeNode}`: the list of intersecting intervals
"""
function findIntersectingIntervals!(query::Interval,currentNode::Union{Nothing, IntervalTreeNode},intersections::Array{IntervalTreeNode})
        #don't do anything if the query is invalid
        if query.min > query.max || query.stepSize <=0
                return
        end

        if currentNode === nothing
                return
        end

        if !( (currentNode.min > query.max) || (currentNode.max < query.min) )
                push!(intersections,currentNode)
        end

        if (currentNode.leftNode !== nothing) && (currentNode.leftNode.subtreeMax >= query.min)
                findIntersectingIntervals!(query, currentNode.leftNode, intersections)
        end

        findIntersectingIntervals!(query, currentNode.rightNode, intersections)
end


"""
    deleteNodeFromIntervalTree!(node::IntervalTreeNode)

Recursively find all intersecting intervals and remove the intersection

# Arguments
- `node::IntervalTreeNode`: the node to be removed from the tree
"""
function deleteNodeFromIntervalTree!(node::IntervalTreeNode)
        #if the node has no children
        if node.leftNode === nothing && node.rightNode === nothing
                #if the node has no children and has parents
                if node.parentNode !== nothing
                        #if the node is a left node
                        if isLeftNodeInIntervalTree(node)
                                #update the parent node to no longer refer to the deleted node
                                node.parentNode.leftNode = nothing
                        else
                                #update the parent node to no longer refer to the deleted node
                                node.parentNode.rightNode = nothing
                        end
                        updateSubtreeMax!(node.parentNode)
                        return findRootOfIntervalTree(node.parentNode)
                #if the node has no children and no parents
                else
                        #return an empty tree
                        return nothing
                end
        #if the node has no right node and has a left node
        elseif node.rightNode === nothing && node.leftNode !== nothing
                #if the node has a parent
                if node.parentNode !== nothing
                        #if the node is a left node
                        if isLeftNodeInIntervalTree(node)
                                #update the left node of the parent to be the left node of the deleted node
                                node.parentNode.leftNode = node.leftNode
                        else
                                #update the right node of the parent to be the left node of the deleted node
                                node.parentNode.rightNode = node.leftNode
                        end
                        #update the parent node of the left node to be the parent of the deleted node
                        node.leftNode.parentNode = node.parentNode
                        #update the subtree value
                        updateSubtreeMax!(node.parentNode)
                        #return the root of the tree
                        return findRootOfIntervalTree(node.parentNode)
                #if the node does not have a parent
                else
                        #then return the node to the left because it is the new root
                        node.leftNode.parentNode = nothing
                        return node.leftNode
                end
        #if the node has no left node and has a right node
        elseif node.leftNode === nothing && node.rightNode !== nothing
                #if the node has a parent
                if node.parentNode !== nothing
                        #if the node is a left node
                        if isLeftNodeInIntervalTree(node)
                                #update the left node of the parent to be the right node of the deleted node
                                node.parentNode.leftNode = node.rightNode
                        else
                                #update the right node of the parent to be the right node of the deleted node
                                node.parentNode.rightNode = node.rightNode
                        end
                        #update the parent node of the rightnode to be the parent of the deleted node
                        node.rightNode.parentNode = node.parentNode
                        #update the subtree value
                        updateSubtreeMax!(node.parentNode)
                        #return the root of the tree
                        return findRootOfIntervalTree(node.parentNode)
                else
                        #then return the node to the right because it is the new root
                        node.rightNode.parentNode = nothing
                        return node.rightNode
                end
        #if the node has no two children
        else
                #find the nodes replacement
                successor = findSuccessorNode(node.rightNode)

                #remove refrences to the replacement from the tree
                deleteNodeFromIntervalTree!(successor)
                #the successors left node becomes the same as the node's left node
                successor.leftNode = node.leftNode
                #the successor has a left child, update that child's parent
                successor.leftNode.parentNode = successor

                #the successors right node becomes the same as the node's right node
                successor.rightNode = node.rightNode
                if successor.rightNode !== nothing
                        #if the successor has a right child, update that child's parent
                        successor.rightNode.parentNode = successor
                end

                #the successor takes the nodes parent
                successor.parentNode = node.parentNode
                #if the successor has a parent
                if successor.parentNode !== nothing
                        #update that parent
                        if isLeftNodeInIntervalTree(node)
                                node.parentNode.leftNode = successor
                        else
                                node.parentNode.rightNode = successor
                        end
                end
                #update the subtreeMax
                updateSubtreeMax!(successor)
                #return the root of the tree
                return findRootOfIntervalTree(successor)
        end
end


"""
    findSuccessorNode(node::IntervalTreeNode)

Find the leftmost node in the subtree,
Should be called on a right node [Ex => findSuccessorNode(node.rightNode)] when being used for deletion

# Arguments
- `node::IntervalTreeNode`: the node whose successor we want to find
"""
function findSuccessorNode(node::IntervalTreeNode)
        if node.leftNode === nothing
                return node
        else
                return findSuccessorNode(node.leftNode)
        end
end

"""
    updateSubtreeMax!(node::Union{Nothing, IntervalTreeNode})

Update the subtreeMax of the tree recursively from the chosen node up through its parents

# Arguments
- `node::Union{Nothing, IntervalTreeNode}`: the node whose parents, grandparents, etc... you want to update
"""
function updateSubtreeMax!(node::Union{Nothing, IntervalTreeNode})
        if node !== nothing
                if node.leftNode === nothing && node.rightNode === nothing
                        node.subtreeMax = node.max
                elseif node.leftNode === nothing
                        node.subtreeMax = max(node.max,node.rightNode.subtreeMax)
                elseif node.rightNode === nothing
                        node.subtreeMax = max(node.max,node.leftNode.subtreeMax)
                else
                        node.subtreeMax = max(node.max,node.leftNode.subtreeMax, node.rightNode.subtreeMax)
                end
                updateSubtreeMax!(node.parentNode)
        end
end

"""
    isLeftNodeInIntervalTree(node::IntervalTreeNode)

Return true if it is a left node and false otherwise

# Arguments
- `node::IntervalTreeNode`: the node to check
"""
function isLeftNodeInIntervalTree(node::IntervalTreeNode)
        if node.parentNode.leftNode === node
                return true
        else
                return false
        end
end

"""
    findRootOfIntervalTree(node::IntervalTreeNode)

Return the root of the tree

# Arguments
- `node::IntervalTreeNode`: the node whose root you want to find
"""
function findRootOfIntervalTree(node::IntervalTreeNode)
        if node.parentNode===nothing
                return node
        else
                return findRootOfIntervalTree(node.parentNode)
        end
end

"""
    findIntervalTreeMax(node::IntervalTreeNode)

Return the maximum value in the domain

# Arguments
- `node::IntervalTreeNode`: any node in the tree
"""
function findIntervalTreeMax(node::IntervalTreeNode)
        return findRootOfIntervalTree(node).subtreeMax
end

"""
    findIntervalTreeMin(node::IntervalTreeNode)

Return the minimum value in the domain

# Arguments
- `node::IntervalTreeNode`: any node in the tree
"""
function findIntervalTreeMin(node::IntervalTreeNode)
        return findSuccessorNode(findRootOfIntervalTree(node)).min
end

"""
    findAllStepSizes!(node::IntervalTreeNode, stepSizes::Set{T}) where T<:Real

Mutate the set to add all step sizes

# Arguments
- `node::IntervalTreeNode`: the root node of the tree
- `stepSizes::Set{T}`: the set to be mutated
"""
function findAllStepSizes!(node::IntervalTreeNode, stepSizes::Set{T}) where T<:Real
        if node !== nothing
                push!(stepSizes, node.stepSizes)
        end
        if node.rightNode !== nothing
                findAllStepSizes!(node.rightNode, stepSizes)
        end
        if node.leftNode !== nothing
                findAllStepSizes!(node.leftNode, stepSizes)
        end
end

"""
    intervalTreeToSet(rootNode::IntervalTreeNode) where T<:Real

Create a Set{T} containing every value from the interval tree

# Arguments
- `rootNode::IntervalTreeNode`: the root node of the tree
"""
function intervalTreeToSet(rootNode::Union{IntervalTreeNode,Nothing}) where T<:Real
        if isnothing(rootNode)
                return Set{Float64}()
        else
                return recursiveIntervalToSet!(Set{Real}(), rootNode)
        end
end

"""
    recursiveIntervalToSet!(domain::Set{Any}, node::Union{Nothing, IntervalTreeNode})

Return a Set{T} containing every value from the interval tree constructed recursively

# Arguments
- `domain::Set{T}`: the set to add the values to
- `rootNode::IntervalTreeNode`: the root node of the tree
"""
function recursiveIntervalToSet!(domain::Set{T}, node::Union{Nothing, IntervalTreeNode}) where T<:Real
        if node !== nothing
                current = node.min
                while current <= node.max
                        push!(domain, current)
                        current += node.stepSize
                end
                recursiveIntervalToSet!(domain, node.rightNode)
                recursiveIntervalToSet!(domain, node.leftNode)
                return domain
        end
        return domain
end

"""
    setToIntervalTree(domain::Union{Array{T},Set{T}}) where T<:Real

Return an IntervalTreeNode containing ranges that represent the input set in a memory efficient way

# Arguments
- `domain::Union{Array{T},Set{T}}`: The value to be put in the set
"""
function setToIntervalTree(domain::Union{Array{T},Set{T}}) where T<:Real
        if isa(domain, Set{T})
                listForm = collect(domain)
        else
                listForm = unique!(domain)
        end

        sort!(listForm)
        rootNode = nothing
        if length(listForm)>=3
                min = listForm[1]
                stepSize = listForm[2]-listForm[1]
                for i in 3:length(listForm)
                        if listForm[i]-listForm[i-1]==stepSize
                                if i<length(listForm)
                                        continue
                                else
                                        if rootNode === nothing
                                                rootNode = IntervalTreeNode(min, listForm[i], roundToEpsilon(stepSize))
                                        else
                                                addNodeToIntervalTree!(rootNode, IntervalTreeNode(min, listForm[i], roundToEpsilon(stepSize)))
                                        end
                                end
                        else
                                if min == listForm[i-1]
                                        stepSize = 1
                                end
                                if rootNode === nothing
                                        rootNode = IntervalTreeNode(min, listForm[i-1], roundToEpsilon(stepSize))
                                else
                                        addNodeToIntervalTree!(rootNode, IntervalTreeNode(min, listForm[i-1], roundToEpsilon(stepSize)))
                                end
                                min = listForm[i]
                                if i < length(listForm)
                                        stepSize = listForm[i+1]-listForm[i]
                                else
                                        addNodeToIntervalTree!(rootNode, IntervalTreeNode(min, min, 1))
                                end
                        end
                end
        elseif length(listForm)==2
                rootNode = IntervalTreeNode(listForm[1], listForm[2], listForm[2]-listForm[1])
        elseif length(listForm)==1
                rootNode = IntervalTreeNode(listForm[1], listForm[1], 1)
        end

        return rootNode
end

"""
    printTree(node::Union{Nothing, IntervalTreeNode})

Print the interval tree from bottom left to parent to right to left ...

# Arguments
- `node::Union{Nothing, IntervalTreeNode}`: The tree being printed
"""
function printTree(node::Union{Nothing, IntervalTreeNode})
    if node === nothing
        return
    end

    if node.leftNode !== nothing
        printTree(node.leftNode)
    end

    print(node, "\n")

    if node.rightNode !== nothing
        printTree(node.rightNode)
    end
end
