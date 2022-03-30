"""
    removePointFromIntervalTree!(rootNode::IntervalTreeNode, value::T) where T<:Real

Remove a value from the tree by splitting any interval containing that point
Return the new root of the tree

# Arguments
- `rootNode::IntervalTreeNode`: the root of the tree
- `value::T`: the value to be removed
"""
function removePointFromIntervalTree!(rootNode::IntervalTreeNode, value::T) where T<:Real
        oldIntervals= IntervalTreeNode[]
        findIntersectingIntervals!(Interval(value,value,1),rootNode,oldIntervals)
        newIntervals = IntervalTreeNode[]
        for interval in oldIntervals
                splitIntervalOnPoint!(value, newIntervals, interval)
                rootNode = deleteNodeFromIntervalTree!(interval)
        end
        for interval in newIntervals
                if isnothing(rootNode)
                        rootNode = interval
                else
                        addNodeToIntervalTree!(rootNode, interval)
                end
        end
        return rootNode
end

"""
    removeIntervalFromIntervalTree!(rootNode::IntervalTreeNode, interval::Interval)

Remove a value from the tree by splitting any interval instersecting that interval
Return the new root of the tree

# Arguments
- `rootNode::IntervalTreeNode`: the root of the tree
- `interval::Interval`: the interval to be removed
"""
function removeIntervalFromIntervalTree!(rootNode::IntervalTreeNode, interval::Interval)
        #don't do anything on a bad interval
        if interval.min > interval.max
                return
        end

        oldIntervals= IntervalTreeNode[]
        secondPass = IntervalTreeNode[]
        findIntersectingIntervals!(interval,rootNode,oldIntervals)
        newIntervals = IntervalTreeNode[]
        for oldInterval in oldIntervals
                if interval.stepSize == oldInterval.stepSize
                        splitIntervalOnInterval!(interval, newIntervals, oldInterval)
                        rootNode = deleteNodeFromIntervalTree!(oldInterval)
                else
                        push!(secondPass, oldInterval)
                end
        end

        for oldInterval in secondPass
                i = interval.min
                while i<=interval.max
                        rootNode = removePointFromIntervalTree!(rootNode, i)
                        i += interval.stepSize
                end
        end

        for newInterval in newIntervals
                addNodeToIntervalTree!(rootNode, newInterval)
        end

        return rootNode
end


"""
    splitIntervalOnPoint!(value::T, newIntervals::Array{IntervalTreeNode}, oldInterval::IntervalTreeNode) where T<:Real

Split an interval on a point, excludes that point from the new intervals, does not remove the old interval

# Arguments
- `value::T`: the value to split on
- `newIntervals::Array{IntervalTreeNode}`: the list of intervals to add the new intervals to
- `oldInterval::IntervalTreeNode`: the interval to split
"""
function splitIntervalOnPoint!(value::T, newIntervals::Array{IntervalTreeNode}, oldInterval::IntervalTreeNode) where T<:Real
        if !isADivisibleByB(value-oldInterval.min, oldInterval.stepSize)
                #add the left side interval
                if value>oldInterval.min
                        push!(newIntervals, IntervalTreeNode(oldInterval.min, value - roundedMod(value - oldInterval.min, oldInterval.stepSize), oldInterval.stepSize))
                end
                #add the right side of the interval
                if value<oldInterval.max
                        push!(newIntervals, IntervalTreeNode(value + roundedMod(oldInterval.max-value, oldInterval.stepSize),oldInterval.max, oldInterval.stepSize))
                end
        else
                #add the left side interval
                if value>oldInterval.min
                        push!(newIntervals, IntervalTreeNode(oldInterval.min, value-oldInterval.stepSize, oldInterval.stepSize))
                end
                #add the right side of the interval
                if value<oldInterval.max
                        push!(newIntervals, IntervalTreeNode(value+oldInterval.stepSize, oldInterval.max, oldInterval.stepSize))
                end
        end
end

"""
    splitIntervalOnInterval!(interval::Interval, newIntervals::Array{IntervalTreeNode}, oldInterval::IntervalTreeNode)

Split an interval on an interval, excludes that interval from the new intervals, does not remove the old interval

# Arguments
- `interval::Interval`: the interval to split on
- `newIntervals::Array{IntervalTreeNode}`: the list of intervals to add the new intervals to
- `oldInterval::IntervalTreeNode`: the interval to split
"""
function splitIntervalOnInterval!(interval::Interval, newIntervals::Array{IntervalTreeNode}, oldInterval::IntervalTreeNode)
        if !isADivisibleByB(abs(interval.min-oldInterval.min),oldInterval.stepSize)
                #add the left side interval
                if interval.min>oldInterval.min
                        push!(newIntervals, IntervalTreeNode(oldInterval.min, interval.min - roundedMod(interval.min - oldInterval.min,oldInterval.stepSize), oldInterval.stepSize))
                end
                #add the right side of the interval
                if interval.max<oldInterval.max
                        push!(newIntervals, IntervalTreeNode(interval.max + roundedMod(oldInterval.max-interval.max,oldInterval.stepSize),oldInterval.max, oldInterval.stepSize))
                end
        else
                #add the left side interval
                if interval.min>oldInterval.min
                        push!(newIntervals, IntervalTreeNode(oldInterval.min, interval.min-oldInterval.stepSize, oldInterval.stepSize))
                end
                #add the right side of the interval
                if interval.max<oldInterval.max
                        push!(newIntervals, IntervalTreeNode(interval.max+oldInterval.stepSize, oldInterval.max, oldInterval.stepSize))
                end
        end
end

"""
    imposeMaximumOnIntervalTree!(rootNode::IntervalTreeNode, max::T) where T<:Real

Remove all values above the input max value from the tree

# Arguments
- `rootNode::IntervalTreeNode`: the root node
- `max::T`: the maximum value to impose
"""
function imposeMaximumOnIntervalTree!(rootNode::IntervalTreeNode, max::T) where T<:Real
        if findIntervalTreeMax(rootNode)>max
                intersections = IntervalTreeNode[]
                findIntersectingIntervals!(Interval(max, findIntervalTreeMax(rootNode),1), rootNode, intersections)
                for intersection in intersections
                        rootNode = deleteNodeFromIntervalTree!(intersection)
                        if intersection.min <= max
                                addNodeToIntervalTree!(rootNode, IntervalTreeNode(intersection.min,max - roundedMod(max-intersection.min,intersection.stepSize),intersection.stepSize))
                        end
                end
        end
        return rootNode
end

"""
    imposeMinimumOnIntervalTree!(rootNode::IntervalTreeNode, min::T) where T<:Real

Remove all values below the input min value from the tree

# Arguments
- `rootNode::IntervalTreeNode`: the root node
- `min::T`: the minimum value to impose
"""
function imposeMinimumOnIntervalTree!(rootNode::IntervalTreeNode, min::T) where T<:Real
        if findIntervalTreeMin(rootNode)<min
                intersections = IntervalTreeNode[]
                findIntersectingIntervals!(Interval(findIntervalTreeMin(rootNode), min,1), rootNode, intersections)
                for intersection in intersections
                        rootNode = deleteNodeFromIntervalTree!(intersection)
                        if intersection.max >= min
                                addNodeToIntervalTree!(rootNode, IntervalTreeNode(min + roundedMod(intersection.max-min,intersection.stepSize),intersection.max ,intersection.stepSize))
                        end
                end
        end
        return rootNode
end

"""
    imposeOpenMaximumOnIntervalTree!(rootNode::IntervalTreeNode, max::T) where T<:Real

Remove all values >= the input max value from the tree

# Arguments
- `rootNode::IntervalTreeNode`: the root node
- `max::T`: the maximum value to impose
"""
function imposeOpenMaximumOnIntervalTree!(rootNode::IntervalTreeNode, max::T) where T<:Real
        if findIntervalTreeMax(rootNode)>=max
                intersections = IntervalTreeNode[]
                findIntersectingIntervals!(Interval(max, findIntervalTreeMax(rootNode),1), rootNode, intersections)
                for intersection in intersections
                        rootNode = deleteNodeFromIntervalTree!(intersection)
                        if intersection.min < max
                                if !isADivisibleByB(max-intersection.min, intersection.stepSize)
                                        addNodeToIntervalTree!(rootNode, IntervalTreeNode(intersection.min,max - roundedMod(max-intersection.min,intersection.stepSize),intersection.stepSize))
                                else
                                        addNodeToIntervalTree!(rootNode, IntervalTreeNode(intersection.min,max - intersection.stepSize,intersection.stepSize))
                                end
                        end
                end
        end
        return rootNode
end

"""
    imposeOpenMinimumOnIntervalTree!(rootNode::IntervalTreeNode, min::T) where T<:Real

Remove all values <= the input min value from the tree

# Arguments
- `rootNode::IntervalTreeNode`: the root node
- `min::T`: the minimum value to impose
"""
function imposeOpenMinimumOnIntervalTree!(rootNode::IntervalTreeNode, min::T) where T<:Real
        if findIntervalTreeMin(rootNode)<=min
                intersections = IntervalTreeNode[]
                findIntersectingIntervals!(Interval(findIntervalTreeMin(rootNode), min,1), rootNode, intersections)
                for intersection in intersections
                        rootNode = deleteNodeFromIntervalTree!(intersection)
                        if intersection.max > min
                                if !isADivisibleByB(intersection.max-min, intersection.stepSize)
                                        addNodeToIntervalTree!(rootNode, IntervalTreeNode(min + roundedMod(intersection.max-min,intersection.stepSize),intersection.max ,intersection.stepSize))
                                else
                                        addNodeToIntervalTree!(rootNode, IntervalTreeNode(min + intersection.stepSize,intersection.max ,intersection.stepSize))
                                end
                        end
                end
        end
        return rootNode
end

"""
    addIntervalToIntervalTreeWithMerging!(rootNode::IntervalTreeNode,newNode::Interval)

Add an interval to the tree, but first remove all nodes it overlaps with and merge them into one node
Returns the root node (which may have changed)

# Arguments
- `rootNode::IntervalTreeNode`: the root node
- `newNode::Interval`: the Interval to add to the tree
"""
function addIntervalToIntervalTreeWithMerging!(rootNode::IntervalTreeNode,newNode::Interval)
        intersections = IntervalTreeNode[]
        findIntersectingIntervals!(Interval(newNode.min,newNode.max,newNode.stepSize),rootNode,intersections)

        for intersection in intersections
                if ((intersection.stepSize == newNode.stepSize) && isADivisibleByB(abs(intersection.min-newNode.min),newNode.stepSize))
                        rootNode = deleteNodeFromIntervalTree!(intersection)
                        newNode.min = min(newNode.min,intersection.min)
                        newNode.max = max(newNode.max, intersection.max)
                end
        end
        if rootNode !== nothing
                addNodeToIntervalTree!(rootNode, IntervalTreeNode(newNode.min, newNode.max, newNode.stepSize))
        else
                rootNode = IntervalTreeNode(newNode.min, newNode.max, newNode.stepSize)
        end
        return rootNode
end

"""
    deepCopyIntervalTreeNodes(rootNode::Union{IntervalTreeNode,Nothing}, isMerging::Bool)

Deep copy an interval tree, but the order of the nodes is not maintained (the new tree will still be a proper ordering though)

# Arguments
- `rootNode::Union{IntervalTreeNode,Nothing}`: the root node
- `isMerging::Bool`: if true, the copied tree will eliminate repeat information as it copies
"""
function deepCopyIntervalTreeNodes(rootNode::Union{IntervalTreeNode,Nothing}, isMerging::Bool)
        if rootNode === nothing
                return nothing
        end

        return recursiveIntervalTreeCopy!(rootNode,IntervalTreeNode(rootNode.min, rootNode.max, rootNode.stepSize),isMerging,true)
end


"""
    recursiveIntervalTreeCopy!(currentNode::Union{Nothing, IntervalTreeNode}, newRoot::IntervalTreeNode,isMerging::Bool,isRoot::Bool)

An internal method for recursively copying an interval tree onto a specified root node

# Arguments
- `currentNode::Union{Nothing, IntervalTreeNode}`: the current node to copy
- `newRoot::IntervalTreeNode`: the root of the copy
- `isMerging::Bool`: if true, the copy will merge overlapping nodes
- `isRoot::Bool`: External calls should always be true, inside the method it is always false
"""
function recursiveIntervalTreeCopy!(currentNode::Union{Nothing, IntervalTreeNode}, newRoot::IntervalTreeNode,isMerging::Bool,isRoot::Bool)
        if isRoot
                if isnothing(currentNode.leftNode) && isnothing(currentNode.rightNode)
                        if isMerging
                                newRoot = addIntervalToIntervalTreeWithMerging!(newRoot, Interval(currentNode.min, currentNode.max,currentNode.stepSize))
                        else
                                addNodeToIntervalTree!(newRoot, IntervalTreeNode(currentNode.min, currentNode.max,currentNode.stepSize))
                        end
                        return newRoot
                else
                        if isMerging
                                newRoot = addIntervalToIntervalTreeWithMerging!(newRoot, Interval(currentNode.min, currentNode.max,currentNode.stepSize))
                        else
                                addNodeToIntervalTree!(newRoot, IntervalTreeNode(currentNode.min, currentNode.max,currentNode.stepSize))
                        end
                        newRoot = recursiveIntervalTreeCopy!(currentNode.leftNode,newRoot, isMerging,false)
                        newRoot = recursiveIntervalTreeCopy!(currentNode.rightNode,newRoot, isMerging,false)
                        return newRoot
                end
        elseif currentNode !== nothing
                if !isMerging
                        addNodeToIntervalTree!(newRoot, IntervalTreeNode(currentNode.min, currentNode.max,currentNode.stepSize))
                else
                        newRoot = addIntervalToIntervalTreeWithMerging!(newRoot, Interval(currentNode.min, currentNode.max,currentNode.stepSize))
                end
                newRoot = recursiveIntervalTreeCopy!(currentNode.leftNode,newRoot, isMerging,false)
                newRoot = recursiveIntervalTreeCopy!(currentNode.rightNode,newRoot, isMerging,false)
        else
                return newRoot
        end
end

"""
    unionIntervalTrees(rootA::IntervalTreeNode, rootB::IntervalTreeNode, isMerging::Bool)

Return the union of two interval trees

# Arguments
- `rootA::IntervalTreeNode`: the root of one of the trees
- `rootB::IntervalTreeNode`: the root of a different tree
- `isMerging::Bool`: if true, the union will merge overlapping nodes
"""
function unionIntervalTrees(rootA::IntervalTreeNode, rootB::IntervalTreeNode, isMerging::Bool)
        return recursiveIntervalTreeCopy!(rootA, deepCopyIntervalTreeNodes(rootB, isMerging), isMerging, true)
end


"""
    intersectIntervalTrees(rootA::Union{Nothing, IntervalTreeNode}, rootB::Union{IntervalTreeNode, Nothing})

Return the intersection of two interval trees
O(m log|n|)  where m is the size of tree A and n is the size of tree B
Significantly faster when the step sizes of overlapping intervals are the same

# Arguments
- `rootA::IntervalTreeNode`: the root of one of the trees
- `rootB::IntervalTreeNode`: the root of a different tree
"""
function intersectIntervalTrees(rootA::Union{Nothing, IntervalTreeNode}, rootB::Union{IntervalTreeNode, Nothing})
        return recursiveIntervalTreeIntersect!(rootA, nothing,rootB)
end

"""
    recursiveIntervalTreeIntersect!(currentNode::Union{Nothing, IntervalTreeNode}, newRoot::Union{Nothing, IntervalTreeNode},rootB::Union{Nothing, IntervalTreeNode})

Internal method that recursively constructs the intersection of two interval trees
O(m log|n|)  where m is the size of tree A and n is the size of tree B
Significantly faster when the step sizes of overlapping intervals are the same

# Arguments
- `currentNode::Union{Nothing, IntervalTreeNode}`: the current node from treeA to test
- `newRoot::Union{Nothing, IntervalTreeNode}`: the resulting interval tree, if calling externally make the argument: nothing
- `rootB::IntervalTreeNode`: the root of a different tree
"""
function recursiveIntervalTreeIntersect!(currentNode::Union{Nothing, IntervalTreeNode}, newRoot::Union{Nothing, IntervalTreeNode},rootB::Union{Nothing, IntervalTreeNode})
        if currentNode === nothing
                return newRoot
        else
                intersections = IntervalTreeNode[]
                findIntersectingIntervals!(Interval(currentNode.min,currentNode.max,currentNode.stepSize),rootB,intersections)
                for intersection in intersections
                        #if they overlap with the same step size, and the steps line up
                        if ((currentNode.stepSize == intersection.stepSize) && isADivisibleByB(abs(intersection.min-currentNode.min),currentNode.stepSize))
                                if newRoot !== nothing
                                        newRoot = addIntervalToIntervalTreeWithMerging!(newRoot, Interval(max(currentNode.min, intersection.min),min(currentNode.max, intersection.max),currentNode.stepSize))
                                else
                                        newRoot = IntervalTreeNode(max(currentNode.min, intersection.min),min(currentNode.max, intersection.max),currentNode.stepSize)
                                end
                        #else if they have different step sizes
                        elseif (currentNode.stepSize != intersection.stepSize)
                                #increment both intervals until you find the lowest intersection point
                                min1 = currentNode.min
                                min2 = intersection.min
                                while (min1 != min2) && (min1 <= min(currentNode.max, intersection.max)) && (min2 <= min(currentNode.max, intersection.max))
                                        if min1 < min2
                                                min1+=currentNode.stepSize
                                                min1 = round(min1, digits = 8)

                                        else
                                                min2+= intersection.stepSize
                                                min2 = round(min2, digits = 8)
                                        end
                                end
                                #if there is a valid starting point
                                if (min1 == min2) && (min1<=min(currentNode.max, intersection.max)) && (min2 <= min(currentNode.max, intersection.max))
                                        #find the gap that allows for overlap (like lCM but with fractions allowed)
                                        gap = max(currentNode.stepSize, intersection.stepSize)
                                        while (!isADivisibleByB(gap, min(currentNode.stepSize, intersection.stepSize)) && (min1+gap <= min(currentNode.max, intersection.max)))
                                                gap += max(currentNode.stepSize, intersection.stepSize)
                                        end
                                        #if there is a valid gap
                                        if min1+gap <= min(currentNode.max, intersection.max)
                                                if newRoot !== nothing
                                                        newRoot = addIntervalToIntervalTreeWithMerging!(newRoot, Interval(min1,min1+(trunc((min(currentNode.max, intersection.max)-min1)/gap)*gap),gap))
                                                else
                                                        newRoot = IntervalTreeNode(min1,min1+(trunc((min(currentNode.max, intersection.max)-min1)/gap)*gap),gap)
                                                end
                                        else # enter in a single point
                                                if newRoot !== nothing
                                                        newRoot = addIntervalToIntervalTreeWithMerging!(newRoot, Interval(min1,min1,min(currentNode.stepSize,intersection.stepSize)))
                                                else
                                                        newRoot = IntervalTreeNode(min1,min1,min(currentNode.stepSize,intersection.stepSize))
                                                end
                                        end
                                end
                        end
                end
                newRoot = recursiveIntervalTreeIntersect!(currentNode.leftNode,newRoot, rootB)
                return recursiveIntervalTreeIntersect!(currentNode.rightNode,newRoot, rootB)
        end
end

"""
    disjunctiveUnion(rootA::IntervalTreeNode, rootB::IntervalTreeNode)

Return the disjunctive union of two interval trees as an ordered touple.
The first elements of the touple contains the elements not in rootA, the second contains the elements not in rootB

# Arguments
- `rootA::IntervalTreeNode`: the root of one of the trees
- `rootB::IntervalTreeNode`: the root of a different tree
"""
function disjunctiveUnion(rootA::IntervalTreeNode, rootB::IntervalTreeNode)
        intersections = intervalTreeToSet(intersectIntervalTrees(rootA,rootB))
        #create disjoint set of rootA
        disjointA = Set{Real}()
        elements = intervalTreeToSet(rootB)
        for e in elements
                if !(e in intersections)
                        push!(disjointA, e)
                end
        end
        #create disjoint set of rootB
        disjointB = Set{Real}()
        elements = intervalTreeToSet(rootA)
        for e in elements
                if !(e in intersections)
                        push!(disjointB, e)
                end
        end

        return (disjointA, disjointB)
end
