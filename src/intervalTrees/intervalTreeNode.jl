mutable struct IntervalTreeNode{R,S,T<:Real}
        min           ::R
        max           ::S
        stepSize      ::T
        subtreeMax    ::Float64
        rightNode     ::Union{Nothing,IntervalTreeNode}
        leftNode      ::Union{Nothing,IntervalTreeNode}
        parentNode    ::Union{Nothing, IntervalTreeNode}
end

mutable struct Interval{R,S,T<:Real}
        min           ::R
        max           ::S
        stepSize      ::T
end

Interval(min,max) = Interval(min,max,1)
IntervalTreeNode(min,max) = IntervalTreeNode(min, max, 1, convert(Float64, max), nothing, nothing, nothing)
IntervalTreeNode(min,max,stepSize) = IntervalTreeNode(min, max, stepSize, convert(Float64, max), nothing, nothing, nothing)

Base.show(io::IO,node::IntervalTreeNode) = Base.print(io, "{ [",node.min,",",node.max,"]"," stepSize:",node.stepSize," subtreeMax:",node.subtreeMax, " }")

Base.:(==)(x::IntervalTreeNode,y::IntervalTreeNode) = x.min==y.min && x.max==y.max && x.stepSize==y.stepSize ? true : false

Base.isless(x::IntervalTreeNode,y::IntervalTreeNode) = (x.min<y.min ? true : (x.min==y.min && x.max<y.max) ? true : false)
