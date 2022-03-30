struct DPModel{U<:ObjectiveFunction}
    variables               ::Array{String}
    objective               ::U
    constraints             ::Dict{Union{DataType,UnionAll}, Array{Any}}
    domains                 ::Dict{String, Union{IntervalTreeNode,Nothing}}
end

"""
    deepCopy(domains::Union{Dict{String, IntervalTreeNode}, Dict{String, Union{Nothing, IntervalTreeNode}}})

Make a deep copy of a domain object

# Arguments
- `domains::Union{Dict{String, IntervalTreeNode}, Dict{String, Union{Nothing, IntervalTreeNode}}}`: the domains to copy
"""
function deepCopy(domains::Union{Dict{String, IntervalTreeNode}, Dict{String, Union{Nothing, IntervalTreeNode}}})
    newDomains = Dict{String, Union{IntervalTreeNode,Nothing}}()
    for (key, value) in domains
        if !isnothing(value)
            newDomains[key] = deepCopyIntervalTreeNodes(value, false)
        end
    end
    return newDomains
end
