mutable struct Sequence{T<:Real}
    sequenceComponents      ::Array{String}
    sequenceValues          ::Array{<:Real}
    value                   ::T
end

Sequence() = Sequence(String[],Float64[],0)

"""
    addToSequence!(seq::Sequence, decVar::String, decision::T, addedValue::R) where {R<:Real,T<:Real}

Add a new part to a stored sequence

# Arguments
- `seq::Sequence`: the sequence to be added to
- `decVar::String`: the sequence member being decided
- `decision::T`: the value of the sequence component
- `addedValue::R`:the cost/benefit of adding the new part to the sequence
"""
function addToSequence!(seq::Sequence, decVar::String, decision::T, addedValue::R) where {R<:Real,T<:Real}
    push!(seq.sequenceComponents, decVar)
    push!(seq.sequenceValues, decision)
    seq.value = seq.value + addedValue
end

"""
    deepCopySequence(seq::Sequence)

Returns a deep copy of a sequence

# Arguments
- `seq::Sequence`: the sequence to be deep copied
"""
function deepCopy(seq::Sequence)
    return Sequence(copy(seq.sequenceComponents),copy(seq.sequenceValues),seq.value)
end

"""
    sequenceToDict(seq::Sequence)

Returns a Dict version of the sequence object

# Arguments
- `seq::Sequence`: the sequence to be translated
"""
function sequenceToDict(seq::Sequence)
    sequence = Dict{String, Float64}()
    for i in 1:length(seq.sequenceComponents)
        sequence[seq.sequenceComponents[i]]=seq.sequenceValues[i]
    end
    return sequence
end

"""
    findBestSequence(list::Array{Sequence},largestValue::Bool)

Returns the best sequence from a list as defined by the smallest or largest value

# Arguments
- `list::Array{Sequence}`: the list to choose a sequence from
- `largestValue::Bool`: return the sequence with the largest value if true, smallest value if false
"""
function findBestSequence(list::Array{Sequence},largestValue::Bool)
    returnSeq = nothing
    if largestValue
        for seq in list
            if isnothing(returnSeq) || seq.value>returnSeq.value
                returnSeq = seq
            end
        end
    else
        for seq in list
            if isnothing(returnSeq) || seq.value<returnSeq.value
                returnSeq = seq
            end
        end
    end
    return returnSeq
end

"""
    getSequenceLength(seq::Sequence)

Returns the length of a sequence

# Arguments
- `seq::Sequence`: the sequence
"""
function getSequenceLength(seq::Sequence)
    return length(seq.sequenceComponents)
end

"""
    getLastDecision(seq::Sequence)

Returns the value of the last decision in the sequence

# Arguments
- `seq::Sequence`: the sequence
"""
function getLastDecision(seq::Sequence)
    return seq.sequenceValues[length(seq.sequenceValues)]
end

"""
    getValues(seq::Sequence)

Get the values in the sequence in order

# Arguments
- `seq::Sequence`: the sequence
"""
function getValues(seq::Sequence)
    return seq.sequenceValues
end
"""
    getValue(seq::Sequence)

Get the value of the sequence

# Arguments
- `seq::Sequence`: the sequence
"""
function getValue(seq::Sequence)
    return seq.value
end
