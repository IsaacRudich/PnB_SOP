function isADivisibleByB(a::S, b::T;digits::Int=8) where {S,T<:Real}
        return ((a%b<=(1/(10^digits))) || (b-(a%b)<=(1/(10^digits))))
end

function roundedMod(a,b)
        return isADivisibleByB(a,b) ? 0 : round(a%b,digits = 8)
end

function roundToEpsilon(a::T;digits::Int=8) where {T<:Real}
        round(a,digits)
end

stop(text="Stopped by Stop Command") = throw(StopException(text))

struct StopException{T}
    S::T
end

function Base.showerror(io::IO, ex::StopException, bt; backtrace=true)
    Base.with_output_color(get(io, :color, false) ? :green : :nothing, io) do io
        showerror(io, ex.S)
    end
end