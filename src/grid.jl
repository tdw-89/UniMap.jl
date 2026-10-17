"""
    parse_int_list("50,100,150") -> sorted unique Vector{Int}
"""
function parse_int_list(s::AbstractString)
    vals = [parse(Int, strip(x)) for x in split(s, ',') if !isempty(strip(x))]
    return sort!(unique!(vals))
end

const ORIENTATIONS = Dict("both" => ['+', '-'], "plus" => ['+'], "minus" => ['-'])

"""
    parse_orientation("both" | "plus" | "minus") -> Vector{Char}

`plus`: R1 from the forward strand at the fragment start; `minus`: R1 = reverse complement
of the fragment end.
"""
function parse_orientation(s::AbstractString)
    o = get(ORIENTATIONS, s, nothing)
    o === nothing && throw(ArgumentError("orientation must be one of: both, plus, minus (got '$s')"))
    return copy(o)
end

function check_lengths(name::AbstractString, v::AbstractVector{<:Integer})
    isempty(v) && throw(ArgumentError("no $name given"))
    all(>(0), v) || throw(ArgumentError("$name must be positive"))
    allunique(v) || throw(ArgumentError("$name must be unique"))
    return sort!(collect(Int, v))
end
