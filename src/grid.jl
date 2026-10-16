"""
    fragment_length_grid(fmin, fmax) -> Vector{Int}

Fragment lengths from `fmin`, with a step that starts at 1 and doubles every two
steps; `fmax` is always the last value. From 50: 50, 51, 52, 54, 56, 60, 64, 72, …, fmax.
"""
function fragment_length_grid(fmin::Integer, fmax::Integer)
    1 <= fmin <= fmax || throw(ArgumentError("need 1 <= fmin <= fmax (got fmin=$fmin, fmax=$fmax)"))
    Fs = Int[fmin]
    step, n = 1, 0
    while (next = Fs[end] + step) < fmax
        push!(Fs, next)
        n += 1
        iseven(n) && (step *= 2)
    end
    Fs[end] == fmax || push!(Fs, fmax)
    return Fs
end

"""
    StepPolicy

Tiling step (bp between consecutive fragment starts) as a function of fragment length.
Rules are `(max_F, step)` pairs in ascending `max_F`; the last rule is a catch-all.
"""
struct StepPolicy
    rules::Vector{Tuple{Int,Int}}
end

const DEFAULT_STEP_POLICY = "150:1,300:2,*:5"

"""
    parse_step_policy("150:1,300:2,*:5") -> StepPolicy

`maxF:step` rules: F ≤ 150 → step 1, F ≤ 300 → step 2, otherwise step 5.
"""
function parse_step_policy(s::AbstractString)
    rules = Tuple{Int,Int}[]
    for tok in split(s, ',')
        parts = split(strip(tok), ':')
        length(parts) == 2 || throw(ArgumentError("bad step rule '$tok' (expected maxF:step)"))
        maxF = strip(parts[1]) == "*" ? typemax(Int) : parse(Int, parts[1])
        step = parse(Int, parts[2])
        step >= 1 || throw(ArgumentError("step must be >= 1 in rule '$tok'"))
        push!(rules, (maxF, step))
    end
    sort!(rules)
    rules[end][1] == typemax(Int) || throw(ArgumentError("step policy needs a catch-all rule, e.g. '*:5'"))
    return StepPolicy(rules)
end

tiling_step(p::StepPolicy, F::Integer) = p.rules[findfirst(r -> F <= r[1], p.rules)][2]

Base.string(p::StepPolicy) =
    join(("$(m == typemax(Int) ? "*" : m):$s" for (m, s) in p.rules), ",")

parse_int_list(s::AbstractString) = sort!(unique!([parse(Int, strip(x)) for x in split(s, ',') if !isempty(strip(x))]))
