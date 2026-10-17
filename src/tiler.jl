const COMPLEMENT = let t = fill(UInt8('N'), 256)
    for (a, b) in zip("ACGTN", "TGCAN")
        t[UInt8(a)+1] = UInt8(b)
    end
    t
end

isacgt(b::UInt8) = b == UInt8('A') || b == UInt8('C') || b == UInt8('G') || b == UInt8('T')

"""
    fragment_starts(locus, F, step, chromlen) -> StepRange

0-based starts of every simulated fragment of length `F` that overlaps the scored
region of `locus`, every `step` bp. Equivalent to tiling the locus padded by `F - 1` on each side.
"""
fragment_starts(locus::Locus, F::Integer, step::Integer, chromlen::Integer) =
    max(0, locus.start0 - F + 1):step:min(locus.end0 - 1, chromlen - F)

"""0-based half-open interval of reference sequence needed to tile `locus` with fragments up to `fmax`."""
sim_region(locus::Locus, fmax::Integer, chromlen::Integer) =
    (max(0, locus.start0 - fmax + 1), min(chromlen, locus.end0 + fmax - 1))

"""
    fragment_has_n(seq, i, F, m)

True if either sequenced end (first `m` or last `m` bases) of the fragment starting at
1-based index `i` of `seq` contains a non-ACGT base. Such fragments are not simulated
and are excluded from the denominator.
"""
function fragment_has_n(seq::AbstractVector{UInt8}, i::Integer, F::Integer, m::Integer)
    @inbounds for k in i:i+m-1
        isacgt(seq[k]) || return true
    end
    @inbounds for k in i+F-m:i+F-1
        isacgt(seq[k]) || return true
    end
    return false
end

"""
    read_name(locus_id, chrom, frag_start0, F, orient)

`locusID|chrom|fragStart0|F|orient`. `orient` is `+` (R1 from the forward strand at the
fragment start) or `-` (R1 = reverse complement of the fragment end).
"""
read_name(locus_id, chrom, frag_start0, F, orient::Char) = "$locus_id|$chrom|$frag_start0|$F|$orient"

function parse_read_name(name::AbstractString)
    f = split(name, '|')
    length(f) == 5 || throw(ArgumentError("not a UniMap read name: $name"))
    return (locus_id=String(f[1]), chrom=String(f[2]), frag_start0=parse(Int, f[3]),
            F=parse(Int, f[4]), orient=f[5][1])
end

const MAX_QNAME = 254

function write_fastq_record(io::IO, name::AbstractString, seq::AbstractVector{UInt8}, qual::AbstractVector{UInt8})
    write(io, UInt8('@'), name, UInt8('\n'))
    write(io, seq)
    write(io, "\n+\n")
    write(io, qual)
    write(io, UInt8('\n'))
end

"""
    TileSummary

One row of `tiling.tsv`: what was simulated for one (locus, L, F). Scoring re-enumerates
fragments with `fragment_starts`/`fragment_has_n` to rebuild denominators.
"""
struct TileSummary
    locus_id::String
    L::Int
    F::Int
    step::Int
    first_start0::Int
    last_start0::Int
    n_starts::Int
    n_skipped_N::Int
    n_pairs::Int
end

"""
    tile_locus!(io1, io2, locus, seq, region_start0, chromlen, L, Fs, step; qual='I', orientations=['+', '-'])

Write the requested orientations of every fragment (all `Fs`, starts every `step` bp) for
one locus and read length `L`. Mates are the first/last `min(L, F)` bases of the fragment
(F < L mimics adapter trimming). `seq` holds the reference from `region_start0` (see `sim_region`).
"""
function tile_locus!(io1::IO, io2::IO, locus::Locus, seq::AbstractVector{UInt8}, region_start0::Integer,
                     chromlen::Integer, L::Integer, Fs::AbstractVector{<:Integer}, step::Integer;
                     qual::Char='I', orientations::AbstractVector{Char}=['+', '-'])
    fwd = Vector{UInt8}(undef, L)
    rev = Vector{UInt8}(undef, L)
    q = fill(UInt8(qual), L)
    rows = TileSummary[]
    for F in Fs
        m = min(L, F)
        starts = fragment_starts(locus, F, step, chromlen)
        fv, rv, qv = view(fwd, 1:m), view(rev, 1:m), view(q, 1:m)
        nskip = 0
        for f0 in starts
            i = f0 - region_start0 + 1
            if fragment_has_n(seq, i, F, m)
                nskip += 1
                continue
            end
            copyto!(fwd, 1, seq, i, m)
            @inbounds for k in 1:m
                rev[k] = COMPLEMENT[seq[i+F-k]+1]
            end
            for o in orientations
                name = read_name(locus.id, locus.chrom, f0, F, o)
                write_fastq_record(io1, name, o == '+' ? fv : rv, qv)
                write_fastq_record(io2, name, o == '+' ? rv : fv, qv)
            end
        end
        n = length(starts)
        push!(rows, TileSummary(locus.id, L, F, step, n > 0 ? first(starts) : -1, n > 0 ? last(starts) : -1,
                                n, nskip, length(orientations) * (n - nskip)))
    end
    return rows
end

fastq_paths(outdir::AbstractString, locus_id::AbstractString, L::Integer) =
    (joinpath(outdir, "fastq", "L$L", "$(locus_id)_R1.fq.gz"),
     joinpath(outdir, "fastq", "L$L", "$(locus_id)_R2.fq.gz"))

function with_gzip(f, path::AbstractString; level::Integer=1)
    io = GzipCompressorStream(open(path, "w"); level=level)
    try
        f(io)
    finally
        close(io)
    end
end

function validate_loci(fa::IndexedFasta, loci::AbstractVector{Locus}, fmax::Integer)
    ids = [l.id for l in loci]
    allunique(ids) || error("duplicate locus IDs: " * join(unique(filter(i -> count(==(i), ids) > 1, ids)), ", "))
    for l in loci
        occursin(r"[\s|/]", l.id) && error("locus ID '$(l.id)' contains whitespace, '|' or '/'")
        haskey(fa, l.chrom) || error("locus $(l.id): '$(l.chrom)' not in reference")
        0 <= l.start0 < l.end0 <= chrom_length(fa, l.chrom) || error("locus $(l.id): bad interval [$(l.start0), $(l.end0))")
        n = length(read_name(l.id, l.chrom, chrom_length(fa, l.chrom), fmax, '+'))
        n <= MAX_QNAME || error("locus $(l.id): read names would be $n characters (max $MAX_QNAME)")
    end
end

function git_commit()
    dir = pkgdir(@__MODULE__)
    try
        rev = readchomp(pipeline(`git -C $dir rev-parse HEAD`; stderr=devnull))
        dirty = !isempty(readchomp(pipeline(`git -C $dir status --porcelain --untracked-files=no`; stderr=devnull)))
        return dirty ? rev * "-dirty" : rev
    catch
        return "unknown"
    end
end

function write_tiling(path::AbstractString, rows::AbstractVector{TileSummary})
    open(path, "w") do io
        println(io, join(fieldnames(TileSummary), '\t'))
        for r in rows
            println(io, join((getfield(r, f) for f in fieldnames(TileSummary)), '\t'))
        end
    end
end

"""
    simulate(fa, loci, outdir; read_lengths, fragment_lengths, step=1, orientation="both", qual='I', checksum=true, provenance=Dict())

Phase 1 tiler. For each (locus, L) writes `fastq/L{L}/{locus}_R1.fq.gz` and `_R2.fq.gz`
containing every F in `fragment_lengths`, with fragment starts every `step` bp, plus
`loci.tsv`, `loci.bed`, `tiling.tsv` and `simulate.toml` (provenance) in `outdir`.
`orientation` is `both`, `plus` or `minus`. Threads over loci.
"""
function simulate(fa::IndexedFasta, loci::AbstractVector{Locus}, outdir::AbstractString;
                  read_lengths::AbstractVector{<:Integer}, fragment_lengths::AbstractVector{<:Integer},
                  step::Integer=1, orientation::AbstractString="both",
                  qual::Char='I', checksum::Bool=true, provenance::AbstractDict=Dict{String,Any}())
    isempty(loci) && error("no loci to simulate")
    Ls = check_lengths("read lengths", read_lengths)
    Fs = check_lengths("fragment lengths", fragment_lengths)
    step >= 1 || throw(ArgumentError("step must be >= 1"))
    orients = parse_orientation(orientation)
    fmax = maximum(Fs)
    validate_loci(fa, loci, fmax)

    for L in Ls
        mkpath(joinpath(outdir, "fastq", "L$L"))
    end
    sha = checksum ? (@info("Computing SHA-256 of $(fa.path)"); file_sha256(fa.path)) : "not computed"

    @info "Tiling $(length(loci)) loci" read_lengths = join(Ls, ",") fragment_lengths = join(Fs, ",") step orientation threads = Threads.nthreads()
    results = Vector{Vector{TileSummary}}(undef, length(loci))
    Threads.@threads :dynamic for k in eachindex(loci)
        locus = loci[k]
        chromlen = chrom_length(fa, locus.chrom)
        rs, re = sim_region(locus, fmax, chromlen)
        seq = fetch_seq(fa, locus.chrom, rs, re)
        rows = TileSummary[]
        for L in Ls
            p1, p2 = fastq_paths(outdir, locus.id, L)
            with_gzip(p1) do io1
                with_gzip(p2) do io2
                    append!(rows, tile_locus!(io1, io2, locus, seq, rs, chromlen, L, Fs, step;
                                              qual, orientations=orients))
                end
            end
        end
        results[k] = rows
        @info "Tiled $(locus.id)" pairs = sum(r -> r.n_pairs, rows) skipped_N = sum(r -> r.n_skipped_N, rows)
    end
    # not named `rows`: that would make the threaded closure's `rows` a shared variable
    all_rows = reduce(vcat, results)

    write_loci(joinpath(outdir, "loci.tsv"), loci)
    write_loci_bed(joinpath(outdir, "loci.bed"), loci)
    write_tiling(joinpath(outdir, "tiling.tsv"), all_rows)
    meta = Dict{String,Any}(
        "unimap" => Dict("version" => string(pkgversion(@__MODULE__)), "git_commit" => git_commit(),
                         "created" => string(now()), "julia" => string(VERSION)),
        "reference" => Dict("fasta" => abspath(fa.path), "sha256" => sha),
        "grid" => Dict("read_lengths" => Ls, "fragment_lengths" => Fs, "step" => step,
                       "orientation" => orientation, "quality_char" => string(qual)),
        "output" => Dict("n_loci" => length(loci), "n_pairs_per_L" => Dict(string(L) => sum(r.n_pairs for r in all_rows if r.L == L) for L in Ls),
                         "read_name_format" => "locusID|chrom|fragStart0|F|orient"),
        "provenance" => provenance,
    )
    open(io -> TOML.print(io, meta), joinpath(outdir, "simulate.toml"), "w")
    @info "Done" outdir total_pairs = sum(r -> r.n_pairs, all_rows)
    return all_rows
end
