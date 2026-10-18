const ST_MISSING = 0x00     # simulated but absent from the BAM
const ST_CORRECT = 0x01
const ST_WRONG = 0x02       # proper pair placed away from its origin
const ST_NOT_PROPER = 0x03  # unaligned, discordant, mixed-mode or incomplete pair
const ST_SKIPPED = 0x04     # N in sequenced bases; never simulated

const SITE_BIN = 1000
const DEFAULT_CUTOFFS = "0,2,5,10,20,30"
const FRAG_MAPQ_RULES = Dict("min" => min, "max" => max)

"""
    LocusFrags

Fate of every simulated fragment of one locus at one read length. Fragments are laid out
by F index, then start, then orientation (in `orients` order), matching `fragment_starts`.
"""
struct LocusFrags
    starts::Vector{StepRange{Int,Int}}
    offsets::Vector{Int}
    status::Vector{UInt8}
    mapq::Vector{UInt8}
    orients::Vector{Char}
end

function LocusFrags(locus::Locus, seq::AbstractVector{UInt8}, region_start0::Integer, chromlen::Integer,
                    L::Integer, Fs::AbstractVector{<:Integer}, step::Integer, orients::AbstractVector{Char})
    no = length(orients)
    starts = [fragment_starts(locus, F, step, chromlen) for F in Fs]
    offsets = [0; cumsum(no .* length.(starts))]
    status = zeros(UInt8, offsets[end])
    for (j, F) in enumerate(Fs), (i, f0) in enumerate(starts[j])
        if fragment_has_n(seq, f0 - region_start0 + 1, F, min(L, F))
            status[offsets[j]+no*(i-1)+1:offsets[j]+no*i] .= ST_SKIPPED
        end
    end
    return LocusFrags(starts, offsets, status, zeros(UInt8, offsets[end]), collect(orients))
end

function frag_slot(lf::LocusFrags, j::Integer, f0::Integer, orient::Char)
    r = lf.starts[j]
    o = findfirst(==(orient), lf.orients)
    (o === nothing || !(f0 in r)) && return 0
    return lf.offsets[j] + length(lf.orients) * div(f0 - first(r), step(r)) + o
end

struct Mate
    flag::UInt16
    refid::Int   # 1-based index into the BAM header; 0 if none
    pos0::Int
    mapq::UInt8
    tlen::Int
end

Mate(rec::BAM.Record) = Mate(BAM.flags(rec), BAM.refid(rec), BAM.position(rec) - 1,
                             BAM.mappingquality(rec), BAM.templength(rec))

"""
    classify_pair(m1, m2, origin, refnames; tolerance=0) -> status

`ST_CORRECT` if the pair is properly paired on the origin chromosome with leftmost
position and fragment length within `tolerance` of the origin and R1 on the expected
strand; `ST_WRONG` for other proper pairs; otherwise `ST_NOT_PROPER`.
"""
function classify_pair(m1::Union{Nothing,Mate}, m2::Union{Nothing,Mate}, origin, refnames; tolerance::Integer=0)
    (m1 === nothing || m2 === nothing) && return ST_NOT_PROPER
    proper = m1.flag & 0x2 != 0 && m1.flag & 0x4 == 0 && m2.flag & 0x4 == 0 && m1.refid == m2.refid > 0
    proper || return ST_NOT_PROPER
    left0 = min(m1.pos0, m2.pos0)
    r1_reverse = m1.flag & 0x10 != 0
    if refnames[m1.refid] == origin.chrom && abs(left0 - origin.frag_start0) <= tolerance &&
       abs(abs(m1.tlen) - origin.F) <= tolerance && r1_reverse == (origin.orient == '-')
        return ST_CORRECT
    end
    return ST_WRONG
end

const Intervals = Dict{String,Vector{Tuple{Int,Int,Int}}}
const NO_INTERVALS = Tuple{Int,Int,Int}[]

function interval_index(ls::AbstractVector{Locus})
    d = Intervals()
    for (k, l) in enumerate(ls)
        push!(get!(Vector{Tuple{Int,Int,Int}}, d, l.chrom), (l.start0, l.end0, k))
    end
    return d
end

"""Target locus with the largest overlap with `[s0, e0)` on `chrom`, or 0."""
function dest_locus(bychrom::Intervals, chrom::AbstractString, s0::Int, e0::Int)
    best, best_ov = 0, 0
    for (ls, le, k) in get(bychrom, chrom, NO_INTERVALS)
        ov = min(le, e0) - max(ls, s0)
        ov > best_ov && ((best, best_ov) = (k, ov))
    end
    return best
end

"""A misplaced proper pair landing on a simulation region: covers region-relative bases `a:b-1`."""
struct Landing
    a::Int32
    b::Int32
    j::Int32   # fragment-length index
    q::UInt8   # fragment MAPQ
end

mutable struct ScanStats
    records::Int
    pairs::Int
    incomplete::Int
    duplicate::Int
end

struct ScoreContext
    loci::Vector{Locus}
    regions::Vector{Locus}
    region_index::Dict{String,Int}
    Fidx::Dict{Int,Int}
    locus_bychrom::Intervals
    region_bychrom::Intervals
    tolerance::Int
    mapq_rule::Function
end

function record_pair!(frags, ctx::ScoreContext, refnames, name::AbstractString, m1, m2,
                      bleed, sites, landings, stats)
    o = parse_read_name(name)
    k = get(ctx.region_index, o.region_id, 0)
    k == 0 && error("read '$name' names region '$(o.region_id)', which is not in regions.tsv")
    j = get(ctx.Fidx, o.F, 0)
    s = j == 0 ? 0 : frag_slot(frags[k], j, o.frag_start0, o.orient)
    s == 0 && error("read '$name' does not match the simulation tiling (different simulate run?)")
    lf = frags[k]
    lf.status[s] == ST_SKIPPED && error("read '$name' should have been skipped (contains N)")
    if lf.status[s] != ST_MISSING
        stats.duplicate += 1
        return
    end
    stats.pairs += 1
    (m1 === nothing || m2 === nothing) && (stats.incomplete += 1)
    st = classify_pair(m1, m2, o, refnames; tolerance=ctx.tolerance)
    lf.status[s] = st
    st == ST_NOT_PROPER && return
    q = ctx.mapq_rule(m1.mapq, m2.mapq)
    lf.mapq[s] = q
    if st == ST_WRONG
        left0 = min(m1.pos0, m2.pos0)
        e0 = left0 + abs(m1.tlen)
        chrom = refnames[m1.refid]
        for (rs, re, r) in get(ctx.region_bychrom, chrom, NO_INTERVALS)
            (rs < e0 && left0 < re) || continue
            push!(landings[r], Landing(max(left0, rs) - rs + 1, min(e0, re) - rs + 1, j, q))
        end
        org = dest_locus(ctx.locus_bychrom, o.chrom, o.frag_start0, o.frag_start0 + o.F)
        org == 0 && return
        d = dest_locus(ctx.locus_bychrom, chrom, left0, e0)
        key = (org, d, o.F, Int(q))
        bleed[key] = get(bleed, key, 0) + 1
        skey = (org, m1.refid, div(left0, SITE_BIN), d, Int(q))
        sites[skey] = get(sites, skey, 0) + 1
    end
    return
end

"""
    scan_bam!(frags, ctx, bam) -> (stats, bleed, sites, landings, refnames)

Single pass over an unsorted (name-grouped) BAM. Secondary/supplementary records are ignored.
`landings[r]` lists misplaced pairs that landed on region `r`.
"""
function scan_bam!(frags::Vector{LocusFrags}, ctx::ScoreContext, bam::AbstractString)
    reader = BAM.Reader(open(bam))
    refnames = copy(reader.refseqnames)
    rec = BAM.Record()
    stats = ScanStats(0, 0, 0, 0)
    bleed = Dict{NTuple{4,Int},Int}()
    sites = Dict{NTuple{5,Int},Int}()
    landings = [Landing[] for _ in ctx.regions]
    cur = ""
    m1::Union{Nothing,Mate} = nothing
    m2::Union{Nothing,Mate} = nothing
    try
        while !eof(reader)
            read!(reader, rec)
            stats.records += 1
            fl = BAM.flags(rec)
            fl & 0x900 == 0 || continue
            name = BAM.tempname(rec)
            if name != cur
                isempty(cur) || record_pair!(frags, ctx, refnames, cur, m1, m2, bleed, sites, landings, stats)
                cur, m1, m2 = name, nothing, nothing
            end
            if fl & 0x40 != 0
                m1 = Mate(rec)
            else
                m2 = Mate(rec)
            end
        end
        isempty(cur) || record_pair!(frags, ctx, refnames, cur, m1, m2, bleed, sites, landings, stats)
    finally
        close(reader)
    end
    stats.duplicate > 0 && error("$bam: $(stats.duplicate) read pairs appear more than once or are not " *
                                 "adjacent; use the unsorted bowtie2 output (do not coordinate-sort it)")
    return stats, bleed, sites, landings, refnames
end

pushrow!(t::NamedTuple, row::NamedTuple) = (foreach(k -> push!(t[k], row[k]), keys(t)); t)
appendcols!(t::NamedTuple, u::NamedTuple) = (foreach(k -> append!(t[k], u[k]), keys(t)); t)

new_summary() = (locus_id=String[], chrom=String[], L=Int32[], F=Int32[], mapq_cutoff=Int32[],
                 n_fragments=Int64[], n_skipped_N=Int64[], n_missing=Int64[], n_not_proper=Int64[],
                 n_kept=Int64[], n_correct_kept=Int64[], n_wrong_kept=Int64[],
                 frac_wrong_kept=Float64[], frac_correct=Float64[], n_ectopic_landed=Int64[],
                 min_denominator=Int32[])
new_bleed() = (origin_locus=String[], dest_locus=String[], L=Int32[], F=Int32[], mapq_cutoff=Int32[],
               n_fragments=Int64[])
new_sites() = (origin_locus=String[], L=Int32[], dest_chrom=String[], dest_bin_start0=Int64[],
               dest_bin_end0=Int64[], dest_locus=String[], fragment_mapq=Int32[], n_fragments=Int64[])

"""
    region_tracks!(ns, nc, ne, region, lf, landings, Fs, cutoffs)

Per-base fragment-span counts for one simulation region (difference arrays):

- `ns[x, F]`: simulated fragments covering x
- `nc[x, F, c]`: of those, placed back correctly with fragment MAPQ ≥ c
- `ne[x, F, c]`: misplaced fragments from any simulated region that landed on x with MAPQ ≥ c
"""
function region_tracks!(ns::AbstractMatrix{UInt16}, nc::AbstractArray{UInt16,3}, ne::AbstractArray{UInt32,3},
                        region::Locus, lf::LocusFrags, landings::AbstractVector{Landing},
                        Fs::AbstractVector{<:Integer}, cutoffs::AbstractVector{<:Integer})
    n, ncut, no = length(region), length(cutoffs), length(lf.orients)
    den = zeros(Int32, n + 1)
    num = zeros(Int32, n + 1, ncut)
    foreign = zeros(Int32, n + 1, ncut)
    for (j, F) in enumerate(Fs)
        fill!(den, 0)
        fill!(num, 0)
        fill!(foreign, 0)
        for (i, f0) in enumerate(lf.starts[j]), o in 1:no
            s = lf.offsets[j] + no * (i - 1) + o
            st = lf.status[s]
            st == ST_SKIPPED && continue
            a = max(f0, region.start0) - region.start0 + 1
            b = min(f0 + F, region.end0) - region.start0 + 1
            den[a] += 1
            den[b] -= 1
            st == ST_CORRECT || continue
            for c in 1:ncut
                lf.mapq[s] >= cutoffs[c] || continue
                num[a, c] += 1
                num[b, c] -= 1
            end
        end
        for ev in landings
            ev.j == j || continue
            for c in 1:ncut
                ev.q >= cutoffs[c] || continue
                foreign[ev.a, c] += 1
                foreign[ev.b, c] -= 1
            end
        end
        ns[:, j] .= view(cumsum!(den, den), 1:n)
        for c in 1:ncut
            nc[:, j, c] .= view(cumsum!(view(num, :, c), view(num, :, c)), 1:n)
            ne[:, j, c] .= view(cumsum!(view(foreign, :, c), view(foreign, :, c)), 1:n)
        end
    end
    return ns, nc, ne
end

"""
    map_back(n_correct, n_simulated)

Pr(map-back) = correct / simulated: the fraction of fragments from a base that come back
to it with MAPQ ≥ c. NaN where nothing was simulated.
"""
map_back(nc, ns) = ns == 0 ? NaN32 : Float32(nc / ns)

"""
    ectopic_rate(n_ectopic, n_simulated)

Ex-map rate | other loci: misplaced fragments from the simulated panel landing on a base
per fragment simulated there. A rate, not a probability; can exceed 1.
"""
ectopic_rate(ne, ns) = ns == 0 ? NaN32 : Float32(ne / ns)

const LOG_RATIO_PSEUDOCOUNT = 0.5

"""
    log2_ratio(n_correct, n_ectopic; pseudocount=0.5)

log2((correct + p) / (ectopic + p)) = log2(map_back / ectopic_rate), the log2 odds that a
kept fragment at the base is native. 0 means as many ectopic as native fragments; NaN
when neither count has anything. With no ectopic fragments the value is a lower bound
that grows with the correct count.
"""
log2_ratio(nc, ne; pseudocount=LOG_RATIO_PSEUDOCOUNT) =
    nc == 0 && ne == 0 ? NaN32 : Float32(log2((nc + pseudocount) / (ne + pseudocount)))

"""Fate counts of fragments in `lf` with length index `j` that overlap `[s0, e0)`, per cutoff."""
function count_fragments(lf::LocusFrags, j::Integer, F::Integer, s0::Integer, e0::Integer,
                         cutoffs::AbstractVector{<:Integer})
    nc, no = length(cutoffs), length(lf.orients)
    nfrag = nskip = nmiss = nnp = 0
    ok, wrong = zeros(Int, nc), zeros(Int, nc)
    for (i, f0) in enumerate(lf.starts[j])
        (f0 < e0 && f0 + F > s0) || continue
        for o in 1:no
            s = lf.offsets[j] + no * (i - 1) + o
            st = lf.status[s]
            if st == ST_SKIPPED
                nskip += 1
                continue
            end
            nfrag += 1
            if st == ST_MISSING
                nmiss += 1
            elseif st == ST_NOT_PROPER
                nnp += 1
            else
                for c in 1:nc
                    lf.mapq[s] >= cutoffs[c] || continue
                    st == ST_CORRECT ? (ok[c] += 1) : (wrong[c] += 1)
                end
            end
        end
    end
    return (; nfrag, nskip, nmiss, nnp, ok, wrong)
end

function push_summary!(summary, id, chrom, L, F, cutoffs, n, ectopic, mind)
    for (c, cut) in enumerate(cutoffs)
        kept = n.ok[c] + n.wrong[c]
        pushrow!(summary, (locus_id=id, chrom=chrom, L=Int32(L), F=Int32(F), mapq_cutoff=Int32(cut),
                           n_fragments=n.nfrag, n_skipped_N=n.nskip, n_missing=n.nmiss, n_not_proper=n.nnp,
                           n_kept=kept, n_correct_kept=n.ok[c], n_wrong_kept=n.wrong[c],
                           frac_wrong_kept=kept == 0 ? NaN : n.wrong[c] / kept,
                           frac_correct=n.nfrag == 0 ? NaN : n.ok[c] / n.nfrag,
                           n_ectopic_landed=ectopic[c], min_denominator=Int32(mind)))
    end
end

"""Misplaced fragments of length index `j` landing on region rows `rows`, per cutoff."""
function count_landings(landings::AbstractVector{Landing}, j::Integer, rows::AbstractUnitRange, cutoffs)
    out = zeros(Int, length(cutoffs))
    for ev in landings
        (ev.j == j && ev.a <= last(rows) && ev.b > first(rows)) || continue
        for (c, cut) in enumerate(cutoffs)
            ev.q >= cut && (out[c] += 1)
        end
    end
    return out
end

"""Rows of `region`'s tracks that belong to `locus`."""
locus_rows(locus::Locus, region::Locus) = locus.start0-region.start0+1:locus.end0-region.start0

"""
    region_summaries!(summary, region, members, lf, landings, ns, L, Fs, cutoffs) -> totals

Per-locus summary rows (fragments overlapping each member locus) and per-(F, cutoff)
totals over the whole region for the genome-wide rows (each fragment counted once).
"""
function region_summaries!(summary, region::Locus, members::AbstractVector{Locus}, lf::LocusFrags,
                           landings::AbstractVector{Landing}, ns::AbstractMatrix{UInt16}, L, Fs, cutoffs)
    totals = Dict{Tuple{Int,Int},NTuple{8,Int}}()
    for (j, F) in enumerate(Fs)
        minds = Int[]
        for l in members
            rows = locus_rows(l, region)
            mind = Int(minimum(view(ns, rows, j)))
            push!(minds, mind)
            push_summary!(summary, l.id, l.chrom, L, F, cutoffs,
                          count_fragments(lf, j, F, l.start0, l.end0, cutoffs),
                          count_landings(landings, j, rows, cutoffs), mind)
        end
        n = count_fragments(lf, j, F, region.start0, region.end0, cutoffs)
        ect = count_landings(landings, j, 1:length(region), cutoffs)
        for (c, cut) in enumerate(cutoffs)
            totals[(F, cut)] = (n.nfrag, n.nskip, n.nmiss, n.nnp, n.ok[c], n.wrong[c], ect[c], minimum(minds))
        end
    end
    return totals
end

function genome_wide_rows!(summary, L, all_totals)
    for key in sort!(collect(keys(first(all_totals))))
        t = [tot[key] for tot in all_totals]
        nfrag, nskip, nmiss, nnp, ok, wrong, ect = (sum(x[i] for x in t) for i in 1:7)
        kept = ok + wrong
        pushrow!(summary, (locus_id="ALL", chrom="", L=Int32(L), F=Int32(key[1]), mapq_cutoff=Int32(key[2]),
                           n_fragments=nfrag, n_skipped_N=nskip, n_missing=nmiss, n_not_proper=nnp,
                           n_kept=kept, n_correct_kept=ok, n_wrong_kept=wrong,
                           frac_wrong_kept=kept == 0 ? NaN : wrong / kept,
                           frac_correct=nfrag == 0 ? NaN : ok / nfrag, n_ectopic_landed=ect,
                           min_denominator=Int32(minimum(x[8] for x in t))))
    end
end

function bleed_rows!(bleed_tbl, sites_tbl, loci, L, cutoffs, bleed, sites, refnames)
    lid(d) = d == 0 ? "elsewhere" : loci[d].id
    agg = Dict{NTuple{4,Int},Int}()
    for ((k, d, F, q), n) in bleed, c in cutoffs
        q >= c && (agg[(k, d, F, c)] = get(agg, (k, d, F, c), 0) + n)
    end
    for key in sort!(collect(keys(agg)))
        k, d, F, c = key
        pushrow!(bleed_tbl, (origin_locus=loci[k].id, dest_locus=lid(d), L=Int32(L), F=Int32(F),
                             mapq_cutoff=Int32(c), n_fragments=agg[key]))
    end
    for key in sort!(collect(keys(sites)))
        k, r, bin, d, q = key
        pushrow!(sites_tbl, (origin_locus=loci[k].id, L=Int32(L), dest_chrom=refnames[r],
                             dest_bin_start0=bin * SITE_BIN, dest_bin_end0=(bin + 1) * SITE_BIN,
                             dest_locus=lid(d), fragment_mapq=Int32(q), n_fragments=sites[key]))
    end
end

function read_tiling(path::AbstractString)
    d = Dict{Tuple{String,Int,Int},Tuple{Int,Int}}()
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue
        f = split(line, '\t')
        d[(String(f[1]), parse(Int, f[2]), parse(Int, f[3]))] = (parse(Int, f[7]), parse(Int, f[8]))
    end
    return d
end

function write_table(path::AbstractString, tbl::NamedTuple)
    isfile(path) && rm(path)
    # Parquet.jl cannot write zero-row tables
    if isempty(first(tbl))
        @info "No rows for $(basename(path)); file not written"
        return
    end
    write_parquet(path, tbl)
end

"""
    score_simulation(simdir; bam_dir, out, cutoffs, tolerance=0, frag_mapq="min", fasta=nothing, read_lengths=nothing)

Phase 2 scoring. Reads `bam_dir/sim_L{L}.bam` for each simulated read length and writes
`mappability.zarr` (per-base tracks, see `create_track_store`), `misplacement_summary.parquet`,
`bleed_matrix.parquet`, `misplaced_sites.parquet` and `score.toml` to `out`.
"""
function score_simulation(simdir::AbstractString; bam_dir::AbstractString=joinpath(simdir, "aln"),
               out::AbstractString=joinpath(simdir, "scores"),
               cutoffs::AbstractVector{<:Integer}=parse_int_list(DEFAULT_CUTOFFS),
               tolerance::Integer=0, frag_mapq::AbstractString="min",
               fasta::Union{Nothing,AbstractString}=nothing,
               read_lengths::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    meta = read_simulate_toml(simdir)
    grid = meta["grid"]
    Fs = Vector{Int}(grid["fragment_lengths"])
    tile_step = Int(grid["step"])
    orients = parse_orientation(grid["orientation"])
    Ls = read_lengths === nothing ? Vector{Int}(grid["read_lengths"]) : collect(Int, read_lengths)
    issubset(Ls, grid["read_lengths"]) || error("read lengths $Ls were not all simulated ($(grid["read_lengths"]))")
    rule = get(FRAG_MAPQ_RULES, frag_mapq, nothing)
    rule === nothing && error("--frag-mapq must be one of: " * join(sort!(collect(keys(FRAG_MAPQ_RULES))), ", "))
    tolerance >= 0 || error("tolerance must be >= 0")
    bams = Dict(L => joinpath(bam_dir, "sim_L$L.bam") for L in Ls)
    for (L, b) in bams
        isfile(b) || error("BAM for L=$L not found: $b")
    end

    fa = IndexedFasta(something(fasta, meta["reference"]["fasta"]))
    loci = read_loci(joinpath(simdir, "loci.tsv"))
    regions_tsv = joinpath(simdir, "regions.tsv")
    isfile(regions_tsv) || error("$regions_tsv not found; re-run `unimap simulate` with this version")
    regions = read_loci(regions_tsv)
    locus_index = Dict(l.id => k for (k, l) in enumerate(loci))
    members = [[loci[locus_index[id]] for id in region_members(r)] for r in regions]
    tiling = read_tiling(joinpath(simdir, "tiling.tsv"))
    checksum = meta["reference"]["sha256"]
    ctx = ScoreContext(loci, regions, Dict(r.id => k for (k, r) in enumerate(regions)),
                       Dict(F => j for (j, F) in enumerate(Fs)), interval_index(loci), interval_index(regions),
                       tolerance, rule)

    fmax = maximum(Fs)
    2fmax <= typemax(UInt16) || error("fragment length $fmax is too large for UInt16 denominators")
    padded = [sim_region(r, fmax, chrom_length(fa, r.chrom)) for r in regions]
    seqs = [fetch_seq(fa, r.chrom, p...) for (r, p) in zip(regions, padded)]

    mkpath(out)
    zpath = joinpath(out, "mappability.zarr")
    store = create_track_store(zpath, loci, Ls, Fs, cutoffs;
                               attrs=Dict("genome_checksum" => checksum, "fragment_mapq_rule" => frag_mapq,
                                          "tolerance" => tolerance, "step" => tile_step,
                                          "orientation" => grid["orientation"],
                                          "precision_sources" => [r.id for r in regions]))
    zlock = ReentrantLock()

    results = Vector{Any}(undef, length(Ls))
    @info "Scoring $(length(loci)) loci in $(length(regions)) regions" read_lengths = join(Ls, ",") cutoffs = join(cutoffs, ",") threads = Threads.nthreads()
    Threads.@threads :dynamic for t in eachindex(Ls)
        L = Ls[t]
        frags = [LocusFrags(r, seqs[k], padded[k][1], chrom_length(fa, r.chrom), L, Fs, tile_step, orients)
                 for (k, r) in enumerate(regions)]
        for (k, r) in enumerate(regions), (j, F) in enumerate(Fs)
            exp = get(tiling, (r.id, L, F), nothing)
            exp === nothing && error("tiling.tsv has no row for ($(r.id), L=$L, F=$F)")
            lf = frags[k]
            got = (length(lf.starts[j]),
                   count(==(ST_SKIPPED), view(lf.status, lf.offsets[j]+1:lf.offsets[j+1])) ÷ length(orients))
            exp == got || error("tiling mismatch for ($(r.id), L=$L, F=$F): simulate wrote (starts, N-skipped) = " *
                                "$exp, rescoring expects $got; was the reference changed?")
        end
        stats, bleed, sites, landings, refnames = scan_bam!(frags, ctx, bams[L])
        summary, bleed_tbl, sites_tbl = new_summary(), new_bleed(), new_sites()
        totals = map(enumerate(regions)) do (k, r)
            ns = Matrix{UInt16}(undef, length(r), length(Fs))
            ncor = Array{UInt16}(undef, length(r), length(Fs), length(cutoffs))
            nect = Array{UInt32}(undef, length(r), length(Fs), length(cutoffs))
            region_tracks!(ns, ncor, nect, r, frags[k], landings[k], Fs, cutoffs)
            tot = region_summaries!(summary, r, members[k], frags[k], landings[k], ns, L, Fs, cutoffs)
            lock(zlock) do
                for l in members[k]
                    i = locus_index[l.id]
                    rows, rr = store.offsets[i]+1:store.offsets[i+1], locus_rows(l, r)
                    store.n_simulated[rows, :, t] = ns[rr, :]
                    store.n_correct[rows, :, :, t] = ncor[rr, :, :]
                    store.n_ectopic[rows, :, :, t] = nect[rr, :, :]
                end
            end
            tot
        end
        genome_wide_rows!(summary, L, totals)
        bleed_rows!(bleed_tbl, sites_tbl, loci, L, cutoffs, bleed, sites, refnames)
        results[t] = (; stats, summary, bleed_tbl, sites_tbl)
        @info "Scored L=$L" bam = bams[L] pairs = stats.pairs incomplete = stats.incomplete
    end
    Zarr.consolidate_metadata(store.group)

    all_summary, all_bleed, all_sites = new_summary(), new_bleed(), new_sites()
    for r in results
        appendcols!(all_summary, r.summary)
        appendcols!(all_bleed, r.bleed_tbl)
        appendcols!(all_sites, r.sites_tbl)
    end
    write_table(joinpath(out, "misplacement_summary.parquet"), all_summary)
    write_table(joinpath(out, "bleed_matrix.parquet"), all_bleed)
    write_table(joinpath(out, "misplaced_sites.parquet"), all_sites)

    align_toml = joinpath(simdir, "align.toml")
    version_txt = joinpath(bam_dir, "bowtie2_version.txt")
    smeta = Dict{String,Any}(
        "unimap" => Dict("version" => string(pkgversion(@__MODULE__)), "git_commit" => git_commit(),
                         "created" => string(now()), "julia" => string(VERSION)),
        "reference" => Dict("fasta" => abspath(fa.path), "sha256" => checksum),
        "simulation" => Dict("dir" => abspath(simdir), "grid" => grid),
        "alignment" => Dict("settings" => isfile(align_toml) ? TOML.parsefile(align_toml) : "unknown",
                            "bowtie2_version" => isfile(version_txt) ? read(version_txt, String) : "unknown",
                            "bams" => Dict(string(L) => abspath(b) for (L, b) in bams)),
        "scoring" => Dict("mapq_cutoffs" => collect(Int, cutoffs), "tolerance" => tolerance,
                          "fragment_mapq_rule" => frag_mapq, "coverage" => "fragment span",
                          "site_bin" => SITE_BIN),
        "bam_stats" => Dict(string(Ls[t]) => Dict("records" => r.stats.records, "pairs" => r.stats.pairs,
                                                   "incomplete_pairs" => r.stats.incomplete)
                            for (t, r) in enumerate(results)),
    )
    open(io -> TOML.print(io, smeta), joinpath(out, "score.toml"), "w")
    @info "Wrote scores" out
    return (; tracks=zpath, summary=all_summary, bleed=all_bleed, sites=all_sites)
end

const CHUNK_POSITIONS = 16_384

"""
    create_track_store(path, loci, Ls, Fs, cutoffs; attrs) -> (; group, n_simulated, n_correct, n_ectopic, offsets)

Zarr (v2) group holding per-base fragment counts for all loci concatenated end to end:

- `n_simulated[position, F, L]` (UInt16): simulated fragments covering the base
- `n_correct[position, F, mapq_cutoff, L]` (UInt16): of those, placed back correctly with MAPQ ≥ cutoff
- `n_ectopic[position, F, mapq_cutoff, L]` (UInt32): misplaced fragments from the simulated
  panel that landed on the base with MAPQ ≥ cutoff
- 1-D coordinate arrays `L`, `F`, `mapq_cutoff`

Derived quantities: `map_back`, `ectopic_rate`, `log2_ratio`. Locus `k` occupies positions
`offsets[k]+1:offsets[k+1]` (attribute `loci`). Dimension names are stored as
`_ARRAY_DIMENSIONS` in C order (what xarray/zarr-python see).
"""
function create_track_store(path::AbstractString, loci::AbstractVector{Locus}, Ls, Fs, cutoffs;
                            attrs::AbstractDict=Dict{String,Any}())
    isdir(path) && rm(path; recursive=true)
    offsets = [0; cumsum(length.(loci))]
    npos, nF, nc, nL = offsets[end], length(Fs), length(cutoffs), length(Ls)
    chunk = min(CHUNK_POSITIONS, npos)
    comp = Zarr.BloscCompressor(cname="zstd", clevel=5, shuffle=2)
    loci_attr = [Dict("locus_id" => l.id, "chrom" => l.chrom, "start0" => l.start0, "end0" => l.end0,
                      "strand" => string(l.strand), "offset" => offsets[k]) for (k, l) in enumerate(loci)]
    gattrs = merge(Dict{String,Any}(attrs), Dict("loci" => loci_attr, "read_lengths" => collect(Int, Ls),
                   "fragment_lengths" => collect(Int, Fs), "mapq_cutoffs" => collect(Int, cutoffs),
                   "n_positions" => npos, "log_ratio_pseudocount" => LOG_RATIO_PSEUDOCOUNT))
    g = zgroup(path; attrs=gattrs)
    dims4 = Dict("_ARRAY_DIMENSIONS" => ["L", "mapq_cutoff", "F", "position"])
    n_simulated = zcreate(UInt16, g, "n_simulated", npos, nF, nL; chunks=(chunk, nF, 1), fill_value=UInt16(0),
                          compressor=comp, attrs=Dict("_ARRAY_DIMENSIONS" => ["L", "F", "position"]))
    n_correct = zcreate(UInt16, g, "n_correct", npos, nF, nc, nL; chunks=(chunk, nF, nc, 1),
                        fill_value=UInt16(0), compressor=comp, attrs=dims4)
    n_ectopic = zcreate(UInt32, g, "n_ectopic", npos, nF, nc, nL; chunks=(chunk, nF, nc, 1),
                        fill_value=UInt32(0), compressor=comp, attrs=dims4)
    for (name, vals) in (("L", Ls), ("F", Fs), ("mapq_cutoff", cutoffs))
        z = zcreate(Int32, g, name, length(vals); attrs=Dict("_ARRAY_DIMENSIONS" => [name]))
        z[:] = Int32.(vals)
    end
    return (; group=g, n_simulated, n_correct, n_ectopic, offsets)
end

"""
    locus_track(zarr, locus_id; L, F, mapq_cutoff)
        -> (; chrom, start0, n_simulated, n_correct, n_ectopic, map_back, ectopic_rate, log2_ratio)

Per-base counts and derived quantities for one locus and channel from a `mappability.zarr`
(path or opened group). Position `i` of the vectors is genomic `start0 + i - 1` (0-based).
"""
function locus_track(g::Zarr.ZGroup, locus_id::AbstractString; L::Integer, F::Integer, mapq_cutoff::Integer)
    haskey(g.arrays, "n_correct") || error("store has no count arrays; re-run `unimap score` with this version")
    loci = g.attrs["loci"]
    k = findfirst(l -> l["locus_id"] == locus_id, loci)
    k === nothing && error("locus '$locus_id' is not in this store")
    idx(key, v) = something(findfirst(==(v), g.attrs[key]), 0)
    iL, iF, ic = idx("read_lengths", L), idx("fragment_lengths", F), idx("mapq_cutoffs", mapq_cutoff)
    (iL == 0 || iF == 0 || ic == 0) && error("(L=$L, F=$F, mapq_cutoff=$mapq_cutoff) is not on the stored grid")
    l = loci[k]
    rows = l["offset"]+1:l["offset"]+l["end0"]-l["start0"]
    ns = g["n_simulated"][rows, iF, iL]
    nc = g["n_correct"][rows, iF, ic, iL]
    ne = g["n_ectopic"][rows, iF, ic, iL]
    pc = get(g.attrs, "log_ratio_pseudocount", LOG_RATIO_PSEUDOCOUNT)
    return (; chrom=l["chrom"], start0=l["start0"], n_simulated=ns, n_correct=nc, n_ectopic=ne,
            map_back=map_back.(nc, ns), ectopic_rate=ectopic_rate.(ne, ns),
            log2_ratio=log2_ratio.(nc, ne; pseudocount=pc))
end

locus_track(path::AbstractString, locus_id::AbstractString; kw...) = locus_track(zopen(path), locus_id; kw...)
