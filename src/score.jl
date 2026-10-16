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
by F index, then start, then orientation (`+`, `-`), matching `fragment_starts`.
"""
struct LocusFrags
    starts::Vector{StepRange{Int,Int}}
    offsets::Vector{Int}
    status::Vector{UInt8}
    mapq::Vector{UInt8}
end

function LocusFrags(locus::Locus, seq::AbstractVector{UInt8}, region_start0::Integer, chromlen::Integer,
                    L::Integer, Fs::AbstractVector{<:Integer}, policy::StepPolicy)
    starts = [fragment_starts(locus, F, tiling_step(policy, F), chromlen) for F in Fs]
    offsets = [0; cumsum(2 .* length.(starts))]
    status = zeros(UInt8, offsets[end])
    for (j, F) in enumerate(Fs), (i, f0) in enumerate(starts[j])
        if fragment_has_n(seq, f0 - region_start0 + 1, F, min(L, F))
            status[offsets[j]+2i-1] = status[offsets[j]+2i] = ST_SKIPPED
        end
    end
    return LocusFrags(starts, offsets, status, zeros(UInt8, offsets[end]))
end

function frag_slot(lf::LocusFrags, j::Integer, f0::Integer, orient::Char)
    r = lf.starts[j]
    f0 in r || return 0
    return lf.offsets[j] + 2 * div(f0 - first(r), step(r)) + (orient == '-' ? 2 : 1)
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

"""Target locus with the largest overlap with `[s0, e0)` on `chrom`, or 0."""
function dest_locus(bychrom::Dict{String,Vector{Tuple{Int,Int,Int}}}, chrom::AbstractString, s0::Int, e0::Int)
    best, best_ov = 0, 0
    for (ls, le, k) in get(bychrom, chrom, Tuple{Int,Int,Int}[])
        ov = min(le, e0) - max(ls, s0)
        ov > best_ov && ((best, best_ov) = (k, ov))
    end
    return best
end

mutable struct ScanStats
    records::Int
    pairs::Int
    incomplete::Int
    duplicate::Int
end

struct ScoreContext
    loci::Vector{Locus}
    locus_index::Dict{String,Int}
    Fidx::Dict{Int,Int}
    bychrom::Dict{String,Vector{Tuple{Int,Int,Int}}}
    tolerance::Int
    mapq_rule::Function
end

function record_pair!(frags, ctx::ScoreContext, refnames, name::AbstractString, m1, m2, bleed, sites, stats)
    o = parse_read_name(name)
    k = get(ctx.locus_index, o.locus_id, 0)
    k == 0 && error("read '$name' names locus '$(o.locus_id)', which is not in loci.tsv")
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
        chrom = refnames[m1.refid]
        d = dest_locus(ctx.bychrom, chrom, left0, left0 + abs(m1.tlen))
        key = (k, d, o.F, Int(q))
        bleed[key] = get(bleed, key, 0) + 1
        skey = (k, m1.refid, div(left0, SITE_BIN), d, Int(q))
        sites[skey] = get(sites, skey, 0) + 1
    end
    return
end

"""
    scan_bam!(frags, ctx, bam) -> (stats, bleed, sites, refnames)

Single pass over an unsorted (name-grouped) BAM. Secondary/supplementary records are ignored.
"""
function scan_bam!(frags::Vector{LocusFrags}, ctx::ScoreContext, bam::AbstractString)
    reader = BAM.Reader(open(bam))
    refnames = copy(reader.refseqnames)
    rec = BAM.Record()
    stats = ScanStats(0, 0, 0, 0)
    bleed = Dict{NTuple{4,Int},Int}()
    sites = Dict{NTuple{5,Int},Int}()
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
                isempty(cur) || record_pair!(frags, ctx, refnames, cur, m1, m2, bleed, sites, stats)
                cur, m1, m2 = name, nothing, nothing
            end
            if fl & 0x40 != 0
                m1 = Mate(rec)
            else
                m2 = Mate(rec)
            end
        end
        isempty(cur) || record_pair!(frags, ctx, refnames, cur, m1, m2, bleed, sites, stats)
    finally
        close(reader)
    end
    stats.duplicate > 0 && error("$bam: $(stats.duplicate) read pairs appear more than once or are not " *
                                 "adjacent; use the unsorted bowtie2 output (do not coordinate-sort it)")
    return stats, bleed, sites, refnames
end

pushrow!(t::NamedTuple, row::NamedTuple) = (foreach(k -> push!(t[k], row[k]), keys(t)); t)
appendcols!(t::NamedTuple, u::NamedTuple) = (foreach(k -> append!(t[k], u[k]), keys(t)); t)

new_tracks() = (genome_checksum=String[], locus_id=String[], chrom=String[], start0=Int32[], end0=Int32[],
                L=Int32[], F=Int32[], mapq_cutoff=Int32[], score=Float64[], min_denominator=Int32[])
new_summary() = (locus_id=String[], chrom=String[], L=Int32[], F=Int32[], mapq_cutoff=Int32[],
                 n_fragments=Int64[], n_skipped_N=Int64[], n_missing=Int64[], n_not_proper=Int64[],
                 n_kept=Int64[], n_correct_kept=Int64[], n_wrong_kept=Int64[],
                 frac_wrong_kept=Float64[], frac_correct=Float64[], min_denominator=Int32[])
new_bleed() = (origin_locus=String[], dest_locus=String[], L=Int32[], F=Int32[], mapq_cutoff=Int32[],
               n_fragments=Int64[])
new_sites() = (origin_locus=String[], L=Int32[], dest_chrom=String[], dest_bin_start0=Int64[],
               dest_bin_end0=Int64[], dest_locus=String[], fragment_mapq=Int32[], n_fragments=Int64[])

"""
    locus_scores!(tracks, summary, locus, lf, L, Fs, cutoffs, checksum)

Fragment-span coverage via difference arrays: the denominator counts every simulated
fragment covering a position; the numerator for cutoff c counts those that mapped
correctly with fragment MAPQ ≥ c. Writes run-length encoded tracks and per-F summaries.
Returns per-(F, cutoff) totals for the genome-wide summary.
"""
function locus_scores!(tracks, summary, locus::Locus, lf::LocusFrags, L::Integer,
                       Fs::AbstractVector{<:Integer}, cutoffs::AbstractVector{<:Integer}, checksum::AbstractString)
    n, nc = length(locus), length(cutoffs)
    den = zeros(Int32, n + 1)
    num = zeros(Int32, n + 1, nc)
    totals = Dict{Tuple{Int,Int},NTuple{7,Int}}()
    for (j, F) in enumerate(Fs)
        fill!(den, 0)
        fill!(num, 0)
        nfrag = nskip = nmiss = nnp = 0
        kept_ok, kept_wrong = zeros(Int, nc), zeros(Int, nc)
        for (i, f0) in enumerate(lf.starts[j]), o in 1:2
            s = lf.offsets[j] + 2(i - 1) + o
            st = lf.status[s]
            if st == ST_SKIPPED
                nskip += 1
                continue
            end
            nfrag += 1
            a = max(f0, locus.start0) - locus.start0 + 1
            b = min(f0 + F, locus.end0) - locus.start0 + 1
            den[a] += 1
            den[b] -= 1
            if st == ST_MISSING
                nmiss += 1
            elseif st == ST_NOT_PROPER
                nnp += 1
            else
                q = lf.mapq[s]
                for c in 1:nc
                    q >= cutoffs[c] || continue
                    if st == ST_CORRECT
                        kept_ok[c] += 1
                        num[a, c] += 1
                        num[b, c] -= 1
                    else
                        kept_wrong[c] += 1
                    end
                end
            end
        end
        cumsum!(den, den)
        mind = n > 0 ? Int(minimum(view(den, 1:n))) : 0
        for c in 1:nc
            col = view(num, :, c)
            cumsum!(col, col)
            ratio(x) = den[x] == 0 ? NaN : col[x] / den[x]
            x = 1
            while x <= n
                sc, mn, y = ratio(x), den[x], x
                while y < n && isequal(ratio(y + 1), sc)
                    y += 1
                    mn = min(mn, den[y])
                end
                pushrow!(tracks, (genome_checksum=checksum, locus_id=locus.id, chrom=locus.chrom,
                                  start0=Int32(locus.start0 + x - 1), end0=Int32(locus.start0 + y),
                                  L=Int32(L), F=Int32(F), mapq_cutoff=Int32(cutoffs[c]), score=sc,
                                  min_denominator=Int32(mn)))
                x = y + 1
            end
            kept = kept_ok[c] + kept_wrong[c]
            pushrow!(summary, (locus_id=locus.id, chrom=locus.chrom, L=Int32(L), F=Int32(F),
                               mapq_cutoff=Int32(cutoffs[c]), n_fragments=nfrag, n_skipped_N=nskip,
                               n_missing=nmiss, n_not_proper=nnp, n_kept=kept, n_correct_kept=kept_ok[c],
                               n_wrong_kept=kept_wrong[c], frac_wrong_kept=kept == 0 ? NaN : kept_wrong[c] / kept,
                               frac_correct=nfrag == 0 ? NaN : kept_ok[c] / nfrag, min_denominator=Int32(mind)))
            totals[(F, cutoffs[c])] = (nfrag, nskip, nmiss, nnp, kept_ok[c], kept_wrong[c], mind)
        end
    end
    return totals
end

function genome_wide_rows!(summary, L, all_totals)
    for key in sort!(collect(keys(first(all_totals))))
        t = [tot[key] for tot in all_totals]
        nfrag, nskip, nmiss, nnp, ok, wrong = (sum(x[i] for x in t) for i in 1:6)
        kept = ok + wrong
        pushrow!(summary, (locus_id="ALL", chrom="", L=Int32(L), F=Int32(key[1]), mapq_cutoff=Int32(key[2]),
                           n_fragments=nfrag, n_skipped_N=nskip, n_missing=nmiss, n_not_proper=nnp,
                           n_kept=kept, n_correct_kept=ok, n_wrong_kept=wrong,
                           frac_wrong_kept=kept == 0 ? NaN : wrong / kept,
                           frac_correct=nfrag == 0 ? NaN : ok / nfrag,
                           min_denominator=Int32(minimum(x[7] for x in t))))
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
    write_parquet(path, tbl)
end

"""
    score_simulation(simdir; bam_dir, out, cutoffs, tolerance=0, frag_mapq="min", fasta=nothing, read_lengths=nothing)

Phase 2 scoring. Reads `bam_dir/sim_L{L}.bam` for each simulated read length and writes
`mappability_tracks.parquet`, `misplacement_summary.parquet`, `bleed_matrix.parquet`,
`misplaced_sites.parquet` and `score.toml` to `out`.
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
    policy = parse_step_policy(grid["step_policy"])
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
    tiling = read_tiling(joinpath(simdir, "tiling.tsv"))
    checksum = meta["reference"]["sha256"]
    bychrom = Dict{String,Vector{Tuple{Int,Int,Int}}}()
    for (k, l) in enumerate(loci)
        push!(get!(Vector{Tuple{Int,Int,Int}}, bychrom, l.chrom), (l.start0, l.end0, k))
    end
    ctx = ScoreContext(loci, Dict(l.id => k for (k, l) in enumerate(loci)),
                       Dict(F => j for (j, F) in enumerate(Fs)), bychrom, tolerance, rule)

    fmax = grid["fmax"]
    regions = [sim_region(l, fmax, chrom_length(fa, l.chrom)) for l in loci]
    seqs = [fetch_seq(fa, l.chrom, r...) for (l, r) in zip(loci, regions)]

    results = Vector{Any}(undef, length(Ls))
    @info "Scoring $(length(loci)) loci" read_lengths = join(Ls, ",") cutoffs = join(cutoffs, ",") threads = Threads.nthreads()
    Threads.@threads :dynamic for t in eachindex(Ls)
        L = Ls[t]
        frags = [LocusFrags(l, seqs[k], regions[k][1], chrom_length(fa, l.chrom), L, Fs, policy)
                 for (k, l) in enumerate(loci)]
        for (k, l) in enumerate(loci), (j, F) in enumerate(Fs)
            exp = get(tiling, (l.id, L, F), nothing)
            exp === nothing && error("tiling.tsv has no row for ($(l.id), L=$L, F=$F)")
            lf = frags[k]
            got = (length(lf.starts[j]), count(==(ST_SKIPPED), view(lf.status, lf.offsets[j]+1:lf.offsets[j+1])) ÷ 2)
            exp == got || error("tiling mismatch for ($(l.id), L=$L, F=$F): simulate wrote (starts, N-skipped) = " *
                                "$exp, rescoring expects $got; was the reference or step policy changed?")
        end
        stats, bleed, sites, refnames = scan_bam!(frags, ctx, bams[L])
        tracks, summary, bleed_tbl, sites_tbl = new_tracks(), new_summary(), new_bleed(), new_sites()
        totals = [locus_scores!(tracks, summary, l, frags[k], L, Fs, cutoffs, checksum) for (k, l) in enumerate(loci)]
        genome_wide_rows!(summary, L, totals)
        bleed_rows!(bleed_tbl, sites_tbl, loci, L, cutoffs, bleed, sites, refnames)
        results[t] = (; stats, tracks, summary, bleed_tbl, sites_tbl)
        @info "Scored L=$L" bam = bams[L] pairs = stats.pairs incomplete = stats.incomplete
    end

    mkpath(out)
    all_tracks, all_summary, all_bleed, all_sites = new_tracks(), new_summary(), new_bleed(), new_sites()
    for r in results
        appendcols!(all_tracks, r.tracks)
        appendcols!(all_summary, r.summary)
        appendcols!(all_bleed, r.bleed_tbl)
        appendcols!(all_sites, r.sites_tbl)
    end
    write_table(joinpath(out, "mappability_tracks.parquet"), all_tracks)
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
    @info "Wrote scores" out track_rows = length(all_tracks.score)
    return (; tracks=all_tracks, summary=all_summary, bleed=all_bleed, sites=all_sites)
end
