using Test
using UniMap
using Random

const GzipDecompressorStream = UniMap.CodecZlib.GzipDecompressorStream

revcomp(s::AbstractString) = String(reverse(map(c -> Dict('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A', 'N' => 'N')[c], s)))

function write_fasta_with_fai(path, seqs::Vector{Pair{String,String}}; width=60)
    open(path, "w") do io
        for (name, s) in seqs
            println(io, ">", name, " synthetic")
            for i in 1:width:length(s)
                println(io, s[i:min(i + width - 1, end)])
            end
        end
    end
    open(path * ".fai", "w") do fai
        offset = 0
        for (name, s) in seqs
            offset += length(">$name synthetic\n")
            println(fai, join((name, length(s), offset, width, width + 1), '\t'))
            offset += length(s) + cld(length(s), width)
        end
    end
end

function read_fastq_gz(path)
    lines = readlines(GzipDecompressorStream(open(path)))
    return [(name=lines[i][2:end], seq=lines[i+1], qual=lines[i+3]) for i in 1:4:length(lines)]
end

const DIR = mktempdir()
const RNG = MersenneTwister(1)
const CHR1 = String(rand(RNG, collect("ACGT"), 200_000))
const CHR2 = let s = collect(String(rand(RNG, collect("ACGT"), 20_000)))
    s[10_001:10_010] .= 'N'
    String(s)
end
const FASTA = joinpath(DIR, "genome.fa")
# soft-masked stretch to check uppercasing
write_fasta_with_fai(FASTA, ["chr1" => CHR1[1:999] * lowercase(CHR1[1000:1100]) * CHR1[1101:end], "chr2" => CHR2])

@testset "UniMap" begin

@testset "length lists and orientation" begin
    @test parse_int_list("150, 50,100,50") == [50, 100, 150]
    @test UniMap.check_lengths("x", [3, 1, 2]) == [1, 2, 3]
    @test_throws ArgumentError UniMap.check_lengths("x", Int[])
    @test_throws ArgumentError UniMap.check_lengths("x", [0, 5])
    @test parse_orientation("both") == ['+', '-']
    @test parse_orientation("plus") == ['+'] && parse_orientation("minus") == ['-']
    @test_throws ArgumentError parse_orientation("+")
end

@testset "indexed FASTA" begin
    fa = IndexedFasta(FASTA)
    @test chrom_length(fa, "chr1") == 200_000
    @test String(fetch_seq(fa, "chr1", 0, 200_000)) == CHR1
    @test String(fetch_seq(fa, "chr1", 59, 121)) == CHR1[60:121]
    @test String(fetch_seq(fa, "chr2", 10_000, 10_010)) == "N"^10
    @test isempty(fetch_seq(fa, "chr1", 5, 5))
    @test_throws ArgumentError fetch_seq(fa, "chr1", 199_990, 200_001)
end

@testset "annotation → loci" begin
    gff = joinpath(DIR, "ann.gff3")
    write(gff, """
    ##gff-version 3
    chr1\tsrc\tregion\t1\t200000\t.\t+\t.\tID=chr1:1..200000;Name=chr1
    chr1\tsrc\tgene\t50001\t56000\t.\t+\t.\tID=gene-geneA;Name=geneA;gene=geneA
    chr1\tsrc\tmRNA\t50001\t55000\t.\t+\t.\tID=rna-TX1;Parent=gene-geneA;transcript_id=TX1
    chr1\tsrc\texon\t50001\t51000\t.\t+\t.\tID=exon-TX1-1;Parent=rna-TX1
    chr1\tsrc\texon\t53001\t55000\t.\t+\t.\tID=exon-TX1-2;Parent=rna-TX1
    chr1\tsrc\tCDS\t50101\t51000\t.\t+\t0\tID=cds-P1;Parent=rna-TX1;protein_id=P1
    chr1\tsrc\tmRNA\t50501\t56000\t.\t+\t.\tID=rna-TX2;Parent=gene-geneA;transcript_id=TX2
    chr1\tsrc\texon\t50501\t56000\t.\t+\t.\tID=exon-TX2-1;Parent=rna-TX2
    chr1\tsrc\tgene\t56501\t60000\t.\t-\t.\tID=gene-geneB;Name=geneB
    chr1\tsrc\tmRNA\t56501\t60000\t.\t-\t.\tID=rna-TX3;Parent=gene-geneB
    chr1\tsrc\texon\t56501\t60000\t.\t-\t.\tID=exon-TX3-1;Parent=rna-TX3
    chr1\tsrc\tgene\t150001\t151000\t.\t+\t.\tID=gene-dup;Name=dupname
    chr1\tsrc\tgene\t160001\t161000\t.\t+\t.\tID=gene-dup2;Name=dupname
    """)
    ann = read_annotation(gff)
    fa = IndexedFasta(FASTA)
    @test ann.genes["gene-geneA"].start0 == 50_000
    @test ann.genes["gene-geneA"].end0 == 56_000
    @test !haskey(ann.genes, "chr1:1..200000")
    loci = @test_logs (:warn,) (:warn,) build_loci(ann, ["P1", "TX2", "geneB", "nope", "dupname"], fa;
                                                    upstream=1000, downstream=2000)
    @test length(loci) == 2
    a, b = loci
    @test (a.id, a.query, a.strand) == ("geneA", "P1,TX2", '+')
    @test (a.start0, a.end0, a.gene_start0, a.gene_end0) == (49_000, 58_000, 50_000, 56_000)
    @test a.upstream_overlaps == "" && a.downstream_overlaps == "geneB"
    @test (b.id, b.strand) == ("geneB", '-')
    @test (b.start0, b.end0) == (56_500 - 2000, 60_000 + 1000)
    @test b.upstream_overlaps == "" && b.downstream_overlaps == "geneA"

    tsv = joinpath(DIR, "loci.tsv")
    write_loci(tsv, loci)
    @test read_loci(tsv) == loci

    gtf = joinpath(DIR, "ann.gtf")
    write(gtf, """
    chr1\tsrc\texon\t1001\t2000\t.\t-\t.\tgene_id "G1"; transcript_id "T1"; gene_name "g1";
    chr1\tsrc\texon\t3001\t4000\t.\t-\t.\tgene_id "G1"; transcript_id "T2"; gene_name "g1";
    chr1\tsrc\tCDS\t3001\t3500\t.\t-\t0\tgene_id "G1"; transcript_id "T2"; protein_id "PR2";
    """)
    l = only(build_loci(read_annotation(gtf), ["PR2"], fa))
    @test (l.id, l.start0, l.end0, l.strand) == ("G1", 0, 5500, '-')
end

@testset "BED loci" begin
    bed = joinpath(DIR, "loci.bed")
    write(bed, "track name=x\nchr1\t100\t200\tmy|locus\t0\t-\nchr2\t5\t50\n")
    l = read_loci(bed)
    @test l[1].id == "my_locus" && l[1].strand == '-' && (l[1].start0, l[1].end0) == (100, 200)
    @test l[2].id == "chr2:5-50"
end

@testset "tiling" begin
    fa = IndexedFasta(FASTA)
    locus = Locus("loc1", "chr1", 100_000, 100_400, '+', 100_000, 100_400, "loc1", "", "")
    out = joinpath(DIR, "sim1")
    Ls = [40, 100]
    Fs = [30, 40, 60, 100, 150]
    rows = simulate(fa, [locus], out; read_lengths=Ls, fragment_lengths=Fs, checksum=false)
    @test length(rows) == length(Ls) * length(Fs)
    @test all(r -> r.n_skipped_N == 0, rows)
    @test isfile(joinpath(out, "simulate.toml")) && isfile(joinpath(out, "tiling.tsv"))

    for L in Ls
        r1 = read_fastq_gz(joinpath(out, "fastq", "L$L", "loc1_R1.fq.gz"))
        r2 = read_fastq_gz(joinpath(out, "fastq", "L$L", "loc1_R2.fq.gz"))
        @test length(r1) == length(r2) == sum(r.n_pairs for r in rows if r.L == L)
        @test all(i -> r1[i].name == r2[i].name, eachindex(r1))
        @test all(r -> !occursin(r"\s", r.name) && length(r.name) <= 254, r1)
        @test all(r -> r.qual == "I"^length(r.seq), r1)

        cov = Dict(F => zeros(Int, length(locus)) for F in Fs)
        ok = true
        for (a, b) in zip(r1, r2)
            n = parse_read_name(a.name)
            m = min(L, n.F)
            frag = CHR1[n.frag_start0+1:n.frag_start0+n.F]
            left, right = frag[1:m], revcomp(frag[end-m+1:end])
            # F < L trimming: mates are the whole fragment
            ok &= length(a.seq) == length(b.seq) == m
            ok &= n.orient == '+' ? (a.seq, b.seq) == (left, right) : (a.seq, b.seq) == (right, left)
            for x in max(n.frag_start0, locus.start0):min(n.frag_start0 + n.F, locus.end0)-1
                cov[n.F][x-locus.start0+1] += 1
            end
        end
        @test ok
        # every scored position is covered by F starts × 2 orientations
        @test all(F -> all(==(2F), cov[F]), Fs)
    end
    @test_throws ArgumentError simulate(fa, [locus], out; read_lengths=Int[], fragment_lengths=Fs, checksum=false)
    @test_throws ArgumentError simulate(fa, [locus], out; read_lengths=Ls, fragment_lengths=Fs, step=0, checksum=false)
end

@testset "single orientation and thinned tiling" begin
    fa = IndexedFasta(FASTA)
    locus = Locus("loc2", "chr1", 120_000, 120_300, '+', 120_000, 120_300, "loc2", "", "")
    both = simulate(fa, [locus], joinpath(DIR, "o_both"); read_lengths=[50], fragment_lengths=[100], checksum=false)
    for (o, c) in (("plus", '+'), ("minus", '-'))
        out = joinpath(DIR, "o_$o")
        r = only(simulate(fa, [locus], out; read_lengths=[50], fragment_lengths=[100], orientation=o, checksum=false))
        @test r.n_pairs == only(both).n_pairs ÷ 2
        reads = read_fastq_gz(joinpath(out, "fastq", "L50", "loc2_R1.fq.gz"))
        @test length(reads) == r.n_pairs && all(x -> parse_read_name(x.name).orient == c, reads)
        @test occursin("orientation = \"$o\"", read(joinpath(out, "simulate.toml"), String))
    end
    @test_throws ArgumentError simulate(fa, [locus], joinpath(DIR, "o_bad"); read_lengths=[50],
                                       fragment_lengths=[100], orientation="sideways", checksum=false)

    r = only(simulate(fa, [locus], joinpath(DIR, "thin"); read_lengths=[50], fragment_lengths=[100], step=10, checksum=false))
    @test r.step == 10 && r.n_starts == length(fragment_starts(locus, 100, 10, 200_000)) == cld(300 + 99, 10)
    @test r.n_pairs == 2r.n_starts
end

@testset "N fragments are skipped" begin
    fa = IndexedFasta(FASTA)
    locus = Locus("nloc", "chr2", 9_900, 10_100, '+', 9_900, 10_100, "nloc", "", "")
    rows = simulate(fa, [locus], joinpath(DIR, "sim2"); read_lengths=[50], fragment_lengths=[60], checksum=false)
    r = only(rows)
    # starts whose first or last 50 bp touch [10000, 10010)
    expected = count(f0 -> (f0 < 10_010 && f0 + 50 > 10_000) || (f0 + 10 < 10_010 && f0 + 60 > 10_000),
                     fragment_starts(locus, 60, 1, 20_000))
    @test r.n_skipped_N == expected > 0
    reads = read_fastq_gz(joinpath(DIR, "sim2", "fastq", "L50", "nloc_R1.fq.gz"))
    @test length(reads) == r.n_pairs && all(x -> !occursin('N', x.seq), reads)
end

@testset "chromosome edge clipping" begin
    locus = Locus("edge", "chr2", 0, 100, '+', 0, 100, "edge", "", "")
    starts = fragment_starts(locus, 80, 1, 20_000)
    @test first(starts) == 0 && last(starts) == 99
end

@testset "threaded tiling summaries" begin
    fa = IndexedFasta(FASTA)
    loci = [Locus("t$k", "chr1", 10_000k, 10_000k + 300 + 100k, '+', 0, 0, "", "", "") for k in 1:6]
    Fs = [50, 75, 150, 200]
    rows = simulate(fa, loci, joinpath(DIR, "sim3"); read_lengths=[40, 100], fragment_lengths=Fs, step=3, checksum=false)
    for l in loci, F in Fs
        n = [r.n_starts for r in rows if r.locus_id == l.id && r.F == F]
        @test length(n) == 2 && n[1] == n[2] == length(fragment_starts(l, F, 3, 200_000))
    end
end

@testset "CLI" begin
    gff = joinpath(DIR, "ann.gff3")
    out = joinpath(DIR, "cli")
    @test UniMap.main(["simulate", "-f", FASTA, "-a", gff, "-g", "geneA,geneB", "-o", out,
                       "-L", "50", "-F", "100,120"]) == 0
    @test isfile(joinpath(out, "fastq", "L50", "geneA_R1.fq.gz"))
    @test isfile(joinpath(out, "fastq", "L50", "geneB_R2.fq.gz"))
    @test length(read_loci(joinpath(out, "loci.tsv"))) == 2
    @test occursin("sha256", read(joinpath(out, "simulate.toml"), String))

    @test UniMap.main(["align-script", "-s", out, "-x", joinpath(DIR, "noindex"), "-t", "2",
                       "--aligner-args=--very-sensitive --no-mixed"]) == 0
    sh = read(joinpath(out, "align.sh"), String)
    @test occursin("ALIGNER_ARGS=(--very-sensitive --no-mixed -X 120)", sh)
    @test occursin("READ_LENGTHS=(50)", sh)
    @test_throws ErrorException write_align_script(out; index="x", aligner_args="--no-unal")
end

@testset "pair classification" begin
    M = UniMap.Mate
    refs = ["chrA", "chrB"]
    plus = (chrom="chrA", frag_start0=100, F=200, orient='+')
    minus = (; plus..., orient='-')
    C, W, N = UniMap.ST_CORRECT, UniMap.ST_WRONG, UniMap.ST_NOT_PROPER
    @test UniMap.classify_pair(M(99, 1, 100, 42, 200), M(147, 1, 250, 42, -200), plus, refs) == C
    @test UniMap.classify_pair(M(83, 1, 250, 42, -200), M(163, 1, 100, 42, 200), minus, refs) == C
    @test UniMap.classify_pair(M(83, 1, 250, 42, -200), M(163, 1, 100, 42, 200), plus, refs) == W
    @test UniMap.classify_pair(M(99, 1, 101, 42, 200), M(147, 1, 251, 42, -200), plus, refs) == W
    @test UniMap.classify_pair(M(99, 1, 101, 42, 200), M(147, 1, 251, 42, -200), plus, refs; tolerance=1) == C
    @test UniMap.classify_pair(M(99, 2, 100, 42, 200), M(147, 2, 250, 42, -200), plus, refs) == W
    @test UniMap.classify_pair(M(97, 1, 100, 42, 200), M(145, 1, 250, 42, -200), plus, refs) == N
    @test UniMap.classify_pair(M(99, 1, 100, 42, 200), nothing, plus, refs) == N
end

have_aligner = all(x -> Sys.which(x) !== nothing, ("bowtie2", "bowtie2-build", "samtools"))
have_aligner || @warn "bowtie2/bowtie2-build/samtools not on PATH; skipping alignment tests"

if have_aligner
@testset "align + score (synthetic genome)" begin
    rng = MersenneTwister(7)
    dir = mktempdir()
    chrA = collect(String(rand(rng, collect("ACGT"), 300_000)))
    chrB = collect(String(rand(rng, collect("ACGT"), 60_000)))
    chrA[150_001:151_000] = chrA[50_001:51_000]             # exact 1 kb duplicate
    near = chrA[80_001:82_000]
    for i in 100:200:2000                                    # SNP every 200 bp
        near[i] = near[i] == 'A' ? 'C' : 'A'
    end
    chrB[30_001:32_000] = near
    fasta = joinpath(dir, "g.fa")
    write_fasta_with_fai(fasta, ["chrA" => String(chrA), "chrB" => String(chrB)])
    run(pipeline(`bowtie2-build -q $fasta $(joinpath(dir, "g"))`; stdout=devnull))

    mk(id, chrom, s, e) = Locus(id, chrom, s, e, '+', s, e, id, "", "")
    loci = [mk("uniq", "chrA", 20_000, 22_000), mk("dupA", "chrA", 49_500, 51_500),
            mk("dupB", "chrA", 149_500, 151_500), mk("nearA", "chrA", 79_500, 82_500),
            mk("nearB", "chrB", 29_500, 32_500)]
    sim = joinpath(dir, "sim")
    Fs = [100, 150, 200, 300, 400]
    simulate(IndexedFasta(fasta), loci, sim; read_lengths=[50], fragment_lengths=Fs, checksum=false)
    script = write_align_script(sim; index=joinpath(dir, "g"), threads=4)
    run(pipeline(`bash $script`; stdout=devnull))
    @test isfile(joinpath(sim, "aln", "sim_L50.bam"))

    res = score_simulation(sim)
    for f in ("misplacement_summary", "bleed_matrix", "misplaced_sites")
        @test isfile(joinpath(sim, "scores", "$f.parquet"))
    end
    zpath = joinpath(sim, "scores", "mappability.zarr")
    @test res.tracks == zpath && isdir(zpath)
    s, b = res.summary, res.bleed
    rows(tbl, pred) = [i for i in eachindex(first(tbl)) if pred(i)]
    cuts = [0, 2, 5, 10, 20, 30]
    function at(id, F, c, x)
        tr = locus_track(zpath, id; L=50, F, mapq_cutoff=c)
        return tr.score[x-tr.start0+1]
    end

    z = UniMap.Zarr.zopen(zpath; consolidated=true)
    @test size(z["score"]) == (sum(length, loci), length(Fs), length(cuts), 1)
    @test z["L"][:] == [50] && z["F"][:] == Fs && z["mapq_cutoff"][:] == cuts
    @test [l["offset"] for l in z.attrs["loci"]] == [0; cumsum(length.(loci))[1:end-1]]

    # unique sequence: every position maps back at every cutoff; denominator constant (= 2F)
    for F in Fs, c in cuts
        tr = locus_track(z, "uniq"; L=50, F, mapq_cutoff=c)
        @test length(tr.score) == 2000 && all(==(1.0f0), tr.score) && all(==(2F), tr.denominator)
    end
    @test_throws ErrorException locus_track(z, "uniq"; L=50, F=99, mapq_cutoff=0)
    @test all(==(0), s.n_missing) && all(==(0), s.n_not_proper)

    # exact duplicate: ~0 inside at cutoff 2, recovered near the edge by long fragments
    @test at("dupA", 100, 2, 50_500) < 0.05
    @test at("dupA", 400, 2, 50_100) > 0.5
    @test at("dupA", 400, 2, 50_100) > at("dupA", 100, 2, 50_100)
    @test at("dupA", 100, 0, 49_600) == 1.0

    # near duplicate: intermediate, with bleed into the copy at low cutoffs
    sA = rows(s, i -> s.locus_id[i] == "nearA" && s.F[i] == 100 && s.mapq_cutoff[i] == 0)
    @test 0.2 < s.frac_correct[only(sA)] < 0.99
    @test s.n_wrong_kept[only(sA)] > 0
    bleed_to(o, d, c) = sum((b.n_fragments[i] for i in rows(b, i -> b.origin_locus[i] == o &&
                             b.dest_locus[i] == d && b.mapq_cutoff[i] == c)); init=0)
    @test bleed_to("nearA", "nearB", 0) > 0
    @test bleed_to("dupA", "dupB", 0) > 0
    @test bleed_to("uniq", "elsewhere", 0) == 0

    g = rows(s, i -> s.locus_id[i] == "ALL")
    @test length(g) == length(Fs) * 6

    # single orientation + thinned tiling: same scores in unique sequence, denominator F/step
    sim1 = joinpath(dir, "sim_minus")
    simulate(IndexedFasta(fasta), loci[1:2], sim1; read_lengths=[50], fragment_lengths=Fs,
             step=2, orientation="minus", checksum=false)
    run(pipeline(`bash $(write_align_script(sim1; index=joinpath(dir, "g"), threads=4))`; stdout=devnull))
    z1 = score_simulation(sim1; cutoffs=[2]).tracks
    for F in Fs
        tr = locus_track(z1, "uniq"; L=50, F, mapq_cutoff=2)
        @test all(==(1.0f0), tr.score) && all(==(F ÷ 2), tr.denominator)
    end
    @test locus_track(z1, "dupA"; L=50, F=100, mapq_cutoff=2).score[50_500-49_500+1] < 0.05
end
end

end
