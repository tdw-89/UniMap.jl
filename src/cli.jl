function cli_settings()
    s = ArgParseSettings(prog="unimap",
        description="Empirical paired-end mappability by simulated read tiling.",
        version=string(pkgversion(@__MODULE__)), add_version=true, commands_are_required=true)

    @add_arg_table! s begin
        "loci"
            help = "build target loci (gene span + strand-aware flanks) from an annotation"
            action = :command
        "simulate"
            help = "tile loci with error-free paired-end fragments and write gzipped FASTQ"
            action = :command
        "align-script"
            help = "write a bash script that aligns the simulated reads with bowtie2"
            action = :command
        "score"
            help = "classify aligned simulated pairs; write per-base tracks (Zarr) and summaries (Parquet)"
            action = :command
        "view"
            help = "interactive 3-D viewer of per-base mappability (GLMakie window)"
            action = :command
    end

    @add_arg_table! s["view"] begin
        "zarr"
            help = "path to scores/mappability.zarr, or a `simulate` directory containing scores/"
            required = true
        "--locus"
            help = "locus to show first"
        "--read-lengths", "-L"
            help = "read lengths to show (default: all in the store)"
        "--plot", "-p"
            help = "base (metric over fragment length x MAPQ cutoff at one base), track (over position x fragment length at one cutoff) or boxplot (distribution across bases per fragment length at one cutoff)"
            default = "base"
        "--cutoff", "-c"
            help = "MAPQ cutoff shown first by the track and boxplot views"
            arg_type = Int
            default = 2
    end

    @add_arg_table! s["align-script"] begin
        "--sim", "-s"
            help = "`unimap simulate` output directory; the script is written to SIM/align.sh"
            required = true
        "--index", "-x"
            help = "bowtie2 index prefix (built from the same FASTA as the simulation)"
            required = true
        "--threads", "-t"
            help = "bowtie2 threads"
            arg_type = Int
            default = 4
        "--aligner-args"
            help = "production bowtie2 arguments as one string; use the = form, e.g. --aligner-args='--very-sensitive -X 500'. -X defaults to the largest simulated fragment length"
            default = ""
        "--bowtie2"
            help = "bowtie2 executable"
            default = "bowtie2"
        "--samtools"
            help = "samtools executable"
            default = "samtools"
    end

    @add_arg_table! s["score"] begin
        "--sim", "-s"
            help = "`unimap simulate` output directory"
            required = true
        "--bam-dir"
            help = "directory containing sim_L{L}.bam, or shards sim_L{L}_*.bam scanned in parallel (default: SIM/aln)"
        "--out", "-o"
            help = "output directory (default: SIM/scores)"
        "--cutoffs"
            help = "comma-separated fragment MAPQ cutoffs"
            default = DEFAULT_CUTOFFS
        "--tolerance"
            help = "allowed bp difference from the origin position and fragment length"
            arg_type = Int
            default = 0
        "--frag-mapq"
            help = "fragment MAPQ from the two mates: min (conservative) or max; real-data counting must use the same rule"
            default = "min"
        "--read-lengths", "-L"
            help = "score only these read lengths (default: all simulated)"
        "--fasta", "-f"
            help = "reference FASTA, if it has moved since `simulate`"
    end

    @add_arg_table! s["loci"] begin
        "--fasta", "-f"
            help = "reference FASTA (with samtools .fai index)"
            required = true
        "--annotation", "-a"
            help = "GFF3 or GTF (optionally .gz)"
            required = true
        "--genes", "-g"
            help = "file with one ID per line, or comma-separated IDs (gene ID/name, transcript ID or protein ID)"
            required = true
        "--upstream"
            help = "flank upstream of the TSS (bp)"
            arg_type = Int
            default = 1500
        "--downstream"
            help = "flank downstream of the TES (bp)"
            arg_type = Int
            default = 1500
        "--out", "-o"
            help = "output loci TSV"
            default = "loci.tsv"
    end

    @add_arg_table! s["simulate"] begin
        "--fasta", "-f"
            help = "reference FASTA (with samtools .fai index); must be the one the aligner index was built from"
            required = true
        "--loci", "-l"
            help = "loci TSV from `unimap loci`, or a BED of loci (used as-is)"
        "--annotation", "-a"
            help = "GFF3 or GTF; alternative to --loci (requires --genes)"
        "--genes", "-g"
            help = "IDs to resolve against --annotation (file or comma-separated)"
        "--upstream"
            help = "flank upstream of the TSS (bp), with --annotation"
            arg_type = Int
            default = 1500
        "--downstream"
            help = "flank downstream of the TES (bp), with --annotation"
            arg_type = Int
            default = 1500
        "--out", "-o"
            help = "output directory"
            required = true
        "--read-lengths", "-L"
            help = "comma-separated read lengths, e.g. 50,100,150"
            required = true
        "--fragment-lengths", "-F"
            help = "comma-separated fragment lengths, e.g. 100,150,200,300,500; the largest should equal the aligner's max fragment length (bowtie2 -X)"
            required = true
        "--step"
            help = "bp between consecutive fragment starts (larger = fewer reads; keep well below the shortest read length)"
            arg_type = Int
            default = 1
        "--orientation"
            help = "which read orientations to simulate per fragment: both, plus (R1 forward at fragment start) or minus (R1 = revcomp of fragment end)"
            default = "both"
        "--qual"
            help = "uniform base quality character"
            default = "I"
        "--no-checksum"
            help = "skip SHA-256 of the reference FASTA"
            action = :store_true
    end
    return s
end

function load_loci(args, fa::IndexedFasta)
    if args["loci"] !== nothing
        args["annotation"] === nothing || error("give either --loci or --annotation/--genes, not both")
        return read_loci(args["loci"])
    end
    (args["annotation"] === nothing || args["genes"] === nothing) &&
        error("give --loci, or both --annotation and --genes")
    @info "Reading annotation $(args["annotation"])"
    ann = read_annotation(args["annotation"])
    return build_loci(ann, read_id_list(args["genes"]), fa; upstream=args["upstream"], downstream=args["downstream"])
end

function cmd_loci(args)
    fa = IndexedFasta(args["fasta"])
    loci = load_loci(merge(args, Dict("loci" => nothing)), fa)
    write_loci(args["out"], loci)
    @info "Wrote $(length(loci)) loci to $(args["out"])"
end

function cmd_simulate(args)
    fa = IndexedFasta(args["fasta"])
    loci = load_loci(args, fa)
    length(args["qual"]) == 1 || error("--qual must be a single character")
    provenance = Dict{String,Any}("command" => "simulate",
        "loci_source" => something(args["loci"], args["annotation"]),
        "genes" => something(args["genes"], ""),
        "upstream" => args["upstream"], "downstream" => args["downstream"])
    simulate(fa, loci, args["out"];
             read_lengths=parse_int_list(args["read-lengths"]),
             fragment_lengths=parse_int_list(args["fragment-lengths"]),
             step=args["step"], orientation=args["orientation"],
             qual=only(args["qual"]), checksum=!args["no-checksum"], provenance)
end

function cmd_align_script(args)
    path = write_align_script(args["sim"]; index=args["index"], threads=args["threads"],
                              aligner_args=args["aligner-args"], bowtie2=args["bowtie2"],
                              samtools=args["samtools"])
    @info "Wrote $path; run it with `bash $path`"
end

function cmd_score(args)
    sim = args["sim"]
    score_simulation(sim; bam_dir=something(args["bam-dir"], joinpath(sim, "aln")),
                     out=something(args["out"], joinpath(sim, "scores")),
                     cutoffs=parse_int_list(args["cutoffs"]), tolerance=args["tolerance"],
                     frag_mapq=args["frag-mapq"], fasta=args["fasta"],
                     read_lengths=args["read-lengths"] === nothing ? nothing : parse_int_list(args["read-lengths"]))
end

function (@main)(argv)
    args = parse_args(argv, cli_settings())
    cmd = args["%COMMAND%"]
    sub = args[cmd]
    if cmd == "loci"
        cmd_loci(sub)
    elseif cmd == "simulate"
        cmd_simulate(sub)
    elseif cmd == "align-script"
        cmd_align_script(sub)
    elseif cmd == "score"
        cmd_score(sub)
    elseif cmd == "view"
        path = sub["zarr"]
        endswith(rstrip(path, '/'), ".zarr") || (path = joinpath(path, "scores", "mappability.zarr"))
        Ls = sub["read-lengths"] === nothing ? nothing : parse_int_list(sub["read-lengths"])
        view_mappability(path; locus=sub["locus"], read_lengths=Ls, plot=Symbol(sub["plot"]),
                         mapq_cutoff=sub["cutoff"])
    end
    return 0
end
