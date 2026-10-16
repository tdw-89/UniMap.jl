"""
    Locus

A target region. `start0`/`end0` (0-based half-open) is the scored region: gene span
plus strand-aware flanks. `gene_start0`/`gene_end0` is the gene span alone.
`upstream_overlaps`/`downstream_overlaps` list other genes overlapping each flank.
"""
struct Locus
    id::String
    chrom::String
    start0::Int
    end0::Int
    strand::Char
    gene_start0::Int
    gene_end0::Int
    query::String
    upstream_overlaps::String
    downstream_overlaps::String
end

Base.length(l::Locus) = l.end0 - l.start0

struct GeneRecord
    id::String
    chrom::String
    strand::Char
    start0::Int
    end0::Int
end

"""
    Annotation

Gene spans (union of each gene's exons, i.e. outermost TSS → outermost TES) and an
alias table mapping gene IDs/names, transcript IDs and protein IDs to gene IDs.
"""
struct Annotation
    genes::Dict{String,GeneRecord}
    aliases::Dict{String,Set{String}}
end

function open_text(path::AbstractString)
    io = open(path)
    return endswith(path, ".gz") ? GzipDecompressorStream(io) : io
end

"""
    read_annotation(path) -> Annotation

Reads GTF (`.gtf[.gz]`) or GFF3 (anything else, e.g. `.gff`, `.gff3[.gz]`).
"""
function read_annotation(path::AbstractString)
    isfile(path) || error("annotation not found: $path")
    base = endswith(path, ".gz") ? path[1:end-3] : path
    io = open_text(path)
    try
        return endswith(lowercase(base), ".gtf") ? parse_gtf(io) : parse_gff3(io)
    finally
        close(io)
    end
end

function gff3_attr(attrs::AbstractString, key::AbstractString)
    for field in eachsplit(attrs, ';')
        field = lstrip(field)
        if startswith(field, key) && ncodeunits(field) > ncodeunits(key) && codeunit(field, ncodeunits(key) + 1) == UInt8('=')
            return SubString(field, ncodeunits(key) + 2)
        end
    end
    return nothing
end

function gtf_attr(attrs::AbstractString, key::AbstractString)
    for field in eachsplit(attrs, ';')
        field = strip(field)
        if startswith(field, key) && ncodeunits(field) > ncodeunits(key) && codeunit(field, ncodeunits(key) + 1) == UInt8(' ')
            return strip(SubString(field, ncodeunits(key) + 2), '"')
        end
    end
    return nothing
end

const Span = Tuple{String,Char,Int,Int}  # chrom, strand, start0, end0

function merge_span!(d::Dict{String,Span}, key::AbstractString, chrom, strand::Char, s0::Int, e0::Int)
    old = get(d, key, nothing)
    if old === nothing
        d[String(key)] = (String(chrom), strand, s0, e0)
    elseif old[1] == chrom
        d[key] = (old[1], old[2], min(old[3], s0), max(old[4], e0))
    end
    return d
end

add_alias!(aliases, alias, target) = push!(get!(Set{String}, aliases, String(alias)), String(target))

const NCBI_ID_PREFIXES = ("gene-", "rna-", "cds-")

function parse_gff3(io::IO)
    parent = Dict{String,String}()
    toplevel = Dict{String,Tuple{Span,String}}()  # id => (span, feature type)
    exon_spans = Dict{String,Span}()               # exon Parent => union of its exons
    feature_aliases = Dict{String,Set{String}}()   # alias => feature IDs
    for line in eachline(io)
        (isempty(line) || line[1] == '#') && continue
        f = split(line, '\t')
        length(f) < 9 && continue
        typ, attrs = f[3], f[9]
        s0, e0 = parse(Int, f[4]) - 1, parse(Int, f[5])
        strand = f[7][1]
        par = gff3_attr(attrs, "Parent")
        par === nothing || (par = first(eachsplit(par, ',')))
        if typ == "exon"
            par === nothing || merge_span!(exon_spans, par, f[1], strand, s0, e0)
            continue
        elseif typ == "CDS"
            pid = gff3_attr(attrs, "protein_id")
            (pid === nothing || par === nothing) || add_alias!(feature_aliases, pid, par)
            continue
        end
        id = gff3_attr(attrs, "ID")
        id === nothing && continue
        if par === nothing
            toplevel[id] = ((String(f[1]), strand, s0, e0), String(typ))
        else
            parent[id] = par
        end
        add_alias!(feature_aliases, id, id)
        for p in NCBI_ID_PREFIXES
            startswith(id, p) && add_alias!(feature_aliases, id[length(p)+1:end], id)
        end
        for key in ("Name", "gene", "gene_id", "transcript_id", "locus_tag")
            v = gff3_attr(attrs, key)
            v === nothing || add_alias!(feature_aliases, v, id)
        end
    end

    function root(id)
        for _ in 1:100
            p = get(parent, id, nothing)
            p === nothing && return id
            id = p
        end
        error("cyclic Parent chain at '$id'")
    end

    gene_spans = Dict{String,Span}()
    for (feat, (chrom, strand, s0, e0)) in exon_spans
        merge_span!(gene_spans, root(feat), chrom, strand, s0, e0)
    end
    for (id, (span, typ)) in toplevel
        haskey(gene_spans, id) || !occursin("gene", lowercase(typ)) || (gene_spans[id] = span)
    end
    genes = Dict(id => GeneRecord(id, s...) for (id, s) in gene_spans)

    aliases = Dict{String,Set{String}}()
    for (alias, feats) in feature_aliases, feat in feats
        g = root(feat)
        haskey(genes, g) && add_alias!(aliases, alias, g)
    end
    return Annotation(genes, aliases)
end

function parse_gtf(io::IO)
    exon_spans = Dict{String,Span}()
    gene_feature_spans = Dict{String,Span}()
    aliases = Dict{String,Set{String}}()
    for line in eachline(io)
        (isempty(line) || line[1] == '#') && continue
        f = split(line, '\t')
        length(f) < 9 && continue
        typ, attrs = f[3], f[9]
        gid = gtf_attr(attrs, "gene_id")
        gid === nothing && continue
        s0, e0 = parse(Int, f[4]) - 1, parse(Int, f[5])
        strand = f[7][1]
        if typ == "exon"
            merge_span!(exon_spans, gid, f[1], strand, s0, e0)
        elseif typ == "gene"
            merge_span!(gene_feature_spans, gid, f[1], strand, s0, e0)
        end
        add_alias!(aliases, gid, gid)
        for key in ("gene_name", "transcript_id", "transcript_name", "protein_id")
            v = gtf_attr(attrs, key)
            v === nothing || add_alias!(aliases, v, gid)
        end
    end
    for (gid, span) in gene_feature_spans
        haskey(exon_spans, gid) || (exon_spans[gid] = span)
    end
    genes = Dict(id => GeneRecord(id, s...) for (id, s) in exon_spans)
    return Annotation(genes, aliases)
end

"""Turn an annotation gene ID into a locus ID safe for read names and file names."""
function locus_id_from_gene(id::AbstractString)
    for p in NCBI_ID_PREFIXES
        startswith(id, p) && (id = id[length(p)+1:end]; break)
    end
    return replace(id, r"[\s|/]" => "_")
end

"""
    read_id_list(s) -> Vector{String}

IDs from a file (one per line; `#` comments allowed) or a comma-separated string.
"""
function read_id_list(s::AbstractString)
    lines = isfile(s) ? readlines(s) : split(s, ',')
    ids = String[]
    for l in lines
        l = strip(first(split(l, '#')))
        isempty(l) || push!(ids, l)
    end
    return unique!(ids)
end

function overlapping_genes(bychrom, chrom, s0, e0, self_id)
    s0 >= e0 && return ""
    hits = [g.id for g in get(bychrom, chrom, GeneRecord[]) if g.id != self_id && g.start0 < e0 && s0 < g.end0]
    return join(sort!(locus_id_from_gene.(hits)), ",")
end

"""
    build_loci(ann, queries, fa; upstream=1500, downstream=1500) -> Vector{Locus}

Resolve `queries` (gene/transcript/protein IDs or gene names) to genes and add
strand-aware flanks (upstream of a minus-strand gene is at higher coordinates).
Flanks are clipped at chromosome ends. Unresolved or ambiguous queries are warned about and skipped.
"""
function build_loci(ann::Annotation, queries::AbstractVector{<:AbstractString}, fa::IndexedFasta;
                    upstream::Integer=1500, downstream::Integer=1500)
    order = String[]
    gene_queries = Dict{String,Vector{String}}()
    missing_q, ambiguous_q = String[], String[]
    for q in queries
        gs = get(ann.aliases, q, nothing)
        if gs === nothing
            push!(missing_q, q)
        elseif length(gs) > 1
            push!(ambiguous_q, "$q => " * join(sort!(collect(gs)), "/"))
        else
            g = only(gs)
            haskey(gene_queries, g) || push!(order, g)
            push!(get!(Vector{String}, gene_queries, g), q)
        end
    end
    isempty(missing_q) || @warn "$(length(missing_q)) ID(s) not found in annotation" missing_q
    isempty(ambiguous_q) || @warn "$(length(ambiguous_q)) ID(s) match multiple genes; skipped" ambiguous_q
    isempty(order) && error("none of the requested IDs were found in the annotation")

    bychrom = Dict{String,Vector{GeneRecord}}()
    for g in values(ann.genes)
        push!(get!(Vector{GeneRecord}, bychrom, g.chrom), g)
    end

    loci = Locus[]
    for gid in order
        g = ann.genes[gid]
        haskey(fa, g.chrom) || error("gene $gid is on '$(g.chrom)', which is not in the reference FASTA")
        len = chrom_length(fa, g.chrom)
        minus = g.strand == '-'
        left, right = minus ? (downstream, upstream) : (upstream, downstream)
        s0, e0 = max(0, g.start0 - left), min(len, g.end0 + right)
        up = minus ? (g.end0, e0) : (s0, g.start0)
        down = minus ? (s0, g.start0) : (g.end0, e0)
        push!(loci, Locus(locus_id_from_gene(gid), g.chrom, s0, e0, minus ? '-' : '+',
                          g.start0, g.end0, join(gene_queries[gid], ","),
                          overlapping_genes(bychrom, g.chrom, up..., gid),
                          overlapping_genes(bychrom, g.chrom, down..., gid)))
    end
    return loci
end

const LOCI_HEADER = ["locus_id", "chrom", "start0", "end0", "strand", "gene_start0", "gene_end0",
                     "query", "upstream_flank_overlaps", "downstream_flank_overlaps"]

function write_loci(path::AbstractString, loci::AbstractVector{Locus})
    open(path, "w") do io
        println(io, join(LOCI_HEADER, '\t'))
        for l in loci
            println(io, join((l.id, l.chrom, l.start0, l.end0, l.strand, l.gene_start0, l.gene_end0,
                              l.query, l.upstream_overlaps, l.downstream_overlaps), '\t'))
        end
    end
    return path
end

function write_loci_bed(path::AbstractString, loci::AbstractVector{Locus})
    open(path, "w") do io
        for l in loci
            println(io, join((l.chrom, l.start0, l.end0, l.id, 0, l.strand), '\t'))
        end
    end
    return path
end

"""
    read_loci(path) -> Vector{Locus}

Reads a loci TSV written by `write_loci`, or a BED file (`.bed`) whose intervals are
used as-is (no flanks added).
"""
function read_loci(path::AbstractString)
    lines = filter(l -> !isempty(strip(l)), readlines(path))
    loci = Locus[]
    if endswith(lowercase(path), ".bed")
        for l in lines
            (startswith(l, "#") || startswith(l, "track") || startswith(l, "browser")) && continue
            f = split(l, '\t')
            chrom, s0, e0 = String(f[1]), parse(Int, f[2]), parse(Int, f[3])
            name = length(f) >= 4 && !isempty(f[4]) ? f[4] : "$chrom:$s0-$e0"
            strand = length(f) >= 6 && f[6] == "-" ? '-' : '+'
            push!(loci, Locus(locus_id_from_gene(name), chrom, s0, e0, strand, s0, e0, name, "", ""))
        end
    else
        header = split(lines[1], '\t')
        header == LOCI_HEADER || error("$path: unrecognised loci header (expected a `unimap loci` TSV or a .bed)")
        for l in lines[2:end]
            f = split(l, '\t'; keepempty=true)
            push!(loci, Locus(f[1], f[2], parse(Int, f[3]), parse(Int, f[4]), f[5][1],
                              parse(Int, f[6]), parse(Int, f[7]), f[8], f[9], f[10]))
        end
    end
    return loci
end
