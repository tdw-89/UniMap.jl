struct FaidxEntry
    length::Int
    offset::Int
    linebases::Int
    linewidth::Int
end

"""
    IndexedFasta(path)

Random access to a FASTA file via its samtools `.fai` index.
"""
struct IndexedFasta
    path::String
    entries::Dict{String,FaidxEntry}
    names::Vector{String}
end

function IndexedFasta(path::AbstractString)
    fai = path * ".fai"
    isfile(path) || error("FASTA not found: $path")
    isfile(fai) || error("FASTA index not found: $fai (create it with `samtools faidx $path`)")
    entries = Dict{String,FaidxEntry}()
    names = String[]
    for line in eachline(fai)
        isempty(line) && continue
        f = split(line, '\t')
        entries[f[1]] = FaidxEntry(parse(Int, f[2]), parse(Int, f[3]), parse(Int, f[4]), parse(Int, f[5]))
        push!(names, f[1])
    end
    return IndexedFasta(String(path), entries, names)
end

function faidx_entry(fa::IndexedFasta, chrom::AbstractString)
    e = get(fa.entries, chrom, nothing)
    e === nothing && error("sequence '$chrom' is not in $(fa.path).fai")
    return e
end

chrom_length(fa::IndexedFasta, chrom::AbstractString) = faidx_entry(fa, chrom).length
Base.haskey(fa::IndexedFasta, chrom::AbstractString) = haskey(fa.entries, chrom)

"""
    fetch_seq(fa, chrom, start0, end0) -> Vector{UInt8}

Bases in the 0-based half-open interval `[start0, end0)`, uppercased.
"""
function fetch_seq(fa::IndexedFasta, chrom::AbstractString, start0::Integer, end0::Integer)
    e = faidx_entry(fa, chrom)
    0 <= start0 <= end0 <= e.length ||
        throw(ArgumentError("interval [$start0, $end0) is outside $chrom (length $(e.length))"))
    n = end0 - start0
    out = Vector{UInt8}(undef, n)
    n == 0 && return out
    byteoffset(p) = e.offset + div(p, e.linebases) * e.linewidth + rem(p, e.linebases)
    first, last = byteoffset(start0), byteoffset(end0 - 1)
    raw = open(fa.path) do io
        seek(io, first)
        read(io, last - first + 1)
    end
    k = 0
    for b in raw
        (b == UInt8('\n') || b == UInt8('\r')) && continue
        k += 1
        out[k] = UInt8('a') <= b <= UInt8('z') ? b - 0x20 : b
    end
    k == n || error("FASTA layout disagrees with $(fa.path).fai for $chrom; re-run samtools faidx")
    return out
end

file_sha256(path::AbstractString) = bytes2hex(open(sha256, path))
