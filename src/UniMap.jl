module UniMap

using ArgParse
import Blosc
using CodecZlib
using Dates
using Parquet: write_parquet
using SHA
using TOML
using XAM: BAM
using Zarr

export parse_int_list, parse_orientation,
       IndexedFasta, fetch_seq, chrom_length,
       Locus, Annotation, read_annotation, build_loci, read_loci, write_loci, read_id_list,
       fragment_starts, fragment_has_n, read_name, parse_read_name, TileSummary, simulate,
       write_align_script, score_simulation, locus_track, map_back, ectopic_rate, log2_ratio,
       view_mappability

include("grid.jl")
include("reference.jl")
include("loci.jl")
include("tiler.jl")
include("align.jl")
include("score.jl")
include("cli.jl")

# GLMakie is loaded only when the viewer is used, so the CLI and HPC runs never pay for it.
function load_viewer()
    isdefined(@__MODULE__, :Viewer) || Base.include(@__MODULE__, joinpath(@__DIR__, "viewer.jl"))
    return invokelatest(getfield, @__MODULE__, :Viewer)
end

"""
    view_mappability(zarr_path; plot=:base, mapq_cutoff=2, read_lengths=nothing, locus=nothing, wait=true)

Open an interactive GLMakie window for `scores/mappability.zarr`, with a locus dropdown
(and prev/next buttons) and a metric dropdown (Pr(map-back), ex-map rate, log2 ratio).

- `plot=:base`: one 3-D plot per read length of the metric over (fragment length x MAPQ
  cutoff) at one base, with a slider over the base.
- `plot=:track`: one 3-D plot per read length of the metric over (position along the locus
  x fragment length) at a fixed MAPQ cutoff (`mapq_cutoff`, changeable in the window).
- `plot=:boxplot`: per read length, box plots across the locus's bases of the metric for
  each fragment length, at a fixed MAPQ cutoff.

With `wait=false` (e.g. from the REPL), the window opens and the viewer is returned instead
of blocking until the window closes.
"""
function view_mappability(path::AbstractString; wait::Bool=true, kw...)
    V = load_viewer()
    return Base.invokelatest(wait ? V.run_viewer : V.open_viewer, path; kw...)
end

end
