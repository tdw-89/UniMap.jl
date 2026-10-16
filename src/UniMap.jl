module UniMap

using ArgParse
using CodecZlib
using Dates
using Parquet: write_parquet
using SHA
using TOML
using XAM: BAM

export fragment_length_grid, StepPolicy, parse_step_policy, tiling_step,
       IndexedFasta, fetch_seq, chrom_length,
       Locus, Annotation, read_annotation, build_loci, read_loci, write_loci, read_id_list,
       fragment_starts, fragment_has_n, read_name, parse_read_name, TileSummary, simulate,
       write_align_script, score_simulation

include("grid.jl")
include("reference.jl")
include("loci.jl")
include("tiler.jl")
include("align.jl")
include("score.jl")
include("cli.jl")

end
