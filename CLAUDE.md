# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

UniMap.jl measures empirical paired-end mappability of target loci by simulation: it tiles loci
with error-free read pairs, aligns them genome-wide with bowtie2, and scores how often fragments
covering each base map back to their origin. `README.md` is the user-facing reference (CLI
options, output schemas, scoring definitions). `docs/probability_model.md` gives the formalism
behind the counts. `docs/simulation_findings.md` records measured parameter choices and costs.

## Commands

```sh
julia --project -e 'using Pkg; Pkg.instantiate()'                 # setup (Julia >= 1.12)
julia --project -e 'using Pkg; Pkg.test(julia_args=["-t", "4"])'  # full test suite
julia -t 4 --project test/runtests.jl                             # faster: run tests without Pkg.test
julia -t 8 --project -m UniMap <command> --help                   # run the CLI (loci|simulate|align-script|score|view)
```

- All tests are in `test/runtests.jl`, one top-level `@testset "UniMap"` with nested testsets, and
  there is no filter mechanism. To iterate on a single testset, run the file's setup lines (the
  helpers and `const` genome at the top) plus that block in a REPL with `using UniMap, Test, Random`.
- The `"align + score (synthetic genome)"` end-to-end testset only runs when `bowtie2`,
  `bowtie2-build` and `samtools` are on `PATH`. Otherwise it is skipped with a warning.
- Use `-t N` with N > 1 so the threaded paths get exercised: `simulate` threads over loci and
  `score` threads over BAM shards, regions and read lengths.

## Architecture

The pipeline has four stages, each a CLI subcommand in `src/cli.jl`, and they communicate **only
through files in the simulation directory**:

1. **`loci.jl`**: GFF3/GTF/BED → `Locus` records (gene span + strand-aware flanks).
   `merge_regions` merges loci that could share a fragment (gap ≤ max F − 2) into *simulation
   regions*. Downstream code tiles and counts regions, then cuts per-locus outputs back out
   (`locus_rows`).
2. **`tiler.jl`** (`simulate`): writes FASTQ per region × read length, plus `loci.tsv`,
   `regions.tsv`, `tiling.tsv` and `simulate.toml` (grid + reference SHA-256). The fragment's
   origin is encoded in the read name (`read_name`/`parse_read_name`:
   `regionID|chrom|fragStart0|F|orient`). This is the only link between a BAM record and where
   it was simulated from.
3. **`align.jl`** (`align-script`): writes `align.sh` (one bowtie2 run per L → unsorted
   `aln/sim_L{L}.bam`) and `align.toml`. UniMap never runs bowtie2 itself.
4. **`score.jl`** (`score_simulation`): does **not** read the FASTQs. It re-enumerates every
   simulated fragment from the reference with `fragment_starts`/`fragment_has_n` into
   `LocusFrags` (a flat status/MAPQ vector indexed by F, then start, then orientation; see
   `frag_slot`) and checks the result against `tiling.tsv`. It then streams the name-grouped BAM
   (`scan_bam!`), classifies each pair (`classify_pair` → `ST_*` codes), and builds per-base
   counts (`region_tracks!`) plus Parquet summaries.

Invariants that span files:

- `simulate` and `score` must enumerate fragments identically. Any change to `fragment_starts`,
  `sim_region`, `fragment_has_n`, orientation order or step handling affects both, and
  `tiling.tsv` validation will fail on old simulation directories.
- Coordinates are 0-based half-open everywhere (`start0`/`end0`) except Julia array indexing.
  Zarr locus `k` occupies positions `offset+1 : offset+(end0-start0)`.
- Tracks store raw **counts** (`n_simulated`, `n_correct`, `n_ectopic`), never ratios. Derived
  quantities (`map_back`, `ectopic_rate`, `log2_ratio` with `LOG_RATIO_PSEUDOCOUNT`) are computed
  on read, in `locus_track` and in the viewer.
- BAMs must stay in bowtie2 output order (mates adjacent). Never sort or filter them. A read
  length's BAM may be split into shards `sim_L{L}_*.bam` (`bam_paths`); they hold disjoint
  pairs and are scanned concurrently into the same `LocusFrags`, with per-shard
  landings/bleed/sites combined by `merge_scans` (which also catches a pair in two shards).
- Zarr writes go through one lock: Zarr.jl's Blosc codec calls `blosc_compress` and
  `blosc_set_compressor`, which use Blosc's global state. Don't remove it; parallelism inside a
  chunk comes from `Blosc.set_num_threads`.
- `n_simulated`/`n_correct` are UInt16, so 2×max F must fit in UInt16 (checked in `score_simulation`).

**Viewer** (`viewer.jl`, a submodule `UniMap.Viewer`): GLMakie is deliberately **not** loaded at
module load. `load_viewer()` in `src/UniMap.jl` `include`s `viewer.jl` on first use and calls it
through `invokelatest`, so the CLI and HPC runs need no display or GLMakie precompile. Keep
viewer-only dependencies out of the top-level `using` list. The viewer reaches the parent module
through `const U = parentmodule(@__MODULE__)`.

**Provenance**: `simulate.toml` and `score.toml` record the UniMap version, git commit, reference
checksum and aligner settings. Keep these updated when adding settings that affect results.
