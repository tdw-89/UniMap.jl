# UniMap.jl

Empirical **paired-end mappability** for target loci, measured by simulation.

UniMap tiles each locus with error-free paired-end fragments over a grid of read
lengths (L) and fragment lengths (F). It aligns them to the **whole genome** with your
production aligner settings, then scores each position by how often the fragments
covering it map back to where they came from. It was built for comparing paralogs and
CNVs, where you need to know which segments of each copy can be measured on their own.

Unlike k-mer tools such as Umap and GenMap, UniMap models your actual aligner, MAPQ
filtering, read length and fragment-length distribution.

## Requirements

- Julia ≥ 1.12
- `bowtie2` and `samtools` on `PATH` (only needed by the alignment step)
- A reference FASTA with a samtools index (`samtools faidx ref.fa`), and a bowtie2
  index built **from that same FASTA**
- An annotation (GFF3 or GTF, optionally gzipped) or a BED file of loci

## Installation

```sh
git clone <repo> UniMap.jl && cd UniMap.jl
julia --project -e 'using Pkg; Pkg.instantiate()'
```

Run commands with:

```sh
julia -t 8 --project=/path/to/UniMap.jl -m UniMap <command> [options]
```

Or install a `unimap` executable to `~/.julia/bin` with `pkg> app develop .` from the
repo root (Julia ≥ 1.12). Then set `JULIA_NUM_THREADS` to control threading. The
examples below use `unimap` as shorthand.

## Quick start

```sh
# 1. Target loci: gene span + 1500 bp flanks (IDs can be gene IDs/names, transcript or protein IDs)
unimap loci -f ref.fa -a annotation.gff -g genes.txt -o loci.tsv

# 2. Simulate read pairs (gzipped FASTQ, one pair of files per locus and read length)
unimap simulate -f ref.fa -l loci.tsv -o sim/ -L 50,100,150 -F 100,150,200,250,300,400,500 \
    --step 5 --orientation plus

# 3. Write the alignment script, then run it (locally or inside a batch job)
unimap align-script -s sim/ -x ref_bt2_index -t 8 --aligner-args='--very-sensitive'
bash sim/align.sh

# 4. Score
unimap score -s sim/
```

Results are written to `sim/scores/`: per-base tracks in a Zarr store and summary tables
in Parquet. Every command takes `--help`.

## How it works

### Fragments and reads

A fragment of length F starting at position `s` gives two mates of length `m = min(L, F)`:

| orientation | R1 | R2 |
|---|---|---|
| `+` | first `m` bases (forward strand) | reverse complement of the last `m` bases |
| `-` | reverse complement of the last `m` bases | first `m` bases |

- **Orientation**: `--orientation both` simulates both `+` and `-` at every start
  position. `plus` or `minus` simulates only one, which halves the read count. The two
  orientations are the same two sequences with the mates swapped. On the zebrafish test
  loci, one-orientation tracks matched each other to a mean per-base difference of
  ≤ 0.0005 at MAPQ ≥ 2.
- When F < L, both mates are trimmed to F, which mirrors adapter trimming.
- All bases get a uniform quality (`I` by default, i.e. Q40).
- Fragments with a non-ACGT base in either mate are skipped and excluded from the
  denominator.
- Each read name records where the fragment came from:
  `locusID|chrom|fragStart0|F|orient`, with 0-based coordinates.

### Grid

Read lengths (`-L`) and fragment lengths (`-F`) are both explicit lists with no defaults.

- **Read lengths**: the lengths your runs round down to, usually their 10th–25th
  percentile post-trimming read length. Each read length is a separate bowtie2 run, so
  list only the ones you need.
- **Fragment lengths**: points spanning your runs' fragment-length histograms. They are
  later used to interpolate each run's histogram. The largest should equal the aligner's
  maximum fragment length (bowtie2 `-X`), because longer fragments can never align as
  proper pairs. All fragment lengths go into the same FASTQ files and bowtie2 run, so
  each one adds reads but not alignment runs.
- **Tiling step** (`--step`, default 1): bp between consecutive fragment starts.
  - Pairs simulated per base per read length = (number of fragment lengths) ×
    (orientations) / step.
  - The fragments covering a base = F × (orientations) / step.
  - A larger step coarsens the per-base average but doesn't bias it. Keep it well below
    the shortest read length (≤ L/5 is a safe rule), so every mate-length stretch of
    sequence is still sampled several times.

### Loci

- `unimap loci` builds each gene's span from its outermost exons across all transcripts.
- It adds **strand-aware** flanks: `--upstream` before the TSS and `--downstream` after
  the TES. For minus-strand genes, upstream is at higher coordinates.
- Flanks are clipped at chromosome ends.
- Flanks overlapping other genes are **flagged**, not trimmed. See the
  `upstream_flank_overlaps` and `downstream_flank_overlaps` columns.
- IDs that aren't found, or that match more than one gene, are reported and skipped.
- A BED file (`-l loci.bed`) is used as-is, with no flanks added.
- Only positions inside the locus are scored. UniMap simulates every fragment that
  overlaps the locus, so edge positions get full coverage.

### Scoring

For position x, read length L, fragment length F and MAPQ cutoff c:

$$m_c(x, L, F) = \frac{\#\{\text{fragments covering } x \text{ that map correctly with fragment MAPQ} \ge c\}}{\#\{\text{simulated fragments covering } x\}}$$

A fragment covers every base from its leftmost to its rightmost position, which is the
same span MACS3 `BAMPE` piles up. Each simulated pair is classified as:

- **correct**: a primary proper pair on the origin chromosome. Its leftmost position and
  fragment length must be within `--tolerance` of the origin (default 0), and R1 must
  be on the expected strand.
- **wrong locus**: a proper pair placed somewhere else. These feed the bleed matrix and
  the misplaced-sites table.
- **not proper**: unaligned, discordant, or only one mate aligned.

**Fragment MAPQ** is the lower MAPQ of the two mates (`--frag-mapq min`, the default; `max`
is also available). Use the same rule when counting real reads.

**Choosing a cutoff**: in `misplacement_summary`, find the lowest cutoff where
`frac_wrong_kept` is below about 1%. In bowtie2, MAPQ 0–1 means the pair has multiple
equally good hits.

## Commands

### `unimap loci`
| option | default | |
|---|---|---|
| `-f/--fasta` | required | reference FASTA (`.fai` needed) |
| `-a/--annotation` | required | GFF3 or GTF (`.gz` ok) |
| `-g/--genes` | required | file with one ID per line, or comma-separated IDs |
| `--upstream`, `--downstream` | 1500, 1500 | flank sizes (bp) |
| `-o/--out` | `loci.tsv` | |

### `unimap simulate`
| option | default | |
|---|---|---|
| `-f/--fasta` | required | must be the FASTA the aligner index was built from |
| `-l/--loci` | | `loci.tsv` or `.bed`; or use `-a` and `-g` (plus flank options) instead |
| `-o/--out` | required | output directory |
| `-L/--read-lengths` | required | comma-separated, e.g. `50,100,150` |
| `-F/--fragment-lengths` | required | comma-separated, e.g. `100,150,200,300,500` |
| `--step` | 1 | bp between fragment starts |
| `--orientation` | `both` | `both`, `plus` or `minus` |
| `--qual` | `I` | |
| `--no-checksum` | | skip the SHA-256 of the FASTA |

Threads over loci (`julia -t N`).

### `unimap align-script`
| option | default | |
|---|---|---|
| `-s/--sim` | required | `simulate` output directory |
| `-x/--index` | required | bowtie2 index prefix |
| `-t/--threads` | 4 | |
| `--aligner-args` | `''` | production bowtie2 arguments. Use the `=` form: `--aligner-args='--very-sensitive -X 500'`. `-X` is set to the largest fragment length if missing. |
| `--bowtie2`, `--samtools` | on `PATH` | executable paths |

Writes `SIM/align.sh`, which runs one bowtie2 job per read length to produce
`SIM/aln/sim_L{L}.bam`. The script is plain bash, so it can also run inside an `sbatch`
job. Re-running it skips read lengths that have already finished. **Don't sort or filter
the BAMs**: scoring needs the two mates of each pair next to each other.

### `unimap score`
| option | default | |
|---|---|---|
| `-s/--sim` | required | `simulate` output directory |
| `--bam-dir` | `SIM/aln` | |
| `-o/--out` | `SIM/scores` | |
| `--cutoffs` | `0,2,5,10,20,30` | fragment MAPQ cutoffs |
| `--tolerance` | 0 | bp |
| `--frag-mapq` | `min` | `min` or `max` |
| `-L/--read-lengths` | all | score a subset |
| `-f/--fasta` | from `simulate.toml` | use if the FASTA has moved |

Threads over read lengths. Denominators are rebuilt from the reference and checked against
`tiling.tsv`, so pairs that are missing from the BAM still count against the score.

## Outputs

**`simulate` output directory**
| file | contents |
|---|---|
| `fastq/L{L}/{locus}_R{1,2}.fq.gz` | simulated pairs, all F values |
| `loci.tsv`, `loci.bed` | loci used (0-based, half-open) |
| `tiling.tsv` | per (locus, L, F): step, start range, starts, N-skipped starts, pairs written |
| `simulate.toml` | reference path and SHA-256, grid, UniMap version and git commit |
| `align.sh`, `align.toml`, `aln/` | written by `align-script` and by running the script |

**`scores/`**
| file | contents |
|---|---|
| `mappability.zarr` | per-base tracks for every locus, concatenated end to end (described below) |
| `misplacement_summary.parquet` | one row per locus × L × F × cutoff, plus genome-wide `ALL` rows: `n_fragments, n_kept, n_correct_kept, n_wrong_kept, frac_wrong_kept, frac_correct, n_not_proper, n_missing, n_skipped_N, min_denominator` |
| `bleed_matrix.parquet` | one row per origin × destination × L × F × cutoff: `origin_locus`, `dest_locus` (another target locus, or `elsewhere`), `n_fragments`. Not written if nothing was misplaced. |
| `misplaced_sites.parquet` | one row per origin × 1 kb destination bin × fragment MAPQ: where wrong-locus pairs landed. Not written if nothing was misplaced. |
| `score.toml` | reference checksum, bowtie2 version and arguments, scoring settings |

### `mappability.zarr`

A Zarr v2 group. All loci are concatenated along one `position` axis, and each
(L, F, cutoff) combination is a channel.

| array | Julia dims | Python/xarray dims | type |
|---|---|---|---|
| `score` | `[position, F, mapq_cutoff, L]` | `(L, mapq_cutoff, F, position)` | Float32, `NaN` where nothing was simulated |
| `denominator` | `[position, F, L]` | `(L, F, position)` | UInt16, simulated fragments covering the base |
| `L`, `F`, `mapq_cutoff` | 1-D | 1-D | grid values |

- The group attributes hold `loci`: one entry per locus with `locus_id, chrom, start0,
  end0, strand, offset`. Locus k covers positions `offset+1 : offset+(end0-start0)` in
  Julia (`offset : offset+(end0-start0)` in Python). Position `i` within a locus is
  genomic `start0 + i - 1`.
- The attributes also hold `genome_checksum`, the grid and the scoring settings.
- Chunks are 16,384 positions × all F × all cutoffs for one L, compressed with
  Blosc/zstd. Reading one locus for one read length therefore touches only a few chunks.

Reading the tracks:

```julia
using UniMap
t = locus_track("sim/scores/mappability.zarr", "pimr158"; L=100, F=240, mapq_cutoff=2)
t.chrom, t.start0, t.score, t.denominator   # per-base vectors

using Zarr                                   # or the whole array
z = zopen("sim/scores/mappability.zarr"; consolidated=true)
z["score"]                                   # [position, F, mapq_cutoff, L]
```
```python
import xarray as xr
ds = xr.open_zarr("sim/scores/mappability.zarr")   # ds.score.sel(L=100, F=240, mapq_cutoff=2)
```

The summary tables are Parquet:

```julia
using Parquet, DataFrames
summary = DataFrame(read_parquet("sim/scores/misplacement_summary.parquet"))
```

## Scale

On 14 zebrafish loci of about 6–27 kb each (gene plus flanks), with 17 fragment lengths
(50–500 bp), both orientations and steps of 1–5 bp:

- Simulation: about 27 pairs per base per read length, about 400k pairs per locus per
  read length, about 7 MB of FASTQ per locus per read length, and about 13 s for 3 read
  lengths on 8 threads.
- Alignment dominates the runtime: bowtie2 took 7–14 min per read length (about 5.6M
  pairs) on 8 threads against the whole genome.
- Scoring all three read lengths took 17 s, and `mappability.zarr` is 12 MB.
- The quick-start settings (7 fragment lengths, `--step 5`, one orientation) give
  1.4 pairs per base per read length, about 20× fewer reads.

## Tests

```sh
julia --project -e 'using Pkg; Pkg.test(julia_args=["-t", "4"])'
```

The unit tests use a synthetic genome. If `bowtie2`, `bowtie2-build` and `samtools` are
on `PATH`, an end-to-end test also simulates, aligns and scores unique, exactly
duplicated and near-duplicated loci.

## Roadmap

- Plots of mappability along a locus
- Per-run weighting by observed fragment-length histograms and read-length quantiles
- SNP/error sensitivity runs, calling of unique segments with RepeatMasker filtering, and
  ATAC and bisulfite modes
