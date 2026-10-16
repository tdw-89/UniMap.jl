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
unimap simulate -f ref.fa -l loci.tsv -o sim/ -L 50,100,150 --fmax 500

# 3. Write the alignment script, then run it (locally or inside a batch job)
unimap align-script -s sim/ -x ref_bt2_index -t 8 --aligner-args='--very-sensitive'
bash sim/align.sh

# 4. Score
unimap score -s sim/
```

Results are written to `sim/scores/` as Parquet files. Every command takes `--help`.

## How it works

### Fragments and reads

A fragment of length F starting at position `s` gives two mates of length `m = min(L, F)`:

| orientation | R1 | R2 |
|---|---|---|
| `+` | first `m` bases (forward strand) | reverse complement of the last `m` bases |
| `-` | reverse complement of the last `m` bases | first `m` bases |

- **Both orientations** are simulated at every start position.
- When F < L, both mates are trimmed to F, which mirrors adapter trimming.
- All bases get a uniform quality (`I` by default, i.e. Q40).
- Fragments with a non-ACGT base in either mate are skipped and excluded from the
  denominator.
- Each read name records where the fragment came from:
  `locusID|chrom|fragStart0|F|orient`, with 0-based coordinates.

### Grid

- **Read lengths**: `-L`, default `36,50,75,100,150`.
- **Fragment lengths**: start at `--fmin` (default 50). The step starts at 1 and doubles
  every two steps, and `--fmax` is always included:
  `50, 51, 52, 54, 56, 60, 64, 72, 80, 96, 112, 144, 176, 240, 304, 432, 500`.
  **Set `--fmax` to the aligner's maximum fragment length (bowtie2 `-X`).** Longer
  fragments can never align as proper pairs.
- **Tiling step**: `--step-policy`, default `150:1,300:2,*:5`. This means step 1 for
  F ≤ 150, step 2 for F ≤ 300, and step 5 above that. A larger step lowers the
  denominator but not the expected score.

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
| `-L/--read-lengths` | `36,50,75,100,150` | |
| `--fmin`, `--fmax` | 50, 500 | |
| `--step-policy` | `150:1,300:2,*:5` | |
| `--qual` | `I` | |
| `--no-checksum` | | skip the SHA-256 of the FASTA |

Threads over loci (`julia -t N`).

### `unimap align-script`
| option | default | |
|---|---|---|
| `-s/--sim` | required | `simulate` output directory |
| `-x/--index` | required | bowtie2 index prefix |
| `-t/--threads` | 4 | |
| `--aligner-args` | `''` | production bowtie2 arguments. Use the `=` form: `--aligner-args='--very-sensitive -X 500'`. `-X` is set to `fmax` if missing. |
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
| file | one row per | key columns |
|---|---|---|
| `mappability_tracks.parquet` | run of constant score | `locus_id, chrom, start0, end0, L, F, mapq_cutoff, score, min_denominator, genome_checksum` |
| `misplacement_summary.parquet` | locus × L × F × cutoff, plus genome-wide `ALL` rows | `n_fragments, n_kept, n_correct_kept, n_wrong_kept, frac_wrong_kept, frac_correct, n_not_proper, n_missing, n_skipped_N, min_denominator` |
| `bleed_matrix.parquet` | origin × destination × L × F × cutoff | `origin_locus, dest_locus` (another target locus, or `elsewhere`), `n_fragments` |
| `misplaced_sites.parquet` | origin × 1 kb destination bin × fragment MAPQ | where wrong-locus pairs landed |
| `score.toml` | | reference checksum, bowtie2 version and arguments, scoring settings |

Reading the tables:

```julia
using Parquet, DataFrames
tracks = DataFrame(read_parquet("sim/scores/mappability_tracks.parquet"))
```
```python
import pandas as pd
tracks = pd.read_parquet("sim/scores/mappability_tracks.parquet")
```

## Scale

On 14 zebrafish loci of about 6–27 kb each (gene plus flanks), with the default F grid:

- Simulation: about 400k pairs per locus per read length, about 7 MB of FASTQ per locus
  per read length, and about 13 s for 3 read lengths on 8 threads.
- Alignment dominates the runtime: bowtie2 handled about 2.1M pairs in 2.3 min on 8
  threads against the whole genome.

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
