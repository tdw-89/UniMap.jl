# Handoff: empirical paired-end mappability by simulated read tiling (Julia)

## Purpose

I'm comparing broad histone-mark ChIP-seq signal (H3K9me2, H3K9me3, H3K27me3, with input DNA controls) between **paralogous genes**, in multiple organisms. Recent paralogs are highly similar, so reads in shared sequence can't be assigned to one copy reliably. The plan is to compare paralogs **only over segments where reads map uniquely and back to where they came from** (unique introns, promoters, flanks). Broad marks blanket the whole locus, so a unique sub-segment can stand in for the gene.

This module estimates **empirical mappability** for every position in a set of target loci. It works by tiling the loci with simulated, error-free paired-end fragments, aligning them to the **whole genome** with the same aligner and settings used for the real data, and scoring how often the fragments covering each position map back correctly.

Existing tools don't cover this. Umap/Bismap and GenMap are single-end, k-mer based, and don't reproduce my aligner, MAPQ, or fragment-length behaviour. Read simulators (wgsim, Mason, ART) place fragments at random rather than tiling. So: **write the tiling and scoring in Julia; use bowtie2 for alignment.**

---

## Decisions already made (don't revisit without asking)

1. **Paired-end simulation, not single reads.** A fragment of length F yields two mates of read length L: R1 = first L bases, R2 = reverse complement of the last L bases. When F < L, trim both mates to F (this mirrors adapter trimming by Trim Galore in the real pipeline).
2. **Systematic tiling, not random sampling.** Every fragment start position (with a configurable step), on **both orientations**: R1 on the forward strand, and R1 on the reverse strand (i.e., R1 = revcomp of the fragment's 3' end).
3. **Error-free reads with uniform high quality** (e.g., all `I`, Q40). Bowtie2's mismatch penalty depends on quality, so quality has to be constant for the baseline to be well defined. Adding errors or SNPs is a later, optional sensitivity run (see Phase 4).
4. **Align to the full reference genome**, never just the locus or the paralog pair. A read can match a third copy, a pseudogene, or a repeat anywhere in the genome.
5. **Use the same aligner settings as production.** These should be configurable by the user
6. **Two-dimensional grid: read length L × fragment length F.**
   - L grid (configurable), default `[36, 50, 75, 100, 150]`.
   - F grid: geometric-ish. Start at `F_min`, step starts at 1 and **doubles every two steps**, and `F_max` is always included. Example from 50: `50, 51, 52, 54, 56, 60, 64, 72, 80, 96, 112, 144, 176, …, F_max`.
   - **`F_max` = the aligner's maximum fragment length** (`-X`; 500 for the ChIP data). Fragments longer than that can't align as proper pairs, in simulation or real data.
7. **Target loci**: each gene = union of its transcripts (outermost TSS → outermost TES) **plus strand-aware flanks**: 1500 bp upstream of the TSS and 1500 bp downstream of the TES. For minus-strand genes, upstream is at higher genomic coordinates. Flag (don't drop) flanks that overlap neighbouring genes.
8. **Padding for simulation: extend each locus by `F_max` on each side.** Scoring uses **fragment-span coverage** (a fragment covers every base from its leftmost to its rightmost position, the usual paired-end ChIP pileup and what MACS3 BAMPE does). A fragment covering position x can start up to `F_max − 1` bp upstream, so positions near the edges only get full coverage with full-`F_max` padding. Only score positions inside the unpadded locus.
9. **Scoring metric** for position x, read length L, fragment length F, MAPQ cutoff c:
   `m_c(x, L, F) = (# simulated fragments covering x that mapped correctly and pass cutoff c) / (# simulated fragments covering x)`.
10. **Choose the MAPQ cutoff after the fact, from simulation results.** Compute scores for a set of cutoffs, default `[0, 2, 5, 10, 20, 30]`. In bowtie2, MAPQ 0–1 means equally good multiple hits. Reads with a second-best hit one mismatch worse get low but nonzero MAPQ, which is how bleed-through between paralogs happens.
11. **Per-run weighting by fragment-length distribution** (Phase 3). Each real sequencing run gets an expected mappability track by weighting the F-grid results with that run's observed fragment-length histogram. Interpolating along F is fine because the result is averaged over a distribution. **Read length is floored to the grid**, using a low quantile (10th–25th percentile) of the run's post-trimming read length, so estimates err on the conservative side.

---

## Inputs

- **Reference FASTA** (+ `.fai`), the exact one the bowtie2 index was built from, plus the **bowtie2 index** prefix. Record a checksum of the FASTA in outputs.
- **Annotation** (GTF/GFF3) and a **list of target gene IDs**, or a pre-built BED of loci. Support both if it's cheap.
- **Grid config**: L list, `F_min`, `F_max`, tiling step policy, MAPQ cutoff list, flank sizes (default 1500/1500).
- **(Phase 3) run table**: one row per SRA run with assay type, aligner (+ version/args), reference, bottom-quartile read length (per mate for PE), paired-end flag, and a path to its fragment-length histogram TSV (`fragment_length<TAB>count`, from a samtools/awk SLURM job I already ran).

## Tiling step policy

Mappability changes slowly along the sequence for long fragments, so use a coarser step for large F to save compute. Make it configurable. Default: step 1 for F ≤ 150, 2 for F ≤ 300, 5 above that. The score is a ratio, so the step only changes coverage depth, not the expected value, as long as every position gets enough fragments. Report the minimum denominator per locus.

---

## Phase 1: Tiler

For each (locus, L) produce one **pair of gzipped FASTQ files containing all F values**. That means one bowtie2 run per (genome, L), so the index is loaded once.

- Coordinates: 0-based half-open internally; be explicit wherever converting.
- **Read names must encode the origin**, with no whitespace (bowtie2 truncates names at whitespace and strips `/1` `/2`). Suggested: `locusID|chrom|fragStart0|F|orient`, where `orient ∈ {+,-}` (`+`: R1 from the forward strand at the fragment start; `-`: R1 = revcomp of the fragment end). L is implied by the file. Check that QNAME stays under 254 characters.
- Skip fragments containing `N` in the sequenced bases. Keep a per-position count of skipped fragments, and leave them out of the denominator.
- Stream output. Don't hold all reads in memory.
- Thread over loci if that's easy. Generating reads is cheap; alignment dominates the cost.

## Phase 2: Align, then score

**Alignment** runs outside Julia, as a SLURM job on the Unity HPC. I'm a non-admin user who cares about fair-share, so keep resource requests modest. Write an sbatch template or command generator that runs:

```
bowtie2 -x <index> -1 R1.fq.gz -2 R2.fq.gz --threads N [production args; default: none, i.e. -X 500] \
  | samtools view -b -o sim_L{L}.bam -
```

Leave it unsorted. Read names stay adjacent, which keeps pairing simple. Don't filter.

**Parsing**: XAM.jl

**Correctness classification per simulated fragment** (pair):
- `correct`: primary, both mates mapped, properly paired (flag 0x2), same chromosome as origin, leftmost aligned position == origin leftmost (allow a tolerance parameter, default 0; reads are error-free and end-to-end), fragment length == F, and the expected orientation.
- `wrong_locus`: properly paired but placed elsewhere. **Record where it landed** (chrom, pos) and, if it overlaps another target locus, which one. This builds a paralog "bleed matrix" saying which copy loses reads to which.
- `not_proper`: unaligned, discordant, or mixed-mode single-mate alignments. These don't count toward the real pileups either, which use proper pairs only.
- **Fragment MAPQ = min(MAPQ of R1, MAPQ of R2).** This is conservative. Make it a named option, and note that real-data counting has to use the same rule.

**Coverage accumulation**: use difference arrays per locus. For each fragment, add +1 at origin start and −1 at origin end, then take a prefix sum. Keep separate arrays for the denominator (all simulated fragments) and for numerators per MAPQ cutoff (correct and fragMAPQ ≥ c). This is O(fragments + positions).

**Outputs (Parquet)**:
1. `mappability_tracks`: run-length encoded intervals of constant score: `genome_checksum, locus_id, chrom, start0, end0, L, F, mapq_cutoff, score, min_denominator`. Only positions inside the unpadded locus.
2. `misplacement_summary`: per (locus, L, F, cutoff): fragments kept, fraction kept that are `wrong_locus`. Plus a genome-wide row per (L, F, cutoff). This is what I'll use to choose the cutoff: the lowest cutoff where misplacement among kept fragments is below ~1%.
3. `bleed_matrix`: (origin locus → destination locus or "elsewhere") counts per (L, F, cutoff).
4. A small metadata record: reference checksum, bowtie2 version and full args, grid config, git commit of this code.

## Phase 3: Per-run weighting

For each run in the run table:
- Pick the L grid point by **floor** of the run's low-quantile read length (for PE, use the shorter mate).
- Truncate the run's fragment histogram at `F_max`, then renormalize.
- Spread each histogram length's mass onto the two neighbouring F grid points with linear interpolation weights. This gives weights `w(F)` that sum to 1.
- Output `m_run,c(x) = Σ_F w(F) · m_c(x, L_run, F)` for each cutoff, as RLE Parquet keyed by run ID.
- Deduplicate: many runs share (reference, aligner+args, L grid point, layout). Compute the per-(L, F) tracks once per such "sim spec" and only redo the weighting per run.

## Phase 4 (optional, later)

- **Polymorphism/error sensitivity**: at one typical (L, F), add substitutions at a configurable rate. Zebrafish are highly polymorphic and lab strains differ from the reference, so SNPs matter more than sequencing errors. Re-score, and report how much the unique segments shrink, especially for young paralogs.
- **Unique-segment calling**: threshold a run's track (e.g., ≥ 0.9) with a minimum segment length, then intersect with a repeat annotation (RepeatMasker). Exclude or flag TE-derived segments, since TEs attract H3K9me3 independently of the gene.
- **ATAC**: same pipeline, with ATAC's own aligner args and `-X`.
- **Bisulfite**: needs in-silico conversion (C→T / G→A per strand) and the bisulfite aligner (Bismark or bwa-meth, whichever produced the data). Design the tiler so a read-transform step can be plugged in, but don't implement it in v1.

---

## Validation (write these as tests)

Build a small synthetic genome (random sequence, a few hundred kb):
1. **Unique random locus** → score ≈ 1 everywhere at every cutoff.
2. **Exact duplicate segment** (e.g., 2 kb copied elsewhere) → score ≈ 0 inside the duplicated interior at cutoff ≥ 2, rising towards its edges as F increases (fragments that span into unique flank get anchored).
3. **Near-duplicate** (copy with a SNP every ~200 bp) → intermediate scores, nonzero `wrong_locus` at low cutoffs. Check the bleed matrix shows origin → copy.
4. **Edge coverage**: the denominator is constant across the whole unpadded locus (confirms the padding is right).
5. **F < L trimming** produces reads of length F, not L.
6. **Orientation round-trip**: `+` and `-` fragments from the same position both classify as correct in unique sequence.
7. **Weighting**: weights sum to 1; a histogram concentrated at one grid point reproduces that grid point's track exactly.

Sanity checks on real data (report, don't assert):
- Mappability should rarely *decrease* with longer L or larger F. Report the fraction of positions where it does.
- Compare against the precomputed **Umap** zebrafish track (https://bismap.hoffmanlab.org/) at the nearest k. Umap is single-end exact matching, so wherever Umap says unique, PE simulation should almost always agree. Note the assembly may differ from my reference.
- Real input-DNA reads in segments the simulation calls unique should be mostly high-MAPQ, proper-pair.