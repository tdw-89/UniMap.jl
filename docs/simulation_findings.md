# Simulation parameter findings

Empirical results from tuning the UniMap simulation grid on zebrafish (GRCz12tu) paralogs.
Written so a person or an LLM can pick up the work without the original chat history.
`README.md` covers motivation, the scoring model and usage. This file records what was
measured, how, and what it implies for parameter choices.

Last updated: 2026-10-18.

## TL;DR: recommended settings

| parameter | recommendation | why (section) |
|---|---|---|
| fragment MAPQ cutoff | **>= 2** (fragment MAPQ = min of the two mates) | removes bowtie2 ties; >= 10 costs ~11 points of map-back for little gain (3) |
| orientation | `--orientation plus` | same result as `both` at MAPQ >= 2, half the reads (2) |
| read lengths | 50, 100, 150 (or whatever the real runs round down to) | each L is a separate bowtie2 run; L=50 is clearly worse (3) |
| fragment lengths | 50,75,100,125,150,175,200,250,300,350,400,450,500 | 25 bp spacing over 100-200 is needed; only 400 can be dropped cheaply (4) |
| `--step` | **10** (15 at most) | halves the reads vs step 5 with small error; >= 20 blurs per-base structure (5) |

Expected full-panel cost with these settings is in section 7.

## 1. Context

- **Goal.** Per-base empirical paired-end mappability for paralog loci, to be combined later with
  each real run's fragment-length histogram and read length into confidence formulas. The user
  wants simple count-based quantities (Pr(map-back), ex-map rate, log2 ratio; see README) and
  avoids loose terms such as "precision"/"recall".
- **Genome and tools.** GRCz12tu RefSeq (`GCF_049306965.2`), bowtie2 2.5.5 end-to-end with
  default settings plus `-X <max F>`, samtools. Index prefix `test/test_data/GRCz12tu`.
- **Gene universe.** `test/test_data/ds_zebrafish.csv` (paralog pairs with dS), mapped to gene
  IDs with `protein_to_gene.csv` and the GFF.
- **Real fragment sizes.** `test/test_data/stats_all_bams.txt` has fragment-length quantiles
  (q10..q90) for ChIP, ATAC and bisulfite runs. Only the ChIP rows (`/ChIP/` in the path) were
  used. "Short", "typical" and "long" runs are the 25th/50th/75th percentile across ChIP runs of
  each quantile:

  | run | q10 | q25 | q50 | q75 | q90 |
  |---|---|---|---|---|---|
  | short | 96 | 120 | 154 | 196 | 242 |
  | typical | 126 | 155 | 196 | 255 | 323 |
  | long | 146 | 177 | 224 | 289 | 361 |

  Each run's per-bp density (piecewise-linear CDF through the quantiles, from q10/2 up to the
  largest F) is projected onto the F grid by linear interpolation. Weights on the 13-point grid
  for the typical run: 50:0.004, 75:0.034, 100:0.040, 125:0.079, 150:0.136, 175:0.150,
  200:0.177, 250:0.173, 300:0.099, 350:0.037, 400:0.028, 450:0.028, 500:0.014. F < 100
  carries 0-5% of the weight for ChIP; ATAC would weight short fragments more.
- **Run-weighted track.** For a base x, read length L and cutoff c:
  `T(x) = sum_F w_F * n_correct[x,F,c,L] / n_simulated[x,F,L]`. Most comparisons below are
  on this track, and on whether it is >= 0.9 ("0.9 flips" = share of bases whose call changes).
- **All simulated reads are error-free** with uniform quality. Every number here is therefore a
  best case for misplacement.

## 2. Early findings on `sim_quick` (14 loci)

`test/test_data/sim_quick`: 14 hand-picked paralogs (`test_loci.tsv`),
`-L 50,100,150 -F 50,100,150,200,250,300,400,500 --step 5 --orientation plus`.

- **bowtie2 ties.** Error-free reads with an equally good hit elsewhere get MAPQ 0 or 1, and
  bowtie2 picks one of the tied hits at random. Essentially all misplacement sits at MAPQ 0-1;
  at MAPQ >= 2 there were zero ectopic fragments in `sim_quick`.
- **Orientation.** `plus`, `minus` and `both` gave the same tracks at MAPQ >= 2, so simulating
  one orientation halves the reads at no cost.
- **Scale (earlier 14-locus runs, loci of 6-27 kb).** With 17 fragment lengths, both
  orientations and steps of 1-5 bp: ~5.6M pairs per read length, bowtie2 7-14 min per read
  length on 8 threads, scoring 17 s, Zarr 12 MB. With the `sim_quick` settings: 332k pairs per
  read length, 1.5 min total alignment, 13 s scoring.

## 3. MAPQ cutoff (random sub-panels)

**Setup** (`test/test_data/subsets/`, script `scripts/subsets.jl`, seed 2026):

- 3 sub-panels x 12 paralog pairs, stratified 4 pairs each by pair dS < 0.1, 0.1-0.5, >= 0.5.
  Genes <= 60 kb, no gene used twice: 72 genes (`genes.txt`, `pairs.tsv`), 72 loci with
  1.5 kb flanks, merged into 71 simulation regions.
- Simulated into `subsets/sim` with
  `-L 50,100,150 -F 50,75,100,125,150,175,200,250,300,350,400,450,500 --step 5 --orientation plus`:
  2,510,930 pairs per read length (~1 Mb of simulated sequence).
- Aligned with `align.sh -t 8` and scored with `--cutoffs 0,1,2,3,5,8,10,20,30,40`.
- Analysis: `scripts/analyze_subsets.jl`. Counts in `misplacement_summary.parquet` are weighted
  by the typical-run F weights.

**bowtie2 MAPQ plateaus.** In these data, results are identical within the groups {0, 1},
{2, 3, 5}, {8, 10, 20, 30} and {40}, so only one value per group matters.

**All three sub-panels pooled** (misplaced % = wrong / kept; map-back = correct / simulated;
last column = share of bases with T >= 0.9):

| L | MAPQ >= | misplaced % | Pr(map-back) | ectopic / kept | bases T >= 0.9 |
|---|---|---|---|---|---|
| 50 | 0 | 24.6 | 0.754 | 0.043 | 52.5% |
| 50 | 2 | 0.098 | 0.686 | 2e-5 | 48.0% |
| 50 | 10 | 0.041 | 0.605 | 1e-5 | 40.4% |
| 50 | 40 | 0.020 | 0.357 | 1e-5 | 21.8% |
| 100 | 0 | 17.4 | 0.826 | 0.034 | 68.0% |
| 100 | 2 | 0.021 | 0.773 | 1e-5 | 64.3% |
| 100 | 10 | 0.005 | 0.660 | <1e-5 | 50.4% |
| 100 | 40 | 0.004 | 0.363 | <1e-5 | 23.9% |
| 150 | 0 | 16.2 | 0.838 | 0.032 | 70.1% |
| 150 | 2 | 0.006 | 0.787 | <1e-5 | 66.8% |
| 150 | 10 | 0.001 | 0.668 | <1e-5 | 53.0% |
| 150 | 40 | 0.001 | 0.354 | <1e-5 | 23.7% |

- **MAPQ >= 2 removes the ties.** Misplacement falls from 16-25% to <= 0.1%, and ectopic
  fragments per kept fragment from 3-4% to ~1e-5.
- **The sub-panels agree.** At MAPQ >= 2 misplacement ranged 0.003-0.14% across panels and
  read lengths (worst: panel 3, L=50).
- **MAPQ >= 10** removes a further 0.015-0.06 points of misplacement but loses ~11 points of
  map-back and 8-14 points of bases with T >= 0.9.
- **MAPQ 40 is unusable**: map-back roughly halves.
- Map-back rises with read length; the step from 50 to 100 is much larger than 100 to 150.

## 4. Fragment-length grid

Greedy backward elimination from the 13-point grid at MAPQ >= 2, always keeping 50 and 500,
removing at each step the F whose loss changes T least (by p99). Errors are |T_reduced -
T_full|, worst case over the short/typical/long runs and the three read lengths.

| kept fragment lengths | max change | p99 change | 0.9 flips |
|---|---|---|---|
| all 13 | 0 | 0 | 0 |
| drop 400 (12) | 0.014 | 0.006 | 0.17% |
| drop 400, 75 (11) | 0.028 | 0.015 | 0.37% |
| drop 400, 75, 450 (10) | 0.033 | 0.015 | 0.38% |
| 50,100,125,150,200,250,300,350,500 (9) | 0.052 | 0.025 | 0.69% |
| 50,100,150,200,250,300,350,500 (8) | 0.068 | 0.030 | 1.39% |
| 50,100,150,200,250,500 (6) | 0.114 | 0.061 | 1.49% |
| 50,100,150,250,500 (5) | 0.237 | 0.090 | 2.45% |
| 50,150,500 (3) | 0.450 | 0.273 | 8.50% |

- Each F costs about the same to simulate and align, so dropping one F saves ~1/13 (~8%).
- Only dropping 400 is nearly free. The 8-point grid used in `sim_quick` is too coarse
  (1.4% of bases flip): most ChIP weight lies in 125-300 bp, which needs ~25 bp spacing.
- Keep 50 even though ChIP gives it little weight: it anchors the low end for short-fragment
  runs (and ATAC).
- Caveat: errors are measured against the 13-point grid itself, not a continuous truth, so
  they understate the true interpolation error of every grid including the full one.

## 5. Tiling step (`--step`)

bowtie2 time is roughly linear in the number of pairs, and pairs per base per read length =
(#F x #orientations) / step, so the step is the largest cost lever.

**Method: emulating coarser steps without re-aligning.** bowtie2 seeds its per-read
randomness from the read name, sequence and qualities, so a read aligns the same way whatever
else is in the run. A step-K tiling (K a multiple of 5) is the subset of step-5 fragments with
`(fragStart0 - first_start(region, F)) % K == 0`, where
`first_start = max(0, region.start0 - F + 1)` (see `fragment_starts` in `src/tiler.jl`).
`scripts/step_subsample.sh`:

1. runs `unimap simulate --step K` into `subsets/stepK` (for `simulate.toml` and `tiling.tsv`;
   the FASTQ is deleted),
2. filters the step-5 BAMs with `samtools view | awk` on the read name into `subsets/stepK/aln`,
3. runs `unimap score` (it checks pair counts against `tiling.tsv`; all matched).

`scripts/analyze_steps.jl` compares T (MAPQ >= 2) with step 5; worst case over runs and L:

| step | pairs vs step 5 | mean change | p99 change | max change | p99 change, 50 bp bins | 0.9 flips | locus-mean change (max) | locus share T >= 0.9 change (max) |
|---|---|---|---|---|---|---|---|---|
| 10 | 50% | 0.003 | 0.026 | 0.079 | 0.020 | 0.58% | 0.004 | 0.029 |
| 15 | 33% | 0.006 | 0.046 | 0.136 | 0.038 | 0.91% | 0.006 | 0.027 |
| 20 | 25% | 0.008 | 0.062 | 0.205 | 0.053 | 1.28% | 0.016 | 0.039 |
| 25 | 20% | 0.012 | 0.100 | 0.348 | 0.076 | 2.27% | 0.027 | 0.050 |
| 50 | 10% | 0.020 | 0.158 | 0.472 | 0.135 | 3.54% | 0.042 | 0.085 |

- **Step 10 halves the cost** with errors similar to the 8-point F grid, which saved far less.
  Locus-level averages barely move.
- Step 15 is the practical upper limit; from step 20 the per-base structure blurs.
- Part of the step-10 vs step-5 difference is step-5's own sampling noise (step 5 is not
  truth either).
- Consistent with the README rule of thumb `step <= L/5` for the shortest read length.

## 6. Other redundancy considered

- **Adaptive refinement** (coarse pass, then dense tiling only where something fails):
  `scripts/adaptive.jl`. A base counts as "perfect" if every coarse fragment covering it (all F,
  all L) maps back at MAPQ >= 2.

  | coarse step | perfect bases | step-5 fragment coverage wrong there | worst step-5 map-back there |
  |---|---|---|---|
  | 10 | 27.2% | 0.007% | 0.80 |
  | 20 | 29.1% | 0.029% | 0.40 |
  | 25 | 30.5% | 0.058% | 0.30 |
  | 50 | 34.3% | 0.174% | 0.10 |

  At most ~30% could be skipped (L=50 with short F fails almost everywhere in paralogs), and the
  skipped regions hide real failures. Not worth the complexity for a paralog-only panel.
- **R1 is shared across F.** With one orientation, a fragment start yields the same R1 for every
  F, so R1 sequences are 13x redundant. bowtie2 paired mode cannot reuse them. Exploiting it
  would need single-end alignment plus our own pairing, which approximates bowtie2's mate
  rescue and concordance choice instead of measuring it. Rejected.
- **Read lengths.** An L=50 read is a prefix of the L=100 read, but they align differently, so
  nothing can be reused.

## 7. Cost and full-panel estimate

Measured on `subsets/sim` (2.51M pairs per read length, 8 threads, whole-genome index):

| L | bowtie2 time | pairs/s |
|---|---|---|
| 50 | 4.6 min | ~9.0k |
| 100 | 6.1 min | ~6.9k |
| 150 | 7.0 min | ~6.0k |

Scoring all three read lengths with 10 cutoffs took 25 s on 3 threads. The earlier `sim_quick`
estimate was ~10.6k pairs/s overall; the sub-panel rate (~7.1k pairs/s) is the better guide.

**Panel size** (genes from `ds_zebrafish.csv`, 1.5 kb flanks):

| panel | genes | locus length | after merging nearby loci |
|---|---|---|---|
| all paralog pairs | 13,222 | 538 Mb | 477 Mb (~35% of the genome) |
| pairs with dS < 1 | 6,852 | 189 Mb | 155 Mb |
| pairs with dS < 0.5 | 6,129 | 169 Mb | 139 Mb |
| pairs with dS < 0.1 | 3,464 | 109 Mb | 96 Mb |

**Pairs and alignment time** (3 read lengths, plus orientation, 8 threads at 7-10.6k pairs/s):

| grid | pairs per base (3 L) | all paralogs | dS < 0.5 |
|---|---|---|---|
| 8 F, step 5 (`sim_quick`) | 4.8 | 2.3 billion, ~60-90 h | 0.68 billion, ~18-27 h |
| 13 F, step 5 | 7.8 | 3.7 billion, ~100-145 h | 1.1 billion, ~30-43 h |
| **13 F, step 10 (recommended)** | 3.9 | **1.9 billion, ~50-73 h** | 0.54 billion, ~14-21 h |

Storage, scaled from `sim_quick` per-base measurements (scale roughly with pairs for FASTQ/BAM):
Zarr ~5.4 B/base (~2.9 GB for all paralogs, 6 cutoffs), FASTQ ~24 B/pair, BAM ~40 B/pair.
FASTQ can be deleted after alignment; BAMs are needed only to re-score. Scoring memory is about
2 bytes per simulated fragment per read length. Alignment, not storage, is the bottleneck; on
a cluster, split `align.sh` by read length and by batches of regions (not yet implemented).

## 8. Caveats

- Reads are error-free; real misplacement at MAPQ >= 2 will be higher. Untested.
- Sub-panels were stratified by dS, not uniformly random; genes > 60 kb were excluded.
- Only bowtie2 with default end-to-end settings was tested. MAPQ plateaus are aligner-specific.
- The ex-map rate and log2 ratio depend on which loci are in the panel (see README); the
  findings above use mainly Pr(map-back), which does not.
- Real F distributions come from one lab's ChIP runs; ATAC and bisulfite were not weighted.

## 9. Open questions / next steps

1. SNP and sequencing-error sensitivity (e.g. 0.1-1% substitutions, realistic quality
   strings): does MAPQ >= 2 still keep misplacement low?
2. Batching in `align-script` (regions x read lengths) for the full panel on a cluster.
3. Re-check the step choice with error-containing reads and with a second phase offset
   (only phase 0 was tested for each coarse step).
4. Per-run weighting by real fragment-length histograms and read-length quantiles inside
   UniMap (currently done in the ad-hoc analysis scripts).

## 10. Reproducing

Data and scripts live under `test/test_data/` (git-ignored; copy elsewhere to keep them).
Run from the repository root unless noted.

```sh
# sub-panels: genes.txt, pairs.tsv (needs ds_zebrafish.csv, protein_to_gene.csv, GFF)
julia --project test/test_data/subsets/scripts/subsets.jl
cd test/test_data
julia -t 8 --project=../.. -m UniMap loci -f GCF_049306965.2_GRCz12tu_genomic.fna \
  -a GCF_049306965.2_GRCz12tu_genomic.gff -g subsets/genes.txt -o subsets/loci.tsv
julia -t 8 --project=../.. -m UniMap simulate -f GCF_049306965.2_GRCz12tu_genomic.fna \
  -l subsets/loci.tsv -o subsets/sim -L 50,100,150 \
  -F 50,75,100,125,150,175,200,250,300,350,400,450,500 --step 5 --orientation plus
julia --project=../.. -m UniMap align-script -s subsets/sim -x GRCz12tu -t 8
bash subsets/sim/align.sh
julia -t 3 --project=../.. -m UniMap score -s subsets/sim --cutoffs 0,1,2,3,5,8,10,20,30,40
cd ../..
julia --project test/test_data/subsets/scripts/analyze_subsets.jl   # sections 3-4
bash test/test_data/subsets/scripts/step_subsample.sh               # section 5 data
julia --project test/test_data/subsets/scripts/analyze_steps.jl     # section 5
julia --project test/test_data/subsets/scripts/adaptive.jl          # section 6
```

The scripts were written as one-offs: paths are hard-coded relative to the repository root, and
`step_subsample.sh` uses an absolute path to the checkout. `subsets.jl` reads gene coordinates
from `scripts/gff_genes.tsv` (gene, chrom, start, end, strand, biotype for every `gene` feature
in the GFF, extracted with awk). Check paths before re-running elsewhere.
