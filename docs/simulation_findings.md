# Simulation parameter findings

Empirical results from tuning the UniMap simulation grid on zebrafish (GRCz12tu) paralogs.
Written so a person or an LLM can pick up the work without the original chat history.
`README.md` covers motivation, the scoring model and usage. The formal model behind the counts
(what is fixed, measured and latent; recall and precision) is in
[probability_model.md](probability_model.md). This file records what was measured, how, and what
it implies for parameter choices.

Last updated: 2026-10-18. Per-base results (sections 3-5) use coverage weights; see
"Weighting" in section 1.

## TL;DR: recommended settings

| parameter | recommendation | why (section) |
|---|---|---|
| fragment MAPQ cutoff | **>= 2** (fragment MAPQ = min of the two mates) | removes bowtie2 ties; >= 10 costs ~11 points of map-back for little gain (3) |
| orientation | `--orientation plus` | same result as `both` at MAPQ >= 2, half the reads (2) |
| read lengths | 50, 100, 150, or whatever the real runs round down to, taken from read lengths before adapter trimming | each L is a separate bowtie2 run; L=50 is clearly worse (3); adapter trimming is already simulated (probability_model section 5) |
| fragment lengths | 50,75,100,125,150,175,200,250,300,350,400,450,500 | 25 bp spacing over 100-300 is needed; 75 and 400 are the cheapest to drop (4) |
| `--step` | **10** (15 at most) | halves the reads vs step 5 with small error; >= 20 blurs per-base structure (5) |

Expected full-panel cost with these settings is in section 8.

## 1. Context

- **Goal.** Per-base empirical paired-end mappability for paralog loci, to be combined later with
  each real run's fragment-length histogram and read length into confidence formulas for
  enrichment of histone marks (especially broad marks). Stored quantities are simple counts
  (Pr(map-back), ex-map rate, log2 ratio; see README); `probability_model.md` defines recall and
  precision in terms of them.
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

  Each run's per-bp density is a piecewise-linear CDF through the quantiles, from q10/2 up to
  the largest F (500). It is projected onto the F grid with linear-interpolation (hat) weights
  that sum to 1. The CDF is forced to reach 1 at 500, so the tail beyond `-X` is folded into
  361-500; all recall figures are conditional on F <= 500 (probability_model section 5).
- **Weighting.** Two weightings, for two estimands (probability_model section 5):
  - *Fragment weights* $p_F$ (share of fragments), for fragment-level summaries from
    `misplacement_summary` (misplaced %, fragment-level map-back).
  - *Coverage weights* $v_F \propto \int p(f)\,f\,h_F(f)\,df$ (share of coverage; fragments of
    length F cover F bases), for per-base tracks.

  | F | 50 | 75 | 100 | 125 | 150 | 175 | 200 | 250 | 300 | 350 | 400 | 450 | 500 |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|---|
  | fragment, typical | .004 | .034 | .040 | .079 | .136 | .150 | .177 | .173 | .099 | .037 | .028 | .028 | .014 |
  | coverage, typical | .001 | .012 | .018 | .048 | .095 | .121 | .170 | .196 | .136 | .059 | .052 | .059 | .032 |

  F < 100 carries 0-5% of the weight for ChIP; ATAC would weight short fragments more. An
  earlier version of this document applied fragment weights to per-base tracks, which
  under-weights long fragments; sections 3-5 now use coverage weights.
- **Run-weighted track.** For a base x, read length L and cutoff c:
  `T(x) = sum_F v_F * n_correct[x,F,c,L] / n_simulated[x,F,L]`, which approximately equals
  `sum_F p_F n_correct / sum_F p_F n_simulated`. Most comparisons below are on this track, and
  on whether it is >= 0.9 ("0.9 flips" = share of bases whose call changes).
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
- Analysis: `scripts/analyze_subsets.jl`.

**bowtie2 MAPQ plateaus.** In these data, results are identical within the groups {0, 1},
{2, 3, 5}, {8, 10, 20, 30} and {40}, so only one value per group matters.

**All three sub-panels pooled.** Misplaced % (wrong / kept), map-back (correct / simulated) and
ectopic / kept are fragment-level, with typical-run fragment weights. The last column is
per-base: the share of bases with T >= 0.9, with typical-run coverage weights.

| L | MAPQ >= | misplaced % | Pr(map-back) | ectopic / kept | bases T >= 0.9 |
|---|---|---|---|---|---|
| 50 | 0 | 24.6 | 0.754 | 0.043 | 52.8% |
| 50 | 2 | 0.098 | 0.686 | 2e-5 | 48.1% |
| 50 | 10 | 0.041 | 0.605 | 1e-5 | 40.6% |
| 50 | 40 | 0.020 | 0.357 | 1e-5 | 22.6% |
| 100 | 0 | 17.4 | 0.826 | 0.034 | 69.8% |
| 100 | 2 | 0.021 | 0.773 | 1e-5 | 66.3% |
| 100 | 10 | 0.005 | 0.660 | <1e-5 | 51.4% |
| 100 | 40 | 0.004 | 0.363 | <1e-5 | 25.4% |
| 150 | 0 | 16.2 | 0.838 | 0.032 | 72.9% |
| 150 | 2 | 0.006 | 0.787 | <1e-5 | 69.9% |
| 150 | 10 | 0.001 | 0.668 | <1e-5 | 54.6% |
| 150 | 40 | 0.001 | 0.354 | <1e-5 | 25.1% |

- **MAPQ >= 2 removes the ties.** Misplacement falls from 16-25% to <= 0.1%, and ectopic
  fragments per kept fragment from 3-4% to ~1e-5.
- **The sub-panels agree.** At MAPQ >= 2 misplacement ranged 0.003-0.14% across panels and
  read lengths (worst: panel 3, L=50).
- **MAPQ >= 10** removes a further 0.015-0.06 points of misplacement but loses ~11 points of
  map-back and 7-15 points of bases with T >= 0.9.
- **MAPQ 40 is unusable**: map-back roughly halves.
- Map-back rises with read length; the step from 50 to 100 is much larger than 100 to 150.

## 4. Fragment-length grid

Greedy backward elimination from the 13-point grid at MAPQ >= 2, always keeping 50 and 500,
removing at each step the F whose loss changes T least (by p99). Errors are |T_reduced -
T_full| with coverage weights, worst case over the short/typical/long runs and the three read
lengths.

| kept fragment lengths | max change | p99 change | 0.9 flips |
|---|---|---|---|
| all 13 | 0 | 0 | 0 |
| drop 75 (12) | 0.012 | 0.006 | 0.15% |
| drop 75, 400 (11) | 0.023 | 0.010 | 0.25% |
| drop 75, 400, 175 (10) | 0.044 | 0.023 | 0.65% |
| 50,100,150,200,250,300,350,450,500 (9) | 0.053 | 0.025 | 0.70% |
| 50,100,150,200,250,300,350,500 (8) | 0.059 | 0.029 | 0.75% |
| 50,100,150,200,250,350,500 (7) | 0.087 | 0.040 | 1.14% |
| 50,100,150,250,350,500 (6) | 0.196 | 0.070 | 1.96% |
| 50,100,150,250,500 (5) | 0.219 | 0.099 | 2.71% |
| 50,150,500 (3) | 0.484 | 0.301 | 5.36% |

- Each F costs about the same to simulate and align, so dropping one F saves ~1/13 (~8%).
- 75 and 400 are the cheapest drops (11 lengths: p99 0.010, 0.25% flips); a third drop doubles
  the error. With fragment weights (earlier version) the order was 400 then 75, so which length
  is cheapest depends on the weighting, but both are cheap under either.
- The 8-point grid used in `sim_quick` is too coarse (0.75% of bases flip): most ChIP coverage
  comes from 150-300 bp fragments, which needs ~25-50 bp spacing.
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

`scripts/analyze_steps.jl` compares T (MAPQ >= 2, coverage weights) with step 5; worst case over
runs and L:

| step | pairs vs step 5 | mean change | p99 change | max change | p99 change, 50 bp bins | 0.9 flips | locus-mean change (max) | locus share T >= 0.9 change (max) |
|---|---|---|---|---|---|---|---|---|
| 10 | 50% | 0.003 | 0.024 | 0.066 | 0.020 | 0.51% | 0.004 | 0.045 |
| 15 | 33% | 0.005 | 0.040 | 0.124 | 0.036 | 0.78% | 0.006 | 0.029 |
| 20 | 25% | 0.007 | 0.055 | 0.193 | 0.048 | 1.09% | 0.016 | 0.058 |
| 25 | 20% | 0.010 | 0.085 | 0.303 | 0.068 | 1.94% | 0.028 | 0.050 |
| 50 | 10% | 0.017 | 0.141 | 0.414 | 0.122 | 3.05% | 0.044 | 0.104 |

- **Step 10 halves the cost** with errors similar to dropping three fragment lengths, which
  saves far less. Locus-mean map-back barely moves; the locus share of bases >= 0.9 is more
  sensitive because it is a threshold count.
- Step 15 is the practical upper limit; from step 20 the per-base structure blurs.
- Part of the step-10 vs step-5 difference is step-5's own sampling noise (step 5 is not
  truth either).
- Only phase 0 was tested for each coarse step. A fixed phase is a systematic subsample:
  unbiased on average over phases, but each base's error depends on the phase.
- When step does not divide F (step 10 with F = 75, 125, 175), N(x) alternates between
  floor(F/step) and ceil(F/step) along the locus; the ratios account for this.
- Consistent with the README rule of thumb `step <= L/5` for the shortest read length.

## 6. Other redundancy considered

- **Adaptive refinement** (coarse pass, then dense tiling only where something fails):
  `scripts/adaptive.jl`. A base counts as "perfect" if every coarse fragment covering it (all F,
  all L) maps back at MAPQ >= 2. No weighting is involved.

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

## 7. Where misplaced fragments go

`scripts/destinations.jl`, on `subsets/sim` (step 5). Destinations come from `bleed_matrix`
(misplaced kept fragments, all F pooled): "same locus" = shifted within the origin locus, "other
locus" = another sub-panel locus, "off panel" = anywhere outside the 72 loci.

| L | MAPQ >= | misplaced kept | same locus | other panel locus | off panel |
|---|---|---|---|---|---|
| 50 | 0 | 639,284 | 11.0% | 2.2% | 86.8% |
| 50 | 2 | 1,876 | 3.3% | 0.3% | 96.4% |
| 50 | 10 | 831 | 5.2% | 0.2% | 94.6% |
| 100 | 0 | 489,188 | 12.5% | 2.2% | 85.3% |
| 100 | 2 | 519 | 3.3% | 0.2% | 96.5% |
| 100 | 10 | 186 | 1.6% | 0.5% | 97.8% |
| 150 | 0 | 451,408 | 12.7% | 2.1% | 85.2% |
| 150 | 2 | 238 | 5.5% | 0.0% | 94.5% |
| 150 | 10 | 72 | 0% | 0% | 100% |

Per-base ectopic counts at MAPQ >= 2 (all loci, all F):

| L | (base, F) cells with E > 0 | max ex-map rate e | 99.99th pct of e |
|---|---|---|---|
| 50 | 0.094% | 0.120 | 0.070 |
| 100 | 0.033% | 0.067 | 0.025 |
| 150 | 0.040% | 0.020 | 0.017 |

- **Most misplacement leaves the panel.** Since the destination sites resemble the panel loci,
  fragments from those sites probably land on the panel at a comparable rate, and that flow is
  not simulated. `n_ectopic`, the ex-map rate and the log2 ratio therefore miss most of the
  real ectopic input (probability_model sections 3 and 6). Some of these sites are paralogs that
  a full panel would include; repeats such as TEs would remain outside.
- **At MAPQ >= 2, ectopic coverage within the panel is negligible**: equal-signal precision is
  ~1 almost everywhere, and the per-base problem is loss (low map-back), not confusion.
- **Shifts within the origin locus** are 11-13% of misplacement at MAPQ 0 and 3-6% at MAPQ 2.
  At tolerance 0 they count against both map-back and the ex-map rate even though the fragment
  still largely covers the right bases.

## 8. Cost and full-panel estimate

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

Dropping 75 and 400 (section 4) would cut these by a further ~15%.

Storage, scaled from `sim_quick` per-base measurements (scale roughly with pairs for FASTQ/BAM):
Zarr ~5.4 B/base (~2.9 GB for all paralogs, 6 cutoffs), FASTQ ~24 B/pair, BAM ~40 B/pair.
FASTQ can be deleted after alignment; BAMs are needed only to re-score. Scoring memory is about
2 bytes per simulated fragment per read length. Alignment, not storage, is the bottleneck; on
a cluster, split `align.sh` by read length and by batches of regions (not yet implemented).

## 9. Caveats

- Reads are error-free; real misplacement at MAPQ >= 2 will be higher. Untested.
- Sub-panels were stratified by dS, not uniformly random; genes > 60 kb were excluded.
- Only bowtie2 with default end-to-end settings was tested. MAPQ plateaus are aligner-specific.
- The ex-map rate and log2 ratio depend on which loci are in the panel (see README and
  section 7); the findings above use mainly Pr(map-back), which does not.
- Real F distributions come from one lab's ChIP runs; ATAC and bisulfite were not weighted.
  Fragments longer than 500 bp were folded into the grid (section 1), so recall figures are
  conditional on F <= 500.
- Coarse steps were tested at one phase only (section 5).

## 10. Open questions / next steps

1. SNP and sequencing-error sensitivity (e.g. 0.1-1% substitutions, realistic quality
   strings): does MAPQ >= 2 still keep misplacement low?
2. Off-panel sources: simulate the off-panel destination sites from `misplaced_sites` as extra
   origins, to measure the ectopic input the panel currently misses (section 7).
3. Batching in `align-script` (regions x read lengths) for the full panel on a cluster.
4. Re-check the step choice with error-containing reads and with other phase offsets.
5. Per-run weighting by real fragment-length histograms and read-length quantiles inside
   UniMap (currently done in the ad-hoc analysis scripts), with coverage weights for per-base
   tracks and fragment weights for per-locus summaries.

## 11. Reproducing

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
julia --project test/test_data/subsets/scripts/destinations.jl      # section 7
```

The scripts were written as one-offs: paths are hard-coded relative to the repository root, and
`step_subsample.sh` uses an absolute path to the checkout. `subsets.jl` reads gene coordinates
from `scripts/gff_genes.tsv` (gene, chrom, start, end, strand, biotype for every `gene` feature
in the GFF, extracted with awk). Check paths before re-running elsewhere.
