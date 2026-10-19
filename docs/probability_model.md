# From simulation counts to probabilities

How UniMap's simulated counts relate to the recall and precision of real ChIP coverage at paralog
loci: what we fix, what the simulation measures, and what stays latent. Measured parameter
choices are in [simulation_findings.md](simulation_findings.md); usage is in the README.

## 1. Notation

Throughout, one read length $L$, one fragment length $F$ and one fragment MAPQ cutoff $c$ are
fixed (a **channel**); every symbol carries an implicit $(L, F, c)$ index unless stated
otherwise. $x$ and $y$ are bases inside the panel.

**Simulated fragments.** A simulated fragment is $\varphi = (s, F, o)$: start $s$ (0-based),
length $F$, orientation $o$. For each simulation region and $F$, starts lie on an arithmetic grid
with spacing `step` and a fixed phase (`fragment_starts` in `src/tiler.jl`), so there is one
simulated fragment per `step` bp of start positions per orientation. Its **true span** is
$[s, s+F)$; its **placed span** runs from the leftmost aligned position over $|\mathrm{TLEN}|$
bases. A fragment **covers** a base if the base lies in the span (as in MACS3 `BAMPE`).

Every simulated fragment ends in exactly one outcome:

| outcome | meaning |
|---|---|
| correct, kept | primary proper pair on the origin chromosome; leftmost position and $\lvert\mathrm{TLEN}\rvert$ within `--tolerance` of the origin; R1 on the expected strand; fragment MAPQ (min of the two mates) $\ge c$ |
| misplaced, kept | any other proper pair with fragment MAPQ $\ge c$, including placements shifted within its own locus |
| below cutoff | proper pair, correct or not, with fragment MAPQ $< c$ |
| not proper | unaligned, discordant, or only one mate aligned |
| missing | absent from the BAM (normally none) |

| symbol | definition | stored as |
|---|---|---|
| $N(x)$ | simulated fragments whose true span covers $x$ | `n_simulated` |
| $C(x)$ | of those, correct and kept | `n_correct` |
| $W(x), B(x), U(x), M(x)$ | of those, misplaced and kept / below cutoff / not proper / missing | per locus only, from `misplacement_summary` |
| $E(y)$ | simulated fragments from any origin that are misplaced and kept and whose placed span covers $y$ | `n_ectopic` |
| $R(x) = C(x)/N(x)$ | Pr(map-back) | `map_back` |
| $e(y) = E(y)/N(y)$ | ex-map rate | `ectopic_rate` |
| $\ell(y) = \log_2\frac{C(y) + 1/2}{E(y) + 1/2}$ | log2 ratio | `log2_ratio` |

$N = C + W + B + U + M$. $C$, $W$ and $B$ depend on $c$; $U$ and $M$ do not. $W(x)$ and $E(y)$
count the same misplaced fragments from opposite ends: at the origin and at the destination.

**Origin units.** To split $E$ by origin, partition the simulated fragments (not bases) into
origin units $u$: single fragments, or windows of fragment starts or midpoints. Then
$E(y) = \sum_u E_u(y)$, where $E_u(y)$ counts the fragments in $u$ that are misplaced, kept and
whose placed span covers $y$; each fragment belongs to exactly one $u$. `bleed_matrix` stores
fragment counts (not coverage) with $u$ = origin locus and destination = the locus with the
largest overlap of the placed span. Per-base $E_u(y)$ is not stored.

## 2. What is fixed, measured and latent

**Fixed a priori (explicit parameters and design choices)**

- Reference genome, used as a stand-in for the sample's genome.
- The panel: which loci are simulated. Only panel origins contribute to $E$.
- Read length(s) $L$ (both mates equal), fragment-length grid $\{F\}$, `--step` and its phase,
  orientation.
- Error-free reads with uniform base quality.
- Aligner and settings (bowtie2 end-to-end), fragment MAPQ rule, cutoff $c$, tolerance.
- **Uniform tiling**: the same number of simulated fragments per bp of starts at every locus. This
  is a design choice that gives every origin the same exposure and makes $\lambda$ (below) well
  defined. It is not an assumption about the sample; equal *real* signal is a separate,
  optional assumption (section 3).

$N(x)$ is therefore set by us. Its mean is $F \cdot o / \text{step}$ ($o$ = number of simulated
orientations). It is exact at every base only when `step` divides $F$; otherwise it alternates
between $\lfloor F/\text{step} \rfloor\,o$ and $\lceil F/\text{step} \rceil\,o$ along the locus
(e.g. 17 and 18 for $F = 175$, step 10). It is lower near non-ACGT bases and chromosome ends.

**Measured by the simulation**

The outcome of every simulated fragment, hence $C$, $E$, $E_u$ (at locus resolution), $W$, $B$,
$U$: a placement kernel from origins to destinations under the fixed parameters.

- At $c \ge 2$ this is a census of the aligner's behaviour on the simulated fragments: bowtie2 is
  deterministic per read, and its random tie-breaking only affects pairs with MAPQ 0-1.
- At $c \in \{0, 1\}$ each tied fragment's placement is one pseudo-random draw, so $C$ and $E$
  carry binomial noise.
- At step > 1 the simulated fragments are a systematic subsample of all starts at one phase. The
  counts are exact for those fragments but estimate the all-starts values with a phase-dependent
  error; they are unbiased only on average over phases.
- Between grid values of $F$, values are interpolated.

**Latent (never observed in real data)**

- $\mu_{F}(s)$: expected number of real fragments of length $F$ starting at $s$, per
  orientation. Define $\lambda_F(s) = \text{step} \cdot \mu_F(s)$, the number of real fragments
  per simulated fragment at that start. It combines enrichment, sequencing depth and local
  biases.
- The true origin of each real read.
- Differences between the sample and the reference: SNPs, indels, copy number, copies absent
  from the reference, collapsed duplications.
- Sequencing errors.
- The run's true fragment-length distribution $p_F$. It is observed only through proper pairs,
  which are truncated at `-X` and filtered by mappability.

**Modelling assumptions** (not measured; each is named where it is used):

- **A1, separability.** $\mu_F(s) = a(s)\,p_F$: one fragment-length distribution at every origin.
  Nucleosome phasing near peaks or PCR length bias can break it locally.
- **A2, local constancy.** $a(s)$ is constant over the starts whose fragments cover a base (an
  $F$-bp window of starts), and over each origin unit $u$ used to split $E$. With the locus-level
  `bleed_matrix`, the latter means constant across a whole source locus.
- **A3, transfer.** A real fragment with $(s, F, o)$ is placed exactly as its simulated
  counterpart. This needs: no sequencing errors; sample = reference; real mates of length $L$;
  placement independent of base qualities; orientation symmetry (the simulation uses `plus`
  only; findings section 2 supports this at $c \ge 2$); end-to-end alignment (soft clipping in
  `--local` mode would make correct placements fail at tolerance 0).

## 3. Recall and precision for one channel

Under A2 and A3, write $\lambda(y)$ for $\lambda_F$ at the starts covering $y$. The expected real
coverage at $y$ in this channel, in expectation over the tiling phase, is:

- truly from $y$ (any outcome): $\lambda(y)\,N(y)$
- truly from $y$, correct and kept: $\lambda(y)\,C(y)$
- placed at $y$ and kept but not native: $\sum_u \lambda_u\,E_u(y) + X(y)$

where $X(y)$ is kept coverage from origins that were not simulated (off the panel, or absent
from the reference). So the expected observed (kept) coverage is

$$\mathbb{E}[O(y)] = \lambda(y)\,C(y) + \sum_u \lambda_u\,E_u(y) + X(y).$$

**Recall** (of the fragments truly covering $y$, the share recovered correctly):

$$\mathrm{Rec}(y) = \frac{\lambda(y)\,C(y)}{\lambda(y)\,N(y)} = R(y).$$

$\lambda(y)$ cancels, so recall is a property of the locus. Without A2,
$\mathrm{Rec}(y) = \sum_\varphi \lambda(\varphi)\,\mathbb{1}[\varphi\ \text{correct, kept}] / \sum_\varphi \lambda(\varphi)$
over the fragments covering $y$, and $R(y)$ is its $\lambda$-uniform case.

**Precision** (of the kept coverage at $y$, the share that is native):

$$\mathrm{Prec}(y) = \frac{\lambda(y)\,C(y)}{\lambda(y)\,C(y) + \sum_u \lambda_u\,E_u(y) + X(y)}.$$

"Native" means: from a fragment that truly covers $y$ and was placed correctly. A fragment
shifted within its own locus whose placed span still covers $y$ is counted in $E$, not $C$, so
at tolerance 0 such shifts count against both recall and precision. For broad marks, a looser
definition of "correct" (a tolerance, or span overlap) may suit better. With tolerance > 0 the
scorer also adds correct fragments over the *simulated* span rather than the placed one, which
differs by at most the tolerance.

**Equal signal, complete panel** ($\lambda_u = \lambda(y)$ for all $u$, and $X = 0$):

$$\mathrm{Prec}_{\mathrm{eq}}(y) = \frac{C(y)}{C(y) + E(y)} = \frac{R(y)}{R(y) + e(y)}.$$

Because $\lambda$ is defined per simulated fragment, this does not require origins to have the
same $N$: $N(y)$ cancels between $R$ and $e$. The pseudocount gives
$1/(1 + 2^{-\ell}) = (C + 1/2)/(C + E + 1)$, which equals $\mathrm{Prec}_{\mathrm{eq}}$ without
it. Where $E = 0$, $\ell = \log_2(2C + 1)$ is a floor that reflects how much was simulated
(more for longer $F$ or smaller step), not the true odds.

**Unequal signal.** Let $\bar\rho(y) = \sum_u \lambda_u E_u(y) / (\lambda(y)\,E(y))$, the
$E$-weighted mean ratio of source signal to the signal at $y$. With $X = 0$:

$$\mathrm{Prec}(y) = \frac{R(y)}{R(y) + \bar\rho(y)\,e(y)}, \qquad
\log_2\text{-odds} = \log_2\frac{C(y)}{E(y)} - \log_2\bar\rho(y).$$

A single source, or several sharing one ratio $\rho$, gives $\bar\rho = \rho$.

## 4. Which way departures from the simulation push

| departure | effect on real recall vs $R$ | effect on real precision vs $\mathrm{Prec}_{\mathrm{eq}}$ |
|---|---|---|
| sequencing errors (random per read) | lower, in expectation | lower, in expectation |
| sample variants vs reference (fixed per locus) | lower at the affected locus | lower at destinations the variant makes more similar |
| real reads longer than the simulated $L$ (e.g. $L$ rounded down) | higher | usually higher |
| $L$ derived from adapter-trimmed lengths (trimming counted twice) | higher | usually higher |
| fragments longer than `-X` (never proper pairs) | lower (recall 0 for them; similar at every locus) | none |
| origins outside the panel, $X > 0$ | none | lower |
| copies in the sample but collapsed in the reference | none | lower; invisible to any simulation |
| copy in the reference but absent from the sample ($\lambda_u = 0$) | none | higher |
| sources more strongly marked than $y$ ($\bar\rho > 1$) | none | lower |
| sources less strongly marked than $y$ ($\bar\rho < 1$) | none | higher |
| $\lambda$ not constant within a fragment span (sharp peaks) | either | either |
| aligner seed heuristics, MAPQ thresholds | either, small | either, small |

Consequences:

- **$R$ is not a bound on real recall under recommended practice.** Errors, variants and long
  fragments push real recall below $R$; rounding $L$ down pushes it above. $R$ is an estimate
  whose biases have opposite signs.
- **$\mathrm{Prec}_{\mathrm{eq}}$ is a reference point, not a bound.** It is an upper bound only
  with respect to omitted sources ($X \ge 0$) and, in expectation, sequencing errors.
- **Errors and variants differ.** Errors are random per read and average out over reads.
  Variants act on every read from a locus, so their effect at one locus does not average out and
  can be large. Neither can break a tie between reference copies that are identical over the
  read: alignment is against the reference, so a mismatch in the read hits both copies. A
  mismatch can only move a read towards another copy's allele or towards a tie; precision can
  rise only through aligner heuristics, MAPQ thresholds or a copy missing from the sample.

## 5. Combining fragment and read lengths for a run

Two different estimands:

- **Base level** (per-base tracks): the probability that a real fragment covering $y$ maps back.
  Fragments of length $F$ cover $F$ bases, so they weigh in proportion to $p_F F$.
- **Fragment level** (per-locus summaries): the share of a locus's real fragments that map back.
  Fragments weigh in proportion to $p_F$.

For the base level, with $p_F$ normalised to $\sum_F p_F = 1$, sum counts before dividing:

$$\tilde C(y) = \sum_F p_F C_F(y), \quad \tilde N(y) = \sum_F p_F N_F(y), \quad \tilde E(y) = \sum_F p_F E_F(y),$$

$$R_p(y) = \frac{\tilde C(y)}{\tilde N(y)}, \qquad e_p(y) = \frac{\tilde E(y)}{\tilde N(y)}, \qquad
\ell_p(y) = \log_2\frac{\tilde C(y) + 1/2}{\tilde E(y) + 1/2}.$$

Since $N_F \propto F$, $R_p = \sum_F v_F R_F$ with **coverage weights** $v_F \propto p_F F$.
Normalising $p_F$ keeps $\tilde C$ and $\tilde E$ on the count scale of one channel, so the
pseudocount keeps its meaning. Never average $\ell$ across channels.

**Projection onto the grid.** If $R(f)$ is interpolated linearly between grid points with hat
functions $h_k$, the consistent coverage weights are

$$v_k = \frac{\int p(f)\, f\, h_k(f)\, df}{\int p(f)\, f\, df},$$

which differ from $p_k F_k$ by a few per cent on an unevenly spaced grid (the hat at 200 with
neighbours 175 and 250 has its centroid at 208).

**Fragments longer than the largest $F$** (= `-X`) never align as proper pairs. Either keep their
coverage share in the denominator, with zero kept, or state that recall is conditional on
$F \le X$. The loss is nearly the same at every locus, so it affects absolute recall but not
comparisons between loci or precision. Fragment-length quantiles from proper pairs are themselves
truncated at `-X`, so the real tail is larger than they show.

**Read lengths.** $N_{F,L}$ does not depend on $L$ (except near non-ACGT bases), so summing
counts with read-length shares $q_L$ is the same as weighting $R_L$ by $q_L$. If $L$ and $F$ are
dependent, use their joint distribution. Derive $L$ from read lengths **before adapter
trimming** (after quality trimming only): the simulation itself shortens reads of fragments
shorter than $L$ to $F$, so deriving $L$ from adapter-trimmed lengths counts that trimming twice
and biases $L$ downward. Unequal R1 and R2 lengths are not modelled.

## 6. What the recommended cutoff does (measured)

From the random sub-panels ([simulation_findings.md](simulation_findings.md), section 7), with
error-free reads:

- **Precision is close to 1 by construction at $c \ge 2$.** Only 0.03-0.09% of (base, $F$) cells
  have $E > 0$; the largest per-base $e$ is 0.12 at $L = 50$ and 0.02 at $L = 150$. Precision
  then moves only for $\bar\rho$ in the thousands, except in those few cells. At the recommended
  cutoff, mappability is mainly a **loss** problem, not a confusion problem: under A1-A3, the
  signal at $y$ can be estimated as $\hat\lambda(y) \propto O(y)/R(y)$ wherever $R(y)$ is not
  small, and masked elsewhere.
- **Most misplaced fragments land off the panel.** Of the misplaced kept fragments from the
  sub-panel loci, 85-87% (at $c = 0$) and 94-97% (at $c = 2$) land outside the panel. Those sites
  resemble the panel loci, so the reverse flow $X(y)$ is probably of comparable size, and it is
  not measured at all. With the current outputs, precision is dominated by what is not
  simulated: $X$, errors and variants.

## 7. Towards estimating $\lambda$ (kernel inversion)

A research direction, not a validated estimator.

- **Units: fragment counts.** Assign each simulated fragment to an origin window by its true
  midpoint and, if kept, to a destination window by its placed midpoint. Define
  $K_{ww'}$ = (simulated fragments from $w'$ kept with placed midpoint in $w$) / (simulated
  fragments from $w'$), counting correct and misplaced fragments alike. Correct fragments and
  same-window misplacements fall on the diagonal; the midpoint convention avoids spill-over of
  correct fragments into neighbouring windows.
- **Observation model.** $O_w$ = real kept (deduplicated) fragments with placed midpoint in $w$,
  and $\mathbb{E}[O_w] = \sum_{w'} K_{ww'} \Lambda_{w'} + X_w$, with $\Lambda_{w'}$ the expected
  number of real fragments originating in $w'$. Use a Poisson likelihood on these counts, or a
  negative binomial for biological variation and PCR duplicates. Do not use summed per-base
  coverage: each fragment adds to about $F$ bases, so coverage sums are overdispersed and
  autocorrelated.
- **Fragment and read lengths.** Under A1, $K = \sum_{F} p_F K_F$ with fragment-level weights
  $p_F$ (not coverage weights), and likewise over $L$.
- **Identifiability** depends on the column structure of $K$, not on $R$ alone. At $c \ge 2$
  (error-free) $K$ is essentially diagonal and the inversion reduces to
  $\hat\Lambda_w = O_w / K_{ww}$; windows with $K_{ww} \approx 0$ are not identified. At $c = 0$,
  identical copies split their fragments evenly through random tie-breaking, their columns are
  nearly collinear, and only $\Lambda_w + \Lambda_{w'}$ is identified. A window with low $R$ can
  still be informative if its fragments land somewhere distinguishable. Where columns are
  collinear, report the identifiable combinations (e.g. sums over copies), not per-window values.
- **Omitted sources.** $X_w$ is an unmodelled column; leaving it out biases $\hat\Lambda$ upward
  at the destinations it feeds.
- **Needs a new output**: a per-window kernel. `bleed_matrix` (locus pairs, fragment counts)
  supports only locus resolution and needs A2 across each whole source locus.

## 8. Open points

1. Error and variant models, so that $R$ becomes a calibrated estimate rather than an estimate
   with known-sign biases.
2. Off-panel sources: simulate the destination sites in `misplaced_sites` (or a genome-wide
   background sample) as extra origins to estimate $X$.
3. Per-window kernel output for section 7.
4. Phase sensitivity of coarse steps (only one phase per step has been tested).
5. Whether "correct" should allow shifts within a locus for broad marks.
6. Sensitivity of the inversion to A1 and A2.
