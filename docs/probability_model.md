# From simulation counts to probabilities

How UniMap's simulated counts relate to recall and precision of real ChIP coverage at paralog
loci: what we fix, what the simulation measures, and what stays latent. Measured parameter
choices are in [simulation_findings.md](simulation_findings.md); usage is in the README.

## 1. Notation

Throughout, one read length $L$, one fragment length $F$ and one fragment MAPQ cutoff $c$ are
fixed; every symbol below carries an implicit $(L, F, c)$ index unless stated otherwise.
Positions $x, y$ are bases inside the panel; $z$ indexes origins (simulated bases or regions).

A simulated fragment is **kept** if its fragment MAPQ (min of the two mates) is $\ge c$, and
**correct** if it is a primary proper pair at its origin (within `--tolerance`, R1 on the
expected strand). A fragment **covers** a base if the base lies in its span (leftmost to
rightmost aligned position, as in MACS3 `BAMPE`).

| symbol | definition | stored as |
|---|---|---|
| $N(x)$ | simulated fragments whose true span covers $x$ | `n_simulated` |
| $C(x)$ | of those, placed correctly and kept | `n_correct` |
| $E(y)$ | simulated fragments, from any origin, placed incorrectly, kept, and whose *placed* span covers $y$ | `n_ectopic` |
| $E_z(y)$ | the part of $E(y)$ coming from origin $z$, so $E(y) = \sum_z E_z(y)$ | per base: not stored; per locus pair: `bleed_matrix` |
| $R(x) = C(x)/N(x)$ | Pr(map-back) | `map_back` |
| $e(y) = E(y)/N(y)$ | ex-map rate | `ectopic_rate` |
| $\ell(y) = \log_2\frac{C(y) + 1/2}{E(y) + 1/2}$ | log2 ratio | `log2_ratio` |

Every simulated fragment covering $x$ ends in exactly one of four outcomes, so
$N(x) = C(x) + W(x) + B(x) + U(x)$: correct and kept ($C$), misplaced and kept ($W$),
below the cutoff wherever placed ($B$), not a proper pair ($U$). Only $C$ and $E$ are stored
per base; $W$ is `n_wrong_kept` in `misplacement_summary` (per locus). Note that $W$ and $E$
count the same misplaced fragments from opposite ends: $W(x)$ at the origin, $E(y)$ at the
destination.

## 2. What is fixed, measured and latent

**Fixed a priori (explicit parameters and modelling choices)**

- Reference genome, used as a stand-in for the sample's genome.
- The panel: which loci are simulated (only they can contribute to $E$).
- Read length(s) $L$, fragment-length grid $\{F\}$, tiling `--step`, orientation.
- Error-free reads with uniform base quality.
- Aligner and its settings; fragment MAPQ rule; cutoff $c$; tolerance.
- Uniform tiling: every base is the start of the same number of fragments per $F$. This is
  the **equal-signal** assumption, and it is why $N(x)$ is (away from N runs) the same at every
  base: $N(x) = F \times (\text{orientations}) / \text{step}$.

$N(x)$ is therefore set by us, not estimated.

**Measured by the simulation**

$C(x)$, $E(y)$, $E_z(y)$ (at locus resolution), $W$, $B$, $U$. Together they describe how the
aligner moves fragments: a placement kernel from origins to destinations under the fixed
parameters. bowtie2 is deterministic per read (only its tie-breaking is random, and ties have
MAPQ 0-1), so at step 1 the counts are a census rather than a sample. At step > 1 they are a
systematic subsample along the locus, and between grid values of $F$ they are interpolated
(see the findings document for the size of both errors).

**Latent (never observed in real data)**

- $\lambda(z)$: the real fragment density at origin $z$ relative to the uniform simulation,
  i.e. real fragments per simulated fragment. It combines enrichment, sequencing depth and
  local biases, and is what an enrichment track ultimately tries to estimate.
- The true origin of each real read.
- Differences between the sample and the reference: SNPs, indels, copy number, copies absent
  from the reference.
- Sequencing errors.
- The run's true fragment-length distribution (observed only through proper pairs, which are
  themselves filtered by mappability).

## 3. Recall and precision under the model

Assume $\lambda$ is roughly constant across one fragment span around each base (reasonable
for broad marks, weaker for sharp peaks). Real fragments then behave like simulated ones,
scaled by $\lambda$. At base $y$, the expected real coverage is:

- truly from $y$: $\lambda(y)\,N(y)$
- truly from $y$ and placed correctly and kept: $\lambda(y)\,C(y)$
- placed at $y$ and kept but from elsewhere: $\sum_z \lambda(z)\,E_z(y)$

so the expected observed (kept) coverage is

$$\mathbb{E}[O(y)] = \lambda(y)\,C(y) + \sum_z \lambda(z)\,E_z(y).$$

**Recall** (of the fragments truly from $y$, the share recovered at $y$):

$$\mathrm{Rec}(y) = \frac{\lambda(y)\,C(y)}{\lambda(y)\,N(y)} = \frac{C(y)}{N(y)} = R(y).$$

$\lambda(y)$ cancels: recall is a property of the locus, and the simulation measures it
directly.

**Precision** (of the kept coverage at $y$, the share truly from $y$):

$$\mathrm{Prec}(y) = \frac{\lambda(y)\,C(y)}{\lambda(y)\,C(y) + \sum_z \lambda(z)\,E_z(y)}
= \frac{1}{1 + \sum_z \frac{\lambda(z)}{\lambda(y)}\,\frac{E_z(y)}{C(y)}}.$$

It depends on the latent signal ratios $\lambda(z)/\lambda(y)$.

**Equal signal** ($\lambda$ constant, the simulation's own assumption). Using
$\sum_z E_z(y) = E(y)$ and dividing by $N(y)$:

$$\mathrm{Prec}_{\mathrm{eq}}(y) = \frac{C(y)}{C(y) + E(y)} = \frac{R(y)}{R(y) + e(y)},
\qquad \mathrm{Prec}_{\mathrm{eq}}(y) \approx \frac{1}{1 + 2^{-\ell(y)}}.$$

The approximation is exact without the pseudocount; $\ell(y)$ is the log2 odds that a kept
fragment at $y$ is native, under equal signal.

**One shared ratio.** If every contributing origin has $\lambda(z) = \rho\,\lambda(y)$:

$$\mathrm{Prec}(y) = \frac{C(y)}{C(y) + \rho\,E(y)} = \frac{R(y)}{R(y) + \rho\,e(y)},
\qquad \log_2\text{-odds} = \log_2\frac{C(y)}{E(y)} - \log_2\rho.$$

## 4. Which way the assumptions push

| departure from the simulation | effect on recall $R$ | effect on precision |
|---|---|---|
| sequencing errors, SNPs vs reference | lower (in expectation) | lower |
| origins outside the panel (TEs, unannotated copies) | none | lower (extra term in the denominator, absent from $E$) |
| paralog more strongly marked ($\rho > 1$) | none | lower |
| paralog less strongly marked ($\rho < 1$) | none | higher |
| real reads longer than the rounded-down $L$ | higher | higher (usually) |

So:

- Error-free $R(y)$ is an upper bound on real recall in expectation, provided the reference
  matches the sample and $L$ is not shorter than the real reads. Errors cannot help a read:
  a mismatch at a base shared by the copies leaves a tie intact, and one at a diverging base
  either favours the paralog or creates a tie; both also lower MAPQ. bowtie2's seed heuristics
  are not strictly monotonic, so this holds in expectation, not at every base. Rounding $L$
  down pushes the other way.
- $\mathrm{Prec}_{\mathrm{eq}}(y)$ is **not** a bound. It is precision at a neutral point
  ($\rho = 1$); real precision can be higher or lower depending on the latent signal at the
  paralogs. It is an upper bound only with respect to the error and off-panel terms.

## 5. Combining fragment and read lengths for a real run

Let $p_F$ be the run's fragment-length distribution (share of *fragments* of length $F$),
projected onto the grid. Fragments of length $F$ cover $F$ bases, so their share of coverage
is proportional to $p_F F$, and $N_F(x) \propto F$ already carries that factor. Counts should
therefore be summed before dividing:

$$R_p(x) = \frac{\sum_F p_F\,C_F(x)}{\sum_F p_F\,N_F(x)},
\qquad e_p(y) = \frac{\sum_F p_F\,E_F(y)}{\sum_F p_F\,N_F(y)}.$$

Equivalently $R_p(x) = \sum_F v_F\,R_F(x)$ with coverage weights
$v_F \propto p_F\,N_F(x) \propto p_F\,F$. The same applies to a mixture of read lengths.
Log ratios must not be averaged; recompute them from the combined counts.

Note: the analysis scripts behind `simulation_findings.md` weighted $R_F$ directly by $p_F$
(fragment weights, not coverage weights), which under-weights long fragments. The relative
comparisons there (cutoffs, grid trimming, step) are unlikely to change qualitatively, but
absolute run-weighted values should be recomputed with $v_F$.

## 6. Towards probabilities for real data

Recall comes straight from the simulation. Precision needs $\lambda$. For broad marks, a
natural route is to estimate $\lambda$ by inverting the kernel over windows $w$ (e.g. 1-5 kb),
where the system is small and well conditioned:

$$\mathbb{E}[O_w] = \sum_{w'} K_{w w'}\,\lambda_{w'},
\qquad K_{w w'} = \text{simulated fragments from } w' \text{ kept with coverage in } w,$$

with $K_{ww}$ built from $C$ and $K_{ww'}$ ($w \ne w'$) from $E_z$ (currently stored per locus
pair in `bleed_matrix`; per-window resolution would need a new output). Solve for
$\lambda \ge 0$ by non-negative least squares or EM with a Poisson likelihood for $O$; then
plug $\hat\lambda$ into section 3 for $\mathrm{Prec}(y)$, carrying the uncertainty of
$\hat\lambda$ into it. Windows where $R$ is near 0 carry no information about $\lambda$ and
should be reported as such rather than imputed.

## 7. Open points

- Error and SNP models, so that $R$ is a calibrated estimate rather than an upper bound.
- How to treat off-panel sources: simulate a background sample of the genome, or flag sites
  from `misplaced_sites` that receive fragments from the panel.
- Per-window ectopic matrices ($E_z$ at window resolution) for the inversion in section 6.
- Sensitivity of the inversion to the assumption that $\lambda$ is locally constant.
