# Bio-tech

<!-- patent-status-notice:start -->
> **Patent pending — rights reserved.** The repository owner confirms filed PCT and/or provisional patent applications in the related invention portfolio. This notice applies to subject matter supported by the applicable pending filings. Public disclosure is not intended to waive patent rights. Existing licences and permissions remain effective. See [PATENT_NOTICE.md](PATENT_NOTICE.md).
<!-- patent-status-notice:end -->

> **Research record:** [Dated snapshot, content hashes and archive links](RESEARCH_PROVENANCE.md) · [Canonical repository governance](GOVERNANCE.md). Cite the exact version used.


Biothermodynamic falsification campaigns built on the rank-2 area-bracket
theorem of **Geometric Completion of Thermodynamic Response**
(arXiv:2603.20773) and the audit discipline of the
[Thermodynamics-Reproducibility](https://github.com/Parveen117/Thermodynamics-Reproducibility)
campaign T01.

## First campaign

`B01_PLUCKER_BIO_RESPONSE`

For any four smooth nondimensional response channels on a two-dimensional
equilibrium chart, the six pairwise area brackets obey

```
P(B) = B12*B34 - B13*B24 + B14*B23 = 0.
```

The bio chart is `(T, c)`: temperature times a chemical control coordinate
(denaturant concentration for protein unfolding; log salt activity for
nucleic-acid melts). Channels are the responses laboratories already
measure: folded fraction (CD/fluorescence), excess-response (DSC proxy),
population balance, and a pure thermal channel.

## What the test detects — and what it does not

The Pluecker audit falsifies the claim that a measurement ensemble lives
on **one** two-dimensional equilibrium chart. It therefore detects the
quality-control failure class:

- hidden third variables — batch aging, scan-history, aggregation, drift
  across the disjoint measurement slots;
- path dependence — heating and cooling branches mixed under a single
  equilibrium description.

It does **not** distinguish two-state from multi-state equilibria: any
equilibrium ensemble whose channels are single-valued smooth functions of
`(T, c)` satisfies the identity regardless of the number of states. It
also does **not** see per-slot constant offsets (they vanish under the
chart gradient); offset-type miscalibration is outside this detection
class and is stated as such.

## Completed rung ladder

| Rung | Design | Result |
| --- | --- | --- |
| B01A | six brackets from one shared gradient pair (algebraic control) | `PASS_CONTROL` |
| B01B | six disjoint slots, identical clean chart | `NOT_FALSIFIED` (z = 0.0) |
| B01C | monotonic batch aging across slots, Tm drift 0.01 K | `NOT_FALSIFIED` (z = 0.0) |
| B01C | monotonic batch aging across slots, Tm drift 2.0 K | `FALSIFIED_MEASUREMENT_CONTRACT` (z = -7.7; zero outside 99.9% bootstrap interval) |
| B01D | heating/cooling branch mix, 2.5 K hysteresis | `FALSIFIED_MEASUREMENT_CONTRACT` (z = +14.7; zero outside 99.9% bootstrap interval) |

Statistic: raw Pluecker residual; SE from 300 noise draws at the declared
measurement contract (per-channel Gaussian noise 1e-4, 25 replicate
averages per stencil node); rejection threshold 5 sigma; 2,000-draw
bootstrap confirmation for every rejecting rung. Certificate:
`results/B01_CERTIFICATE.json` (pinned; CI regenerates and diffs it).

## Design defects found and fixed (recorded, not hidden)

- **DD-1 — offset contamination has no teeth.** The first hidden-variable
  rung injected per-slot additive constants; constants vanish under the
  chart gradient, so every bracket was untouched and the rung could not
  falsify anything. Fix: contamination must move the response *surface*
  (Tm drift). A regression test pins the fact that offsets remain
  invisible, because that is a property of the theorem, not a bug.
- **DD-2 — single-scalar channels collapse to rank 1.** The first channel
  set factored entirely through the stability scalar `g(T, c)`; all
  gradients were parallel, all brackets vanished identically, and the
  audit was vacuous. Fix: a chart-spanning channel set (the DSC-type
  excess channel depends on `T` beyond `g`; one pure thermal channel).
  A rank regression test guards this.
- **Design note — normalization variance.** Normalizing the residual by
  the per-draw max bracket squared inflated the noise variance and hid
  real signals; the statistic now uses the raw residual with a fixed
  measurement contract.

## Second campaign

`B02_CUT_SQUARE_RESPONSE`

Consumes the Thermodynamic Cut-Square Response Decomposition Theorem
(Recognition-Kernel-Framework, `theorum/thermodynamics/10`, T3.1/T4.1):
for a declared decoder `c` with `beta = |c|_H^2 <= 1`, every response
event `y` splits exactly as

```
|y|_H^2 = |<c, y>_H|^2  +  (1 - beta)|y|_H^2  +  sum of six pairwise
          (captured)       (reserve)             transverse squares
```

with the theorem's whitening contract obeyed: pairwise squares are
computed only after the declared diagonal noise-precision metric —
raw channels with unlike units are never wedged.

The bio decoder is the clean two-state response direction (beta = 0.95
declared); events are channel responses to a declared protocol step.
B02 upgrades B01 from detection to localization: B01 says the
single-chart contract is broken; B02's six transverse squares say where
the response energy that the two-state decoder cannot capture is
sitting, channel pair by channel pair.

| Rung | Design | Result |
| --- | --- | --- |
| B02A | exact T3.1 identity over a 36-event grid | `PASS_CONTROL` (max rel. residual < 1e-12) |
| B02B | T4.1 ledger closure, beta = 0.95 | `PASS_CONTROL` |
| B02C | clean aligned event | transverse share < 1e-12, captured share = beta |
| B02D | contamination localization | model drifts (Tm aging, cooling branch, m-value): dominant pair 23; instrument thermal-readout gain: dominant pair 34 |

Honest finding recorded in the certificate: the three model-level
drifts share one transverse fingerprint because they all act through
the stability scalar — localization discriminates contamination
CLASSES (model-side vs instrument-side), not individual model drifts.
Certificate: `results/B02_CERTIFICATE.json` (pinned; CI diffs it).

## Next rung — B01E experimental (open)

Data contract for the first real-data rung, declared before any data is
fitted:

1. Public thermal-unfolding or melt datasets providing at least two
   channels over a genuine `(T, c)` grid (temperature x denaturant, or
   temperature x salt), e.g. paired CD + DSC unfolding studies or qPCR
   melt panels with salt titration.
2. Six disjoint source slots (different publications, samples, or
   batches), one bracket per slot; no shared gradient pair.
3. Declared per-channel uncertainty from the source; 5 sigma threshold
   fixed in advance; bootstrap audit for any rejection.
4. Outcomes are reported whichever way they fall: `NOT_FALSIFIED` is a
   consistency certificate for the ensemble; `FALSIFIED_MEASUREMENT_CONTRACT`
   localizes a hidden-variable or path-dependence failure — the
   quality-control product claim.

## Claim boundary

`NOT_FALSIFIED` on synthetic rungs establishes only that the audit is
implemented correctly and has demonstrated falsification power at
declared contamination levels. No claim about any real biological
dataset is made until B01E is executed against its pre-declared
contract. No medical or diagnostic claim is made anywhere in this
repository.

## Reproduce

```bash
python -m pip install numpy pytest
python -m pytest tests -q
python campaigns/B01_PLUCKER_BIO_RESPONSE/run_b01.py
```
