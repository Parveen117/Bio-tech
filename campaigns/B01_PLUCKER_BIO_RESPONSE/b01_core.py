"""B01 core: Pluecker bracket audit for biothermodynamic response channels.

Framework source
----------------
arXiv:2603.20773 (Geometric Completion of Thermodynamic Response), rank-2
area-bracket theorem: for any four smooth nondimensional response channels
on a two-dimensional equilibrium chart, the six pairwise area brackets
B_ij assemble into a skew 4x4 matrix of rank at most two, hence

    P(B) = B12*B34 - B13*B24 + B14*B23 = 0.

The audit discipline (common-gradient control vs independently estimated
brackets, delta-method sigma, five-sigma rejection) follows campaign T01
of the Thermodynamics-Reproducibility repository.

Bio chart
---------
The two-dimensional chart is (T, c): absolute temperature and a chemical
control coordinate (denaturant concentration for protein unfolding; for
nucleic-acid melts c may be log salt activity). Channels are smooth
nondimensional response functions on that chart.

Honest claim boundary (stated up front)
---------------------------------------
The Pluecker test falsifies the *single two-dimensional equilibrium
chart* description of a measurement ensemble. It does NOT distinguish
two-state from multi-state equilibria: any equilibrium ensemble whose
channels are single-valued smooth functions of (T, c) satisfies the
identity regardless of the number of states. What the test catches is
exactly the quality-control failure class: hidden third variables
(batch history, scan-rate memory, aggregation, drift) and
path-dependence (hysteresis), i.e. data that does not live on one chart.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass

import numpy as np

R_GAS = 8.31446261815324  # J / (mol K)

REJECTION_SIGMA = 5.0


# ---------------------------------------------------------------------------
# Two-state equilibrium generator (linear extrapolation model)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TwoStateModel:
    """Two-state folding equilibrium on the (T, c) chart.

    dG(T, c) = dH0 * (1 - T/Tm) + dCp * (T - Tm - T*ln(T/Tm)) - m * c
    theta(T, c) = 1 / (1 + exp(-dG / (R T)))   (folded fraction)
    """

    dH0: float = 250e3     # J/mol at Tm
    Tm: float = 333.15     # K
    dCp: float = 6.0e3     # J/(mol K)
    m: float = 4.0e3       # J/(mol M)

    def dG(self, T, c):
        T = np.asarray(T, dtype=float)
        return (self.dH0 * (1.0 - T / self.Tm)
                + self.dCp * (T - self.Tm - T * np.log(T / self.Tm))
                - self.m * np.asarray(c, dtype=float))

    def theta(self, T, c):
        return 1.0 / (1.0 + np.exp(-self.dG(T, c) / (R_GAS * np.asarray(T))))


def channels_two_state(model: TwoStateModel, T, c):
    """Four smooth nondimensional channels on the (T, c) chart.

    Chart-spanning by construction (design defect DD-2 fix): x2 depends on
    T beyond the stability scalar g, and x4 is a pure thermal channel, so
    the four gradients genuinely span the two chart directions. An earlier
    channel set factored entirely through the single scalar g(T, c); all
    brackets then vanished identically (rank <= 1) and the test had no
    falsification power. Recorded in README.
    """
    T = np.asarray(T, dtype=float)
    c = np.asarray(c, dtype=float)
    g = model.dG(T, c) / (R_GAS * T)
    th = 1.0 / (1.0 + np.exp(-g))
    dH_T = model.dH0 + model.dCp * (T - model.Tm)
    x1 = th                                        # folded fraction (CD/fluorescence)
    x2 = th * (1.0 - th) * dH_T / (R_GAS * T) * 0.01   # excess-response (DSC) proxy
    x3 = np.tanh(g / 2.0)                          # signed population balance
    x4 = np.log(T / 300.0)                         # nondimensional thermal channel
    return np.stack([x1, x2, x3, x4], axis=0)


# ---------------------------------------------------------------------------
# Brackets
# ---------------------------------------------------------------------------

def gradients_from_stencil(channel_fn, T0, c0, dT, dc):
    """Central-difference chart gradient of each channel at (T0, c0)."""
    Xp_T = channel_fn(T0 + dT, c0)
    Xm_T = channel_fn(T0 - dT, c0)
    Xp_c = channel_fn(T0, c0 + dc)
    Xm_c = channel_fn(T0, c0 - dc)
    u = (Xp_T - Xm_T) / (2.0 * dT)   # d/dT, shape (4,)
    v = (Xp_c - Xm_c) / (2.0 * dc)   # d/dc
    return u, v


PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def brackets_common(u, v):
    """All six brackets from ONE shared gradient pair (algebraic control)."""
    return {f"B{i+1}{j+1}": u[i] * v[j] - u[j] * v[i] for i, j in PAIRS}


def pluecker_residual(B):
    return B["B12"] * B["B34"] - B["B13"] * B["B24"] + B["B14"] * B["B23"]


def pluecker_scale(B):
    vals = np.array([B[k] for k in B], dtype=float)
    s = float(np.max(np.abs(vals)))
    return s * s if s > 0 else 1.0


def independent_brackets(slot_fns, T0, c0, dT, dc):
    """Six brackets, each estimated from its own disjoint slot.

    slot_fns : list of six channel functions f(T, c) -> (4,) array, one per
    bracket slot, representing disjoint measurement ensembles (different
    samples / publications / batches). Each slot contributes only its own
    bracket, so no shared gradient pair forces P(B)=0 algebraically.
    """
    B = {}
    for (i, j), fn in zip(PAIRS, slot_fns):
        u, v = gradients_from_stencil(fn, T0, c0, dT, dc)
        B[f"B{i+1}{j+1}"] = u[i] * v[j] - u[j] * v[i]
    return B


N_REPLICATE_AVG = 25   # declared replicate averaging per stencil node


def noise_standard_error(slot_fns, T0, c0, dT, dc, noise_sigma,
                         n_rep=300, seed=7):
    """Empirical SE of the RAW Pluecker residual under per-channel Gaussian
    measurement noise with declared replicate averaging. Raw (unnormalized)
    residual is used for the z statistic; the earlier per-draw
    max-bracket-squared normalization inflated the variance (recorded as
    design note in README)."""
    rng = np.random.default_rng(seed)
    eff = noise_sigma / math.sqrt(N_REPLICATE_AVG)
    res = []
    for _ in range(n_rep):
        B = {}
        for (i, j), fn in zip(PAIRS, slot_fns):
            def noisy(T, c, fn=fn):
                x = fn(T, c)
                return x + rng.normal(0.0, eff, size=x.shape)
            u, v = gradients_from_stencil(noisy, T0, c0, dT, dc)
            B[f"B{i+1}{j+1}"] = u[i] * v[j] - u[j] * v[i]
        res.append(pluecker_residual(B))
    res = np.asarray(res)
    return float(np.mean(res)), float(np.std(res, ddof=1))


def verdict(z):
    return "NOT_FALSIFIED" if abs(z) < REJECTION_SIGMA else "FALSIFIED_MEASUREMENT_CONTRACT"


def sha256_of(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()
