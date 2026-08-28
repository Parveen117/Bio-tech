"""B02 core: exact cut-square decomposition of biothermodynamic response.

Framework source (consumed, not rederived)
------------------------------------------
Recognition-Kernel-Framework, theorum/thermodynamics/10
(Thermodynamic Cut-Square Response Decomposition Theorem):

  T3.1 (exact identity)   beta * |y|_H^2 - |L_c(y)|^2
                            = sum_{i<j} |C_i Y_j - C_j Y_i|^2
  T4.1 (contractive form)  for beta <= 1,
                          |y|_H^2 = |L_c(y)|^2
                                    + (1 - beta) |y|_H^2
                                    + sum_{i<j} |C_i Y_j - C_j Y_i|^2,
                          every term nonnegative,
  with C = H^{1/2} c, Y = H^{1/2} y, L_c(y) = <c, y>_H,
  beta = |c|_H^2.

Whitening contract from the theorem's Section 2, obeyed here: raw
response coordinates with unlike physical units are never wedged;
all pairwise squares are computed in metric-whitened coordinates
after the declared positive metric H.

Bio carrier
-----------
Events y are local response vectors of the four B01 channels to a
declared protocol step (dT, dc) at a chart point. The declared decoder
c is the clean two-state response direction for the same protocol,
scaled to beta = 0.95 < 1 so the contractive reserve term of T4.1 is
active. H is the declared diagonal noise-precision metric.

What B02 adds beyond B01: the six transverse squares LOCALIZE a
contract violation to channel pairs. B01 says the single-chart claim is
broken; B02 says where the response energy that the two-state decoder
cannot capture is sitting.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "B01_PLUCKER_BIO_RESPONSE"))
from b01_core import R_GAS, TwoStateModel, channels_two_state  # noqa: E402

CHANNEL_NOISE = np.array([1e-4, 1e-4, 1e-4, 1e-4])
BETA_DECLARED = 0.95
PAIR_LABELS = ["12", "13", "14", "23", "24", "34"]


def metric_H():
    """Declared positive response metric: diagonal noise precision."""
    return np.diag(1.0 / CHANNEL_NOISE ** 2)


def response_event(model: TwoStateModel, T0, c0, dT, dc, wT=1.0, wc=1.0):
    """Local response of the four channels to the declared protocol step
    (wT * dT, wc * dc) at chart point (T0, c0), by central differences."""
    xp = channels_two_state(model, T0 + dT, c0)
    xm = channels_two_state(model, T0 - dT, c0)
    yp = channels_two_state(model, T0, c0 + dc)
    ym = channels_two_state(model, T0, c0 - dc)
    return wT * (xp - xm) / 2.0 + wc * (yp - ym) / 2.0


def declared_decoder(base: TwoStateModel, T0, c0, dT, dc):
    """Clean two-state response direction, scaled to beta = BETA_DECLARED."""
    H = metric_H()
    c_raw = response_event(base, T0, c0, dT, dc)
    nrm = math.sqrt(float(c_raw @ H @ c_raw))
    return c_raw * (math.sqrt(BETA_DECLARED) / nrm)


def cut_square_ledger(c, y):
    """Exact T3.1 / T4.1 accounting for one event.

    Returns the three nonnegative parts of |y|_H^2 (captured, reserve,
    transverse), the six pairwise transverse squares in whitened
    coordinates, and the two identity residuals which must vanish to
    machine precision.
    """
    H = metric_H()
    Hs = np.sqrt(H)          # H is diagonal
    C = Hs @ c
    Y = Hs @ y
    beta = float(C @ C)
    total = float(Y @ Y)
    L = float(C @ Y)
    captured = L * L
    reserve = (1.0 - beta) * total
    squares = {}
    for lbl, (i, j) in zip(PAIR_LABELS,
                           [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]):
        squares[lbl] = float((C[i] * Y[j] - C[j] * Y[i]) ** 2)
    transverse = float(sum(squares.values()))
    resid_T31 = beta * total - captured - transverse
    resid_T41 = total - captured - reserve - transverse
    scale = max(beta * total, 1.0)
    return {
        "beta": beta,
        "total_response_energy": total,
        "captured": captured,
        "reserve": reserve,
        "transverse": transverse,
        "transverse_pair_squares": squares,
        "identity_residual_T31_rel": resid_T31 / scale,
        "ledger_residual_T41_rel": resid_T41 / max(total, 1.0),
        "captured_share": captured / total if total > 0 else 0.0,
        "transverse_share": transverse / total if total > 0 else 0.0,
        "dominant_pair": max(squares, key=squares.get),
    }


def sha256_of(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()
