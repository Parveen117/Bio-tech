"""Run rungs B01A-B01C and write the pinned certificate.

B01A  common-gradient algebraic control              -> PASS_CONTROL
B01B  clean disjoint slots (independent, same chart) -> NOT_FALSIFIED
B01C  monotonic batch aging across slots (Tm drift)  -> small: NOT_FALSIFIED
                                                        large: FALSIFIED
B01D  heating/cooling branch mix (hysteresis)        -> FALSIFIED

Statistic: raw Pluecker residual P(B); SE from 300 noise draws at the
declared channel noise (1e-4 per channel, 25 replicate averages per
stencil node); rejection at 5 sigma; bootstrap confirmation for every
rung whose point estimate rejects.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from b01_core import (PAIRS, REJECTION_SIGMA, TwoStateModel, brackets_common,
                      channels_two_state, gradients_from_stencil,
                      independent_brackets, noise_standard_error,
                      pluecker_residual, sha256_of, verdict)

T0, C0, DT, DC = 330.0, 1.5, 0.5, 0.1
NOISE = 1e-4

results = {}
base = TwoStateModel()


def clean(T, c):
    return channels_two_state(base, T, c)


# ---- B01A: common-gradient algebraic control ----
u, v = gradients_from_stencil(clean, T0, C0, DT, DC)
rA = pluecker_residual(brackets_common(u, v))
results["B01A_COMMON_GRADIENT_CONTROL"] = {
    "raw_residual": rA,
    "status": "PASS_CONTROL" if abs(rA) < 1e-18 else "FAIL_CONTROL",
    "note": "Six brackets from one shared gradient pair: P(B)=0 is an "
            "algebraic identity here; implementation control only.",
}

# ---- B01B: clean independent slots ----
fns = [clean] * 6
r = pluecker_residual(independent_brackets(fns, T0, C0, DT, DC))
mu, se = noise_standard_error(fns, T0, C0, DT, DC, NOISE)
results["B01B_CLEAN_INDEPENDENT_SLOTS"] = {
    "raw_residual": r, "noise_standard_error": se,
    "z_score": r / se, "rejection_threshold": REJECTION_SIGMA,
    "status": verdict(r / se),
}


# ---- B01C: monotonic batch aging (hidden third variable) ----
def aging_slots(amp_K):
    hs = np.linspace(0.0, 1.0, 6) * amp_K
    return [(lambda T, c, h=h:
             channels_two_state(TwoStateModel(Tm=base.Tm + h), T, c))
            for h in hs]


for label, amp in (("SMALL", 0.01), ("LARGE", 2.0)):
    fns = aging_slots(amp)
    r = pluecker_residual(independent_brackets(fns, T0, C0, DT, DC))
    mu, se = noise_standard_error(fns, T0, C0, DT, DC, NOISE)
    results[f"B01C_BATCH_AGING_{label}"] = {
        "aging_amplitude_K": amp,
        "raw_residual": r, "noise_standard_error": se,
        "z_score": r / se, "rejection_threshold": REJECTION_SIGMA,
        "status": verdict(r / se),
    }


# ---- B01D: hysteresis branch mix ----
def cooling(T, c):
    return channels_two_state(TwoStateModel(Tm=base.Tm - 2.5), T, c)


fns = [clean, cooling, clean, cooling, clean, cooling]
r = pluecker_residual(independent_brackets(fns, T0, C0, DT, DC))
mu, se = noise_standard_error(fns, T0, C0, DT, DC, NOISE)
results["B01D_HYSTERESIS_BRANCH_MIX"] = {
    "branch_Tm_shift_K": 2.5,
    "raw_residual": r, "noise_standard_error": se,
    "z_score": r / se, "rejection_threshold": REJECTION_SIGMA,
    "status": verdict(r / se),
}

# ---- bootstrap confirmation for rejecting rungs ----
rng = np.random.default_rng(11)
for key, fns_b in (("B01C_BATCH_AGING_LARGE", aging_slots(2.0)),
                   ("B01D_HYSTERESIS_BRANCH_MIX",
                    [clean, cooling, clean, cooling, clean, cooling])):
    if results[key]["status"] != "FALSIFIED_MEASUREMENT_CONTRACT":
        continue
    draws = []
    eff = NOISE / np.sqrt(25)
    for _ in range(2000):
        B = {}
        for (i, j), fn in zip(PAIRS, fns_b):
            def noisy(T, c, fn=fn):
                x = fn(T, c)
                return x + rng.normal(0.0, eff, size=x.shape)
            uu, vv = gradients_from_stencil(noisy, T0, C0, DT, DC)
            B[f"B{i+1}{j+1}"] = uu[i] * vv[j] - uu[j] * vv[i]
        draws.append(pluecker_residual(B))
    lo, hi = np.percentile(draws, [0.05, 99.95])
    results[key]["bootstrap"] = {
        "n_draws": 2000, "p999_interval": [float(lo), float(hi)],
        "zero_inside_99p9": bool(lo <= 0.0 <= hi),
    }

results["_meta"] = {
    "chart": "(T, c): temperature x chemical control coordinate",
    "chart_point": {"T_K": T0, "c_M": C0},
    "stencil": {"dT_K": DT, "dc_M": DC},
    "channel_noise_sigma": NOISE, "replicate_averages": 25,
    "framework_source": "arXiv:2603.20773 rank-2 area-bracket theorem",
    "audit_discipline": "Thermodynamics-Reproducibility campaign T01",
}
results["_pin"] = sha256_of({k: v for k, v in results.items() if k != "_pin"})

out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "..", "results")
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, "B01_CERTIFICATE.json"), "w") as f:
    json.dump(results, f, indent=2, sort_keys=True)

for k, vr in results.items():
    if k.startswith("_"):
        continue
    line = f"{k}: {vr['status']}"
    if "z_score" in vr:
        line += f"  z={vr['z_score']:+.1f}"
    if "bootstrap" in vr:
        line += f"  zero_in_99.9%CI={vr['bootstrap']['zero_inside_99p9']}"
    print(line)
print("pin:", results["_pin"])
