"""Run rungs B02A-B02D and write the pinned certificate.

B02A  exact identity audit (T3.1) over an event grid  -> PASS_CONTROL
B02B  contractive ledger closure (T4.1)               -> PASS_CONTROL
B02C  decoder adequacy on clean two-state events      -> transverse ~ 0
B02D  contamination localization                      -> transverse jumps,
      dominant channel pair reported per contamination
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from b02_core import (BETA_DECLARED, TwoStateModel, cut_square_ledger,
                      declared_decoder, response_event, sha256_of)

T0, C0, DT, DC = 330.0, 1.5, 0.5, 0.1
base = TwoStateModel()
c = declared_decoder(base, T0, C0, DT, DC)

results = {}

# event grid: protocol weights and nearby chart points
grid = [(T0 + a, C0 + b, wT, wc)
        for a in (-2.0, 0.0, 2.0)
        for b in (-0.3, 0.0, 0.3)
        for (wT, wc) in ((1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, -0.5))]

# ---- B02A + B02B: identity and ledger closure over the grid ----
maxT31, maxT41 = 0.0, 0.0
for (Tg, cg, wT, wc) in grid:
    y = response_event(base, Tg, cg, DT, DC, wT, wc)
    led = cut_square_ledger(c, y)
    maxT31 = max(maxT31, abs(led["identity_residual_T31_rel"]))
    maxT41 = max(maxT41, abs(led["ledger_residual_T41_rel"]))
results["B02A_EXACT_IDENTITY_T31"] = {
    "events": len(grid), "max_relative_residual": maxT31,
    "status": "PASS_CONTROL" if maxT31 < 1e-12 else "FAIL_CONTROL",
}
results["B02B_CONTRACTIVE_LEDGER_T41"] = {
    "events": len(grid), "beta_declared": BETA_DECLARED,
    "max_relative_residual": maxT41,
    "status": "PASS_CONTROL" if maxT41 < 1e-12 else "FAIL_CONTROL",
}

# ---- B02C: decoder adequacy on clean events (same protocol as decoder) ----
y_clean = response_event(base, T0, C0, DT, DC)
led = cut_square_ledger(c, y_clean)
results["B02C_CLEAN_DECODER_ADEQUACY"] = {
    "captured_share": led["captured_share"],
    "transverse_share": led["transverse_share"],
    "note": "decoder is the clean response direction scaled to beta="
            f"{BETA_DECLARED}; captured share equals beta and transverse "
            "vanishes exactly on the aligned event",
    "status": "PASS_CONTROL"
    if led["transverse_share"] < 1e-12
    and abs(led["captured_share"] - BETA_DECLARED) < 1e-9
    else "FAIL_CONTROL",
}

# ---- B02D: contamination localization ----
contaminations = {
    "BATCH_AGING_2K": TwoStateModel(Tm=base.Tm + 2.0),
    "COOLING_BRANCH_2p5K": TwoStateModel(Tm=base.Tm - 2.5),
    "M_VALUE_DRIFT_20PC": TwoStateModel(m=base.m * 1.2),
}
loc = {}
events = {name: response_event(mdl, T0, C0, DT, DC)
          for name, mdl in contaminations.items()}
# instrument-side contamination: 20% gain error in the thermal readout
# channel (channel 4) -- a measurement distortion, not a model drift
y_gain = response_event(base, T0, C0, DT, DC).copy()
y_gain[3] *= 1.2
events["THERMAL_READOUT_GAIN_20PC"] = y_gain
for name, y in events.items():
    led = cut_square_ledger(c, y)
    loc[name] = {
        "captured_share": led["captured_share"],
        "transverse_share": led["transverse_share"],
        "dominant_pair": led["dominant_pair"],
        "pair_squares_normalized": {
            k: v / led["transverse"] for k, v in
            led["transverse_pair_squares"].items()
        } if led["transverse"] > 0 else {},
    }
results["B02D_CONTAMINATION_LOCALIZATION"] = {
    "cases": loc,
    "clean_transverse_share": results["B02C_CLEAN_DECODER_ADEQUACY"]["transverse_share"],
    "status": "PASS_CONTROL"
    if all(v["transverse_share"] > 1e-6 for v in loc.values())
    else "FAIL_CONTROL",
    "note": "every contamination moves response energy into the transverse "
            "squares; the dominant pair localizes where the two-state "
            "decoder fails. Channel key: 1 folded fraction, 2 excess "
            "response, 3 population balance, 4 thermal. Honest finding: "
            "the three model-level drifts (Tm aging, cooling branch, "
            "m-value) share one fingerprint (dominant pair 23) because "
            "they all act through the stability scalar; localization "
            "discriminates contamination CLASSES, not individual "
            "model-level drifts. The instrument-side thermal-readout "
            "gain case shows a distinct fingerprint on the thermal "
            "pairs.",
}

results["_meta"] = {
    "framework_source": "Recognition-Kernel-Framework "
                        "theorum/thermodynamics/10 (T3.1, T4.1)",
    "whitening_contract": "pairwise squares computed only after declared "
                          "diagonal noise-precision metric, per theorem "
                          "Section 2",
    "chart_point": {"T_K": T0, "c_M": C0},
    "protocol_step": {"dT_K": DT, "dc_M": DC},
    "beta_declared": BETA_DECLARED,
}
results["_pin"] = sha256_of({k: v for k, v in results.items() if k != "_pin"})

out = os.path.join(HERE, "..", "..", "results")
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, "B02_CERTIFICATE.json"), "w") as f:
    json.dump(results, f, indent=2, sort_keys=True)

for k, v in results.items():
    if k.startswith("_"):
        continue
    print(k, "->", v["status"])
    if k == "B02D_CONTAMINATION_LOCALIZATION":
        for n, d in v["cases"].items():
            print(f"   {n}: transverse={d['transverse_share']:.3e} "
                  f"dominant_pair={d['dominant_pair']}")
print("pin:", results["_pin"])
