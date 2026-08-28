import json
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CAMP = os.path.join(ROOT, "campaigns", "B02_CUT_SQUARE_RESPONSE")
sys.path.insert(0, CAMP)

from b02_core import (BETA_DECLARED, TwoStateModel, cut_square_ledger,
                      declared_decoder, response_event)

T0, C0, DT, DC = 330.0, 1.5, 0.5, 0.1
BASE = TwoStateModel()
C = declared_decoder(BASE, T0, C0, DT, DC)


def test_exact_identity_T31_random_events():
    rng = np.random.default_rng(3)
    for _ in range(50):
        y = rng.normal(0, 1e-2, 4)
        led = cut_square_ledger(C, y)
        assert abs(led["identity_residual_T31_rel"]) < 1e-12


def test_ledger_T41_closes_and_parts_nonnegative():
    y = response_event(TwoStateModel(Tm=BASE.Tm + 1.0), T0, C0, DT, DC)
    led = cut_square_ledger(C, y)
    assert abs(led["ledger_residual_T41_rel"]) < 1e-12
    assert led["captured"] >= 0 and led["reserve"] >= 0 and led["transverse"] >= 0
    assert led["beta"] <= 1.0 + 1e-12


def test_clean_event_has_zero_transverse():
    y = response_event(BASE, T0, C0, DT, DC)
    led = cut_square_ledger(C, y)
    assert led["transverse_share"] < 1e-12
    assert abs(led["captured_share"] - BETA_DECLARED) < 1e-9


def test_whitening_contract_matters():
    # raw-coordinate wedging (no metric) must NOT reproduce the identity
    # terms; guards the theorem's Section 2 contract.
    y = response_event(TwoStateModel(Tm=BASE.Tm + 2.0), T0, C0, DT, DC)
    led = cut_square_ledger(C, y)
    raw = sum((C[i] * y[j] - C[j] * y[i]) ** 2
              for i in range(4) for j in range(i + 1, 4))
    assert abs(raw - led["transverse"]) > 1e-16


def test_instrument_vs_model_fingerprints_differ():
    y_model = response_event(TwoStateModel(Tm=BASE.Tm + 2.0), T0, C0, DT, DC)
    y_inst = response_event(BASE, T0, C0, DT, DC).copy()
    y_inst[3] *= 1.2
    assert cut_square_ledger(C, y_model)["dominant_pair"] != \
        cut_square_ledger(C, y_inst)["dominant_pair"]


def test_certificate_pin_reproduces():
    r = subprocess.run([sys.executable, os.path.join(CAMP, "run_b02.py")],
                       capture_output=True, text=True, check=True)
    with open(os.path.join(ROOT, "results", "B02_CERTIFICATE.json")) as f:
        cert = json.load(f)
    assert cert["_pin"] in r.stdout
    for k in ("B02A_EXACT_IDENTITY_T31", "B02B_CONTRACTIVE_LEDGER_T41",
              "B02C_CLEAN_DECODER_ADEQUACY",
              "B02D_CONTAMINATION_LOCALIZATION"):
        assert cert[k]["status"] == "PASS_CONTROL"
