import json
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CAMP = os.path.join(ROOT, "campaigns", "B01_PLUCKER_BIO_RESPONSE")
sys.path.insert(0, CAMP)

from b01_core import (PAIRS, TwoStateModel, brackets_common,
                      channels_two_state, gradients_from_stencil,
                      independent_brackets, pluecker_residual)

T0, C0, DT, DC = 330.0, 1.5, 0.5, 0.1
BASE = TwoStateModel()


def clean(T, c):
    return channels_two_state(BASE, T, c)


def test_common_gradient_is_algebraic_identity():
    u, v = gradients_from_stencil(clean, T0, C0, DT, DC)
    assert abs(pluecker_residual(brackets_common(u, v))) < 1e-18


def test_channels_span_chart():
    # DD-2 regression: gradients must not all be parallel (rank must be 2).
    u, v = gradients_from_stencil(clean, T0, C0, DT, DC)
    M = np.stack([u, v])
    s = np.linalg.svd(M, compute_uv=False)
    assert s[1] / s[0] > 1e-3, "channel gradients collapsed to rank 1"


def test_clean_independent_slots_near_zero():
    B = independent_brackets([clean] * 6, T0, C0, DT, DC)
    assert abs(pluecker_residual(B)) < 1e-15


def test_contamination_moves_residual():
    # DD-1 regression: gradient-affecting contamination must move P(B);
    # additive offsets must not be the falsification mechanism.
    def shifted(T, c):
        return channels_two_state(TwoStateModel(Tm=BASE.Tm + 2.0), T, c)
    fns = [clean, shifted, clean, shifted, clean, shifted]
    B = independent_brackets(fns, T0, C0, DT, DC)
    assert abs(pluecker_residual(B)) > 1e-8


def test_offset_contamination_is_invisible():
    # documents the DD-1 fact: per-slot constant offsets leave brackets
    # unchanged, so they are OUTSIDE this test's detection class.
    def offset(T, c):
        return channels_two_state(BASE, T, c) + 0.05
    B_clean = independent_brackets([clean] * 6, T0, C0, DT, DC)
    B_off = independent_brackets([offset] * 6, T0, C0, DT, DC)
    for k in B_clean:
        assert abs(B_clean[k] - B_off[k]) < 1e-15


def test_certificate_pin_reproduces():
    r = subprocess.run([sys.executable, os.path.join(CAMP, "run_b01.py")],
                       capture_output=True, text=True, check=True)
    with open(os.path.join(ROOT, "results", "B01_CERTIFICATE.json")) as f:
        cert = json.load(f)
    assert cert["_pin"] in r.stdout
    assert cert["B01A_COMMON_GRADIENT_CONTROL"]["status"] == "PASS_CONTROL"
    assert cert["B01B_CLEAN_INDEPENDENT_SLOTS"]["status"] == "NOT_FALSIFIED"
    assert cert["B01C_BATCH_AGING_SMALL"]["status"] == "NOT_FALSIFIED"
    assert cert["B01C_BATCH_AGING_LARGE"]["status"] == "FALSIFIED_MEASUREMENT_CONTRACT"
    assert cert["B01D_HYSTERESIS_BRANCH_MIX"]["status"] == "FALSIFIED_MEASUREMENT_CONTRACT"
