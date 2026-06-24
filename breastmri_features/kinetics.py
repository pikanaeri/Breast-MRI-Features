"""Dynamic (kinetic) DCE-MRI descriptors inside a tumour mask.

Signal-enhancement-ratio and enhancement-curve descriptors follow:
    Agliozzo S, et al. CAD for DCE breast MRI of mass-like lesions. Med Phys 2012;39:1704-1715.
    Dalmis MU, et al. A CAD system for breast DCE-MRI. Med Phys 2016;43:6260-6273.

Inputs:  phases — sequence [pre, early, late], each [D, H, W]
         mask   — binary array [D, H, W]
"""
from __future__ import annotations
import numpy as np

_EPS = 1e-6


def _rim_ratio(enh, mask) -> float:
    """Peripheral / core mean-enhancement ratio (rim enhancement detector)."""
    from scipy import ndimage
    m = np.asarray(mask) > 0
    core = ndimage.binary_erosion(m, iterations=2)
    rim = m & ~core
    if rim.sum() == 0 or core.sum() == 0:
        return float("nan")
    return float(enh[rim].mean() / (enh[core].mean() + _EPS))


def signal_enhancement(phases, mask) -> dict:
    """Return kinetic descriptors averaged inside the mask.

    init_pct_mean    initial (early vs pre) percent enhancement
    delayed_pct_mean delayed-phase percent change (late vs early) — the washout/persistence slope
    ser_mean         signal-enhancement ratio (early-pre)/(late-pre)
    enh_cov          coefficient of variation of enhancement (internal heterogeneity)
    enh_cov_frac     fraction of voxels with meaningful enhancement
    rim_ratio        peripheral/core enhancement ratio (rim enhancement)
    """
    pre, early, late = (np.asarray(p, float) for p in phases)
    m = np.asarray(mask) > 0
    if m.sum() == 0:
        return {k: float("nan") for k in
                ("init_pct_mean", "delayed_pct_mean", "ser_mean", "enh_cov",
                 "enh_cov_frac", "rim_ratio")}
    pre_m, early_m, late_m = pre[m], early[m], late[m]
    init = (early_m - pre_m) / (pre_m + _EPS) * 100.0
    delayed = (late_m - early_m) / (early_m - pre_m + _EPS) * 100.0
    ser = (early_m - pre_m) / (late_m - pre_m + _EPS)
    enh = early_m - pre_m
    cov = float(np.std(enh) / (abs(np.mean(enh)) + _EPS))
    frac = float(np.mean(enh > 0.1 * (np.abs(pre_m).mean() + _EPS)))
    return {
        "init_pct_mean": float(np.mean(init)),
        "delayed_pct_mean": float(np.median(delayed)),
        "ser_mean": float(np.mean(ser)),
        "enh_cov": cov,
        "enh_cov_frac": frac,
        "rim_ratio": _rim_ratio(early - pre, m),
    }
