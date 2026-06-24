"""Background parenchymal enhancement (BPE) quantification.

Volume/intensity-based BPE over fibroglandular tissue follows:
    Wei D, et al. Fully automatic quantification of FGT and BPE. Med Phys 2021;48:238-252.
    Zhu Y, et al. Automated BPE as a biomarker of treatment response. Breast Cancer Res Treat 2026.
A signal-enhancement threshold defines "enhancing" parenchyma; categories use cohort cut-offs.
"""
from __future__ import annotations
import numpy as np

_EPS = 1e-6


def bpe_percent(phases, parenchyma_mask, enh_threshold: float = 0.0) -> float:
    """Mean percent enhancement of non-tumour parenchyma (early vs pre).

    parenchyma_mask should cover fibroglandular tissue with the tumour excluded. If only a breast
    mask is available, exclude the tumour mask first. Returns mean % enhancement over voxels whose
    enhancement exceeds `enh_threshold` (set e.g. to the 55% level of Zhu 2026 for a volume-based
    measure; 0.0 keeps all parenchyma).
    """
    pre, early = np.asarray(phases[0], float), np.asarray(phases[1], float)
    m = np.asarray(parenchyma_mask) > 0
    if m.sum() == 0:
        return float("nan")
    enh = (early[m] - pre[m]) / (pre[m] + _EPS) * 100.0
    sel = enh[enh > enh_threshold * 100.0]
    return float(sel.mean()) if sel.size else 0.0


def bpe_categories(values, q=(1 / 3, 2 / 3)):
    """Bin a vector of continuous BPE values into 3 radiologist-style categories by cohort tertiles.

    Returns (labels, (t1, t2)). minimal-mild / moderate / marked, matching the cohort-cutoff
    approach of Zhu 2026.
    """
    v = np.asarray(values, float)
    t1, t2 = np.nanquantile(v, q[0]), np.nanquantile(v, q[1])
    labels = np.where(v < t1, "minimal_mild", np.where(v < t2, "moderate", "marked"))
    return labels, (float(t1), float(t2))
