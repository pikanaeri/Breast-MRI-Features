"""Rule-based mapping of measured features to BI-RADS / kinetic categories, with a per-label
ground-truth DERIVATION confidence (distance of the measured value from the decision threshold,
normalised to [0, 1]; 0 = on the threshold/ambiguous, 1 = unambiguous).

The thresholds are the only tunable parameters of the generator and live in one table. The
measurements being thresholded come from `shape` and `kinetics`; their methods are cited there.
BI-RADS-aligned morphological/kinetic descriptors follow:
    Gilhuijs KGA, Giger ML, Bick U. Med Phys 1998;25:1647-1654.
    Huang YH, et al. Comput Methods Programs Biomed 2013;112:508-517.
    Agliozzo S, et al. Med Phys 2012;39:1704-1715.
    Dalmis MU, et al. Med Phys 2016;43:6260-6273.
Categorisation thresholds are investigator-applied (see README).
"""
from __future__ import annotations

# ---- decision thresholds (single source of truth) ----------------------------------------------
THRESHOLDS = {
    "shape_sphericity": 0.82,     # sphericity >= -> round/oval
    "margin_spiculation": 0.45,   # spiculation index >= -> spiculated
    "margin_irregular": 1.15,     # boundary irregularity below -> circumscribed
    "nme_solidity": 0.50,         # solidity below (or >=3 foci) -> non-mass enhancement
    "curve_persistent": 10.0,     # delayed % >= -> persistent
    "curve_washout": -10.0,       # delayed % <= -> washout
    "init_slow": 50.0,            # initial % below -> slow
    "init_fast": 100.0,           # initial % above -> fast
    "rim_ratio": 1.30,            # peripheral/core ratio >= -> rim enhancement
    "hetero_cov": 0.60,           # enhancement CoV >= -> heterogeneous
}
# per-concept scale used to normalise the threshold margin into a confidence
_SCALE = {"shape": 0.12, "mass_vs_nme": 0.15, "margin_spic": 0.20, "margin_irreg": 0.30,
          "curve": 10.0, "init": 25.0, "rim": 0.30, "hetero": 0.25}

T = THRESHOLDS


def _conf(dist, scale):
    return round(float(min(1.0, abs(dist) / scale)), 3)


# ---- categorical maps + their derivation confidence ---------------------------------------------
def shape(sphericity):
    return "round_oval" if sphericity >= T["shape_sphericity"] else "irregular"

def shape_conf(sphericity):
    return _conf(sphericity - T["shape_sphericity"], _SCALE["shape"])


def margin(spiculation, irregularity):
    if spiculation >= T["margin_spiculation"]:
        return "spiculated"
    return "circumscribed" if irregularity <= T["margin_irregular"] else "irregular"

def margin_conf(spiculation, irregularity):
    if spiculation >= T["margin_spiculation"]:
        return _conf(spiculation - T["margin_spiculation"], _SCALE["margin_spic"])
    return _conf(irregularity - T["margin_irregular"], _SCALE["margin_irreg"])


def mass_vs_nme(solidity, n_foci):
    return "nme" if (solidity < T["nme_solidity"] or n_foci >= 3) else "mass"

def mass_vs_nme_conf(solidity, n_foci):
    return 1.0 if n_foci >= 3 else _conf(solidity - T["nme_solidity"], _SCALE["mass_vs_nme"])


def multifocality(n_foci):
    return "multifocal" if n_foci >= 2 else "unifocal"

def multifocality_conf(n_foci):
    return 1.0 if (n_foci <= 1 or n_foci >= 3) else 0.5


def kinetic_curve(delayed_pct):
    if delayed_pct >= T["curve_persistent"]:
        return "persistent"
    return "washout" if delayed_pct <= T["curve_washout"] else "plateau"

def kinetic_curve_conf(delayed_pct):
    return _conf(min(abs(delayed_pct - T["curve_persistent"]),
                     abs(delayed_pct - T["curve_washout"])), _SCALE["curve"])


def initial_enhancement_rate(init_pct):
    if init_pct < T["init_slow"]:
        return "slow"
    return "fast" if init_pct > T["init_fast"] else "medium"

def initial_enhancement_rate_conf(init_pct):
    return _conf(min(abs(init_pct - T["init_slow"]), abs(init_pct - T["init_fast"])), _SCALE["init"])


def internal_enhancement(rim_ratio, enh_cov):
    if rim_ratio >= T["rim_ratio"]:
        return "rim"
    return "heterogeneous" if enh_cov >= T["hetero_cov"] else "homogeneous"

def internal_enhancement_conf(rim_ratio, enh_cov):
    if rim_ratio >= T["rim_ratio"]:
        return _conf(rim_ratio - T["rim_ratio"], _SCALE["rim"])
    return _conf(enh_cov - T["hetero_cov"], _SCALE["hetero"])
