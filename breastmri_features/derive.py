"""Top-level: derive every quantitative + BI-RADS/kinetic feature from one tumour mask and its
3-phase DCE acquisition, with a ground-truth derivation confidence on each categorical concept.

Each continuous feature can optionally be overridden by a vendor sheet value (for traceability),
recorded in `provenance`.

Method references per feature are listed in `breastmri_features.CITATIONS` and in the README
References table (IBSI/Zwanenburg 2020; PyRadiomics/van Griethuysen 2017; RECIST/Eisenhauer 2009;
scikit-image/van der Walt 2014; Gilhuijs 1998; Huang 2013; Agliozzo 2012; Dalmis 2016; Wei 2021;
Zhu 2026).
"""
from __future__ import annotations
import numpy as np
from . import shape as _shape, kinetics as _kin, birads as _b


def derive_features(mask, phases, spacing, *, sheet: dict | None = None,
                    clean_components: int = 2) -> dict:
    """Derive all features.

    Parameters
    ----------
    mask     : binary array [D, H, W] (axial along axis 0)
    phases   : sequence [pre, early, late], each [D, H, W]
    spacing  : (z, y, x) voxel size in mm
    sheet    : optional {longest_diameter, tumor_volume, sphericity, bpe} vendor values to PREFER
               for the continuous concepts (kept traceable); algorithm used where absent.
    clean_components : keep the N largest connected components first (speckle removal)

    Returns a flat dict of features, `<concept>_gtconf` derivation confidences, and `provenance`.
    """
    mask = _shape.keep_largest_components(mask, clean_components)
    sp = tuple(map(float, spacing))
    sheet = sheet or {}
    prov = {}

    # continuous (algorithmic, vendor-preferred where given)
    alg = {
        "longest_diameter_mm": _shape.max_2d_diameter_mm(mask, sp),
        "tumor_volume_mm3": _shape.voxel_volume_mm3(mask, sp),
        "sphericity": _shape.mesh_sphericity(mask, sp),
    }
    cont = {}
    for k, alias in [("longest_diameter_mm", "longest_diameter"),
                     ("tumor_volume_mm3", "tumor_volume"), ("sphericity", "sphericity")]:
        if alias in sheet and sheet[alias] is not None and not _isnan(sheet[alias]):
            cont[alias] = float(sheet[alias]); prov[alias] = "vendor_sheet"
        else:
            cont[alias] = alg[k]; prov[alias] = "algorithm"

    # intermediate shape descriptors that FEED the categorical BI-RADS concepts below
    # (not standalone concepts themselves; shape.elongation_flatness is available in the library
    #  but is not part of the derived concept set, so it is not computed here).
    solidity = _shape.solidity(mask)                     # -> mass_vs_nme
    spic = _shape.spiculation_index(mask, sp)            # -> margin (boundary descriptor)
    n_foci = _shape.n_components(mask)                   # -> mass_vs_nme, multifocality

    # kinetics
    k = _kin.signal_enhancement(phases, mask)

    out = dict(cont)
    out.update({
        "solidity": solidity, "spiculation": spic, "n_foci": n_foci,
        "shape": _b.shape(cont["sphericity"]),
        "margin": _b.margin(spic),
        "mass_vs_nme": _b.mass_vs_nme(solidity, n_foci),
        "multifocality": _b.multifocality(n_foci),
        "kinetic_curve": _b.kinetic_curve(k["delayed_pct_mean"]),
        "initial_enhancement_rate": _b.initial_enhancement_rate(k["init_pct_mean"]),
        "internal_enhancement": _b.internal_enhancement(k["rim_ratio"], k["enh_cov"]),
        # ground-truth derivation confidences
        "shape_gtconf": _b.shape_conf(cont["sphericity"]),
        "margin_gtconf": _b.margin_conf(spic),
        "mass_vs_nme_gtconf": _b.mass_vs_nme_conf(solidity, n_foci),
        "multifocality_gtconf": _b.multifocality_conf(n_foci),
        "kinetic_curve_gtconf": _b.kinetic_curve_conf(k["delayed_pct_mean"]),
        "initial_enhancement_rate_gtconf": _b.initial_enhancement_rate_conf(k["init_pct_mean"]),
        "internal_enhancement_gtconf": _b.internal_enhancement_conf(k["rim_ratio"], k["enh_cov"]),
        "provenance": prov,
    })
    out.update({f"kin_{kk}": vv for kk, vv in k.items()})
    return out


def _isnan(x):
    try:
        return np.isnan(float(x))
    except (TypeError, ValueError):
        return True
