"""Minimal example: derive all features from a synthetic lesion (no data needed)."""
import numpy as np
from breastmri_features import derive_features, CITATIONS

# synthetic 3-phase DCE volume + a roughly spherical enhancing lesion
D = H = W = 48
zz, yy, xx = np.ogrid[:D, :H, :W]
mask = ((zz - 24) ** 2 + (yy - 24) ** 2 + (xx - 24) ** 2) <= 9 ** 2
pre = np.random.RandomState(0).rand(D, H, W) * 50 + 50
early = pre + mask * 120        # strong early enhancement inside the lesion
late = early - mask * 30        # washout in the delayed phase

feats = derive_features(mask.astype(np.uint8), [pre, early, late], spacing=(2.0, 0.8, 0.8))

for k in ["longest_diameter", "tumor_volume", "sphericity", "shape", "margin",
          "mass_vs_nme", "kinetic_curve", "internal_enhancement"]:
    conf = feats.get(f"{k}_gtconf")
    print(f"{k:24s} = {feats[k]!s:16s}" + (f"  (gt-conf {conf})" if conf is not None else ""))
print("\nprovenance:", feats["provenance"])
print("\nsphericity sources:", CITATIONS["sphericity"][0])
