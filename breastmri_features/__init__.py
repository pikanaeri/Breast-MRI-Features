"""breastmri_features — reproducible, citation-backed derivation of quantitative and BI-RADS /
kinetic descriptors from breast DCE-MRI tumour segmentations.

    from breastmri_features import derive_features
    feats = derive_features(mask, [pre, early, late], spacing=(2.0, 0.8, 0.8))

See CITATIONS for the method references behind each feature.
"""
from .derive import derive_features
from . import shape, kinetics, bpe, birads
from .citations import CITATIONS

__version__ = "0.1.0"
__all__ = ["derive_features", "shape", "kinetics", "bpe", "birads", "CITATIONS"]
