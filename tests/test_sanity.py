"""Sanity tests: known shapes give expected feature values."""
import numpy as np
from breastmri_features import shape, derive_features


def _sphere(r=12, n=40):
    z, y, x = np.ogrid[:n, :n, :n]
    return (((z - n // 2) ** 2 + (y - n // 2) ** 2 + (x - n // 2) ** 2) <= r ** 2).astype(np.uint8)


def test_sphere_sphericity_near_one():
    s = shape.mesh_sphericity(_sphere(), (1.0, 1.0, 1.0))
    assert 0.85 <= s <= 1.05, s            # marching-cubes sphere ~0.9-1.0


def test_volume_matches_voxel_count():
    m = _sphere()
    assert shape.voxel_volume_mm3(m, (2.0, 0.5, 0.5)) == m.sum() * 0.5


def test_diameter_positive_and_3d_ge_2d():
    m = _sphere()
    sp = (2.0, 0.8, 0.8)
    d2, d3 = shape.max_2d_diameter_mm(m, sp), shape.max_3d_diameter_mm(m, sp)
    assert d2 > 0 and d3 >= d2 - 1e-6


def test_derive_returns_all_concepts():
    m = _sphere()
    pre = np.ones((40, 40, 40)) * 80.0
    feats = derive_features(m, [pre, pre + m * 100, pre + m * 70], spacing=(2.0, 0.8, 0.8))
    for key in ["shape", "margin", "mass_vs_nme", "multifocality", "kinetic_curve",
                "initial_enhancement_rate", "internal_enhancement"]:
        assert key in feats and feats[f"{key}_gtconf"] is not None


def test_components_filter():
    m = _sphere(6, 40)
    m[2:4, 2:4, 2:4] = 1                    # add a tiny speckle
    assert shape.n_components(shape.keep_largest_components(m, 1)) == 1


if __name__ == "__main__":
    for fn in [test_sphere_sphericity_near_one, test_volume_matches_voxel_count,
               test_diameter_positive_and_3d_ge_2d, test_derive_returns_all_concepts,
               test_components_filter]:
        fn(); print("PASS", fn.__name__)
