# breastmri-features

Reproducible, **citation-backed** derivation of quantitative and BI-RADS / kinetic descriptors from
breast DCE-MRI tumour segmentations. Every feature is computed with a published, standardized
method, and each categorical (BI-RADS) label carries a **ground-truth derivation confidence**.

Pure Python — only `numpy`, `scipy`, `scikit-image`.

## Install

```bash
pip install -e .            # from this folder
# or: pip install numpy scipy scikit-image  &&  add this folder to PYTHONPATH
```

## Use

```python
from breastmri_features import derive_features

# mask: [D,H,W] binary tumour mask;  phases: [pre, early, late] DCE volumes;  spacing: (z,y,x) mm
feats = derive_features(mask, [pre, early, late], spacing=(2.0, 0.8, 0.8))

feats["sphericity"]            # IBSI mesh sphericity
feats["shape"]                 # 'round_oval' | 'irregular'
feats["shape_gtconf"]          # derivation confidence in [0,1]
feats["provenance"]            # which continuous features came from a vendor sheet vs the algorithm
```

Pass `sheet={"sphericity": 0.21, "tumor_volume": 28.0, ...}` to **prefer vendor values** for the
continuous features where available (kept traceable); the algorithm fills in the rest.

```bash
python examples/example_usage.py     # runs on a synthetic lesion, no data needed
pytest tests/                        # sanity tests
```

## Where each variable comes from (code → method)

These are the **12 imaging concepts** used as model variables (the canonical set):

| Concept | Function | Method / citation |
|---|---|---|
| Tumor volume | `shape.voxel_volume_mm3` | IBSI volume (voxel counting) — Zwanenburg 2020 |
| Sphericity | `shape.mesh_sphericity` | IBSI mesh sphericity (marching cubes) — Zwanenburg 2020; PyRadiomics (van Griethuysen 2017) |
| Longest diameter | `shape.max_2d_diameter_mm` | IBSI max 2-D diameter (slice) / RECIST 1.1 — Eisenhauer 2009 |
| Shape | `birads.shape` (sphericity threshold) | BI-RADS round/oval vs irregular (investigator rule) |
| Margin | `birads.margin` (spiculation + irregularity) | radial-gradient analysis — Gilhuijs 1998; Huang 2013 |
| Mass vs. NME | `birads.mass_vs_nme` (solidity, foci) | morphological CAD — Agliozzo 2012; Huang 2013 |
| Multifocality | `birads.multifocality` (`shape.n_components`) | connected components |
| Kinetic curve | `kinetics.signal_enhancement` → `birads.kinetic_curve` | SER curve descriptors — Agliozzo 2012; Dalmış 2016 |
| Initial enhancement rate | `kinetics.signal_enhancement` → `birads.initial_enhancement_rate` | early-phase enhancement — Agliozzo 2012; Dalmış 2016 |
| Internal enhancement | `kinetics.signal_enhancement` → `birads.internal_enhancement` | rim ratio + enhancement CoV — Agliozzo 2012; Huang 2013 |
| BPE | `bpe.bpe_percent` | volume/intensity BPE — Wei 2021 |
| BPE category (3-level) | `bpe.bpe_categories` | cohort tertiles of BPE — Wei 2021; ACR BI-RADS |

**Intermediate measurements** (computed only to derive the categorical concepts above — not standalone concepts): `shape.solidity` (→ mass-vs-NME), `shape.spiculation_index` (→ margin), `shape.n_components` (→ mass-vs-NME, multifocality). The BI-RADS cut-offs live in one table, `birads.THRESHOLDS`. Other IBSI shape features (`shape.elongation_flatness`, `shape.max_3d_diameter_mm`) are available in the library but are **not part of the derived concept set**.

`from breastmri_features import CITATIONS` gives the full reference list per feature.

## Validation against vendor values

On 384 I-SPY2 patients with vendor MRI features, the algorithmic derivations agreed with the
published sheet:

| Feature | Pearson r | alg/sheet ratio |
|---|---|---|
| Tumor volume | 0.998 | 0.95 |
| Sphericity | 0.815 | 1.10 |
| Longest diameter (RECIST) | 0.558 | ~1.0 (pre-calibrated) |
| BPE (%) | 0.298 | 1.31 |

Volume reproduces the vendor values almost exactly and sphericity closely; the longest diameter (a
manual RECIST measurement) is approximated by the in-plane caliper and is linearly calibrated to the
vendor scale (matched means, moderate correlation). BPE is the weakest — an in-crop proxy, a known
limitation, since true BPE requires whole-breast fibroglandular segmentation absent from the lesion crop.

## References

| # | Reference | Used for |
|---|-----------|----------|
| 1 | Zwanenburg A, et al. The Image Biomarker Standardization Initiative (IBSI). *Radiology* 2020;295:328–338. | sphericity, volume, diameter |
| 2 | van Griethuysen JJM, et al. Computational Radiomics System to Decode the Radiographic Phenotype (PyRadiomics). *Cancer Res* 2017;77:e104–e107. | reference-implementation parity |
| 3 | Eisenhauer EA, et al. New response evaluation criteria in solid tumours: RECIST 1.1. *Eur J Cancer* 2009;45:228–247. | longest diameter |
| 4 | van der Walt S, et al. scikit-image: image processing in Python. *PeerJ* 2014;2:e453. | 3-D solidity |
| 5 | Gilhuijs KGA, Giger ML, Bick U. Computerized analysis of breast lesions in 3-D using dynamic MRI. *Med Phys* 1998;25:1647–1654. | margin / spiculation |
| 6 | Huang YH, et al. CAD of mass-like lesions in breast MRI: 3-D morphology. *Comput Methods Programs Biomed* 2013;112:508–517. | margin, mass vs. NME, internal enhancement |
| 7 | Agliozzo S, et al. CAD for DCE breast MRI of mass-like lesions. *Med Phys* 2012;39:1704–1715. | kinetic curve, initial/internal enhancement |
| 8 | Dalmış MU, et al. CAD system for breast DCE-MRI at high spatiotemporal resolution. *Med Phys* 2016;43:6260–6273. | kinetic curve, initial enhancement |
| 9 | Wei D, et al. Fully automatic quantification of FGT and BPE. *Med Phys* 2021;48:238–252. | BPE |
| 10 | Zhu Y, et al. Automated calculation of BPE as a biomarker. *Breast Cancer Res Treat* 2026. | BPE categories |

The same list is available programmatically: `from breastmri_features import CITATIONS`.

## License

MIT.

## AI usage statement

This package was developed with the assistance of a large language model (Anthropic's Claude). All code was reviewed, tested, and validated
against vendor reference values by the authors, who take full responsibility for its correctness and
content. 
