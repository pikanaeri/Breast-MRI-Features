"""Method references behind each derived feature. `CITATIONS[feature]` -> list of references."""

CITATIONS = {
    "sphericity": ["Zwanenburg A, et al. The Image Biomarker Standardization Initiative. "
                   "Radiology 2020;295:328-338.",
                   "van Griethuysen JJM, et al. Computational Radiomics System to Decode the "
                   "Radiographic Phenotype. Cancer Res 2017;77:e104-e107."],
    "tumor_volume": ["Zwanenburg A, et al. IBSI. Radiology 2020 (Volume, voxel counting)."],
    "longest_diameter": ["Zwanenburg A, et al. IBSI. Radiology 2020 (Maximum 2-D diameter, slice).",
                         "RECIST 1.1: Eisenhauer EA, et al. Eur J Cancer 2009;45:228-247."],
    "solidity": ["van der Walt S, et al. scikit-image. PeerJ 2014;2:e453 (regionprops solidity)."],
    "shape": ["IBSI sphericity thresholded to BI-RADS round/oval vs irregular (investigator rule)."],
    "margin": ["Gilhuijs KGA, Giger ML, Bick U. Computerized analysis of breast lesions in 3-D "
               "using dynamic MRI. Med Phys 1998;25:1647-1654.",
               "Huang YH, et al. Differential analysis of 3-D morphology. "
               "Comput Methods Programs Biomed 2013;112:508-517."],
    "mass_vs_nme": ["Solidity / foci rule following 3-D morphological CAD descriptors "
                    "(Agliozzo 2012; Huang 2013)."],
    "kinetic_curve": ["Agliozzo S, et al. CAD for DCE breast MRI. Med Phys 2012;39:1704-1715.",
                      "Dalmis MU, et al. CAD system for breast DCE-MRI. Med Phys 2016;43:6260-6273."],
    "initial_enhancement_rate": ["Agliozzo 2012; Dalmis 2016 (early-phase enhancement)."],
    "internal_enhancement": ["Rim / heterogeneity descriptors (Agliozzo 2012; Huang 2013)."],
    "bpe": ["Wei D, et al. Fully automatic FGT and BPE quantification. Med Phys 2021;48:238-252.",
            "Zhu Y, et al. Automated BPE as a biomarker. Breast Cancer Res Treat 2026."],
}
