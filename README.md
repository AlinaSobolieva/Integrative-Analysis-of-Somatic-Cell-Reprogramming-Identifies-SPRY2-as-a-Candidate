# SPRY2-reprogramming-transcriptomic-analysis

Reproducible R-based integrative transcriptomic and network analysis of somatic cell reprogramming toward pluripotency (*Mus musculus*).

## Overview
This repository contains the analysis code, data processing pipeline, and results for investigating stage-specific transcriptional remodeling during somatic cell reprogramming (MEF → 48h OSKM → pre-iPSC → ESC) and characterizing *SPRY2* as an exploratory candidate gene.

## Repository Structure
```
SPRY2-reprogramming-transcriptomic-analysis/
├── README.md
├── LICENSE
├── data/
│   └── README.md
├── scripts/
│   ├── 01_import_cleaning.R
│   ├── 02_quality_control.R
│   ├── 03_differential_expression.R
│   ├── 04_functional_enrichment.R
│   ├── 05_SPRY2_analysis.R
│   └── 06_network_analysis.R
├── results/
│   ├── figures/
│   │   ├── main/
│   │   │   ├── Figure_1_QC.png
│   │   │   ├── Figure_2_DE.png
│   │   │   ├── Figure_3_Enrichment.png
│   │   │   └── Figure_4_SPRY2.png
│   │   └── supplementary/
│   │       ├── Supplementary_Figure_S1.png
│   │       └── Supplementary_Figure_S2.png
│   └── tables/
│       ├── Supplementary_Table_S1.xlsx
│       ├── Supplementary_Table_S2.xlsx
│       └── Supplementary_Table_S3.xlsx
└── manuscript/
    └── README.md
```

## Authors
- **Alina Sobolieva** ([ORCID: 0009-0006-6536-7155](https://orcid.org/0009-0006-6536-7155))
- **Khrystyna Malysheva** ([ORCID: 0000-0001-8396-6799](https://orcid.org/0000-0001-8396-6799))
- Stepan Gzhytskyi National University of Veterinary Medicine and Biotechnologies Lviv

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
