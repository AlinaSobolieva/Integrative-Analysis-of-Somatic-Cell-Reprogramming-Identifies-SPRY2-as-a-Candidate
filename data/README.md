# Data

This directory contains the input data used for the transcriptomic analysis of somatic cell reprogramming toward pluripotency.

## Dataset

The analysis is based on publicly available RNA-seq data from the NCBI Gene Expression Omnibus (GEO), specifically the RNA-seq SubSeries **GSE90894**, which forms part of the **GSE90895 SuperSeries**.

The primary expression matrix used throughout the analysis is:
- `GSE90894_RPKM_mRNAseq_table.xlsx` — RPKM expression matrix used for data preprocessing, quality control, differential expression analysis, functional enrichment, and SPRY2-focused analyses.

The main reprogramming trajectory analyzed in this study comprises:
**MEFs → 48 h → pre-iPSC → ESCs**

## Data organization

The `data/` directory separates original input data from datasets generated during analysis:

```text
data/
├── raw/
│   └── GSE90894/
│       └── GSE90894_RPKM_mRNAseq_table.xlsx
│
└── processed/
```
