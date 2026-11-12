# MSigDB Gene Sets

This directory contains instructions for obtaining the Molecular Signatures
Database (MSigDB) gene sets required to reproduce the gene set enrichment
analyses.

## Required Files

The analyses use two human gene-set collections from MSigDB version 2023.2.Hs:

| File                           | Collection          |
|--------------------------------|---------------------|
| `h.all.v2023.2.Hs.symbols.gmt` | Hallmark gene sets  |
| `c2.cp.v2023.2.Hs.symbols.gmt` | Canonical gene sets |

## Download Instructions

1. Visit the official MSigDB website: https://www.gsea-msigdb.org/gsea/msigdb/
2. Access the archived human gene-set collections for version 2023.2.Hs.
3. Download both required collections in GMT format using gene symbols.
4. Place the downloaded files in this directory without changing their filenames.

The expected directory structure is:

```text
input/msigdb_v2023.2.Hs_GMTs/
├── README.md
├── h.all.v2023.2.Hs.symbols.gmt
└── c2.cp.v2023.2.Hs.symbols.gmt
```

## Licensing

MSigDB gene sets are distributed by their respective rights holders and are
subject to the licensing terms provided by MSigDB.

The gene-set files are not redistributed in this repository. Users must obtain
them directly from the official source and comply with the applicable terms.

MSigDB licensing information:
https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp

## Reproducibility

These gene sets are required for the gene set enrichment analyses performed by
`scripts/limmas.R`.

The exact MSigDB release and filenames are specified above to ensure that
analyses can be reproduced using the same reference collections.
