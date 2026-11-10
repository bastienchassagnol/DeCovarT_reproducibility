# Raw data (immutable)

Do not edit files in this folder. Record the download command that
created each object. Large RDS files are gitignored; track them with DVC
when they exceed GitHub limits.

## Suppinger gastruloids (GEO)

Paper: Suppinger et al., Cell Stem Cell 2023.
DOI: <https://doi.org/10.1016/j.stem.2023.04.018>
GEO series: GSE229513 (single-cell), GSE229386 (bulk).

Both files are tracked with DVC in this repository. The pointers come
from the sibling repository `GastroDeconv2FateMap` and share its remote
(`mmg_cluster`). From **this** repository root:

```bash
dvc pull data/raw/GSE229513_gastruloidsobject.rds.dvc
dvc pull data/raw/GSE229386_AllHTSeqCountsWithGeneNames.txt.gz.dvc
```

| File | Content | Size |
|---|---|---:|
| `GSE229513_gastruloidsobject.rds` | Authors' Seurat 4.1.0 object: RNA counts (27,027 genes × 71,005 cells, Salmon/Alevin), SCT and integrated assays (2,944 integration genes) | about 14.6 GB |
| `GSE229386_AllHTSeqCountsWithGeneNames.txt.gz` | Bulk RNA-seq counts, one row per Ensembl id, column `gene_name`, early-Wnt (`_Early_`) and standard-Wnt (`_CTL_`) arms at 48, 72 and 96 h | about 1.1 MB |

Pull the RDS on a machine whose `data/raw` is not cloud-synced. The
HTSeq table is the bulk count matrix, not the single-cell counts.
