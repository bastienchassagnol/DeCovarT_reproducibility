# Raw data (immutable)

Do not edit files in this folder. Record the download command that
created each object. Large RDS files are gitignored; track them with DVC
when they exceed GitHub limits.

## Suppinger gastruloids (GEO)

Paper: Suppinger et al., Cell Stem Cell 2023.
DOI: https://doi.org/10.1016/j.stem.2023.04.018
GEO series: GSE229513 (single-cell), GSE229386 (bulk).

Example (adjust names to the GEO supplementary file list):

```bash
# cd data/raw
# curl -L -C - -o GSE229513_gastruloidsobject.rds \
#   "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE229513&format=file&file=..."
# curl -L -C - -o GSE229386_AllHTSeqCountsWithGeneNames.txt.gz \
#   "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE229386&format=file&file=GSE229386_AllHTSeqCountsWithGeneNames.txt.gz"
```

Confirm filenames on the GEO landing pages before running `curl`.
