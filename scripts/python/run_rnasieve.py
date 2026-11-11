# /// script
# requires-python = ">=3.11,<3.12"
# dependencies = [
#     "rnasieve==0.1.4",
#     "numpy==1.24.4",
#     "scipy==1.10.1",
#     "pandas==2.0.3",
#     "cvxpy==1.6.0",
# ]
# ///
"""Run the official RNA-Sieve implementation on one reference and its bulks.

Usage: uv run --script scripts/python/run_rnasieve.py <in_dir> <out_dir>

<in_dir> holds `bulk.csv` (genes x bulk samples, first column gene ids),
`celltypes.csv` (column `celltype`, in output order) and one
`reference_<k>.csv` per cell type k = 1..J (genes x reference cells, same
gene order as the bulk). All values are raw counts of the panel genes.

The fit is `rnasieve.preprocessing.model_from_raw_counts()` followed by
`RNASieveModel.predict()`, jointly over all bulk columns (the package
default). Preprocessing keeps the package defaults for the cell trimming
(2% per type) and the sparsity filter (a gene is kept when it is detected
in more than 20% of the cells of at least one type). The library-size
normalisation is switched off: the input is a gene panel, so the column
sum is not the library of the cell, and a cell without any count on the
panel would give a division by zero. The pinned versions are those of the
DeCovarT competitor runners (`scripts/benchmark_metrics/python/uv.lock`).

Writes <out_dir>/{estimates,lower,upper}.csv (cell types x bulk columns)
and <out_dir>/genes_kept.csv. `lower`/`upper` are empty unless
`--intervals` is given; they are then the released
`compute_marginal_confidence_intervals()`, which rnasieve 0.1.4 sets to
[0, 1] when its averaged Hessian is ill-conditioned.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from rnasieve.preprocessing import model_from_raw_counts


def fit(in_dir: Path, level: float, intervals: bool) -> dict[str, pd.DataFrame]:
    """Fit RNA-Sieve on the files of `in_dir`; return the output tables."""
    bulk = pd.read_csv(in_dir / "bulk.csv", index_col=0)
    celltypes = pd.read_csv(in_dir / "celltypes.csv")["celltype"].tolist()
    raw_counts = {}
    for k, celltype in enumerate(celltypes, start=1):
        reference = pd.read_csv(in_dir / f"reference_{k}.csv", index_col=0)
        if list(reference.index) != list(bulk.index):
            raise ValueError(f"reference_{k}.csv genes differ from bulk.csv genes")
        raw_counts[celltype] = reference.to_numpy(dtype=float)

    model, psis = model_from_raw_counts(
        raw_counts,
        bulk.to_numpy(dtype=float),
        normalization=False,
    )
    alpha_hats = model.predict(psis)
    labels = list(alpha_hats.columns)

    def frame(values) -> pd.DataFrame:
        out = pd.DataFrame(np.asarray(values).T, index=labels, columns=bulk.columns)
        return out.loc[celltypes]

    tables = {"estimates": frame(alpha_hats.to_numpy())}
    missing = np.full(alpha_hats.shape, np.nan)
    tables["lower"] = frame(missing)
    tables["upper"] = frame(missing)
    if intervals:
        cis = model.compute_marginal_confidence_intervals(sig=1.0 - level)
        tables["lower"] = frame([[ci[0] for ci in row] for row in cis])
        tables["upper"] = frame([[ci[1] for ci in row] for row in cis])
    tables["genes_kept"] = pd.DataFrame({"n_genes": [psis.shape[0]]})
    return tables


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("in_dir", type=Path)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--level", type=float, default=0.95)
    parser.add_argument("--intervals", action="store_true")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for name, table in fit(args.in_dir, args.level, args.intervals).items():
        table.to_csv(args.out_dir / f"{name}.csv", index=name != "genes_kept")
    print("RNA-Sieve: done.", flush=True)


if __name__ == "__main__":
    main()
