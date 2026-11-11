# /// script
# requires-python = ">=3.11,<3.12"
# dependencies = [
#     "redeconv==1.2.1",
#     "pandas==2.0.3",
# ]
# ///
"""Run the official ReDeconv implementation on one signature and its bulks.

Usage: uv run --script scripts/python/run_redeconv.py <in_dir> <out_dir>

<in_dir> holds two tab-separated files written by
`deconvolute_redeconv()` (R/04_01_stats_deconvolution.R):

- `signature.tsv`: first column the gene, then one `<type>_mean` and one
  `<type>_std` column per cell type, the format read by
  `redeconv.__ReDeconv_P.read_sig_gene_mean_std()`;
- `bulk.tsv`: first column the gene, then one column per bulk sample, the
  format read by `read_bulk_RNA_seq_column_for_sample()`.

The default entry point `ReDeconv()` is called (MLMD, six reweighting
loops). Writing the signature from supplied moments bypasses ReDeconv's
own signature-gene selection (t-tests, fold change, top genes by
coefficient of variation): every gene of the file is used. ReDeconv
reports no confidence interval.

Writes <out_dir>/estimates.csv (cell types x bulk columns) and keeps
ReDeconv's own table in <out_dir>/percentages.tsv. Adapted from the
DeCovarT competitor runner
(`scripts/benchmark_metrics/python/run_redeconv.py`), same package
version.
"""

import argparse
import contextlib
import importlib
import io
import os
from pathlib import Path

import pandas as pd

os.environ.setdefault("MPLBACKEND", "Agg")


def fit(in_dir: Path, out_dir: Path) -> pd.DataFrame:
    """Fit ReDeconv on the files of `in_dir`; return types x samples."""
    redeconv_p = importlib.import_module("redeconv.__ReDeconv_P")
    result_file = out_dir / "percentages.tsv"
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        redeconv_p.ReDeconv(
            str(in_dir / "signature.tsv"),
            str(in_dir / "bulk.tsv"),
            str(result_file),
        )
    result = pd.read_csv(result_file, sep="\t", index_col=0)
    samples = pd.read_csv(in_dir / "bulk.tsv", sep="\t", index_col=0, nrows=0)
    return result.T.loc[:, list(samples.columns)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("in_dir", type=Path)
    parser.add_argument("out_dir", type=Path)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    fit(args.in_dir, args.out_dir).to_csv(args.out_dir / "estimates.csv")


if __name__ == "__main__":
    main()
