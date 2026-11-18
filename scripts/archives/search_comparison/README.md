# Search comparison (archived 2026-11-17)

`03_04_search_stage.R` is the former stage `search` of
`scripts/03_04_cv_covariance.R`. It compared, for the shrinkage level of each
type × time × gene set, a refit grid of 10 values, the exact path on 200 values
(one spectral decomposition per fold), `optimize()` on the path score and 10
random values; for the two-parameter arm (panels ≤ 1,500 genes) a 5 × 8 grid
against 40 random pairs.

The path, which is a grid scored at the cost of one fit, chose the best value at
no loss, and random search did no better. With one or two hyperparameters the
pipeline therefore keeps grid search only. The result stays in
`output/grn/cv/search_comparison.csv` (read by chapter 3,
`#sec-cov-search-results`). The header of the R file says how to rerun it.
