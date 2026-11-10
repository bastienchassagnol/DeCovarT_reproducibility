
● Both sections are rewritten, but I have not rendered the book, so the new
  table chunk and cross-references are untested in Quarto (the table's R code
  runs on its own).

  Your message stopped after scenario 1, so I took scenario 2 to be the runs of
  this book, where every method gets the same gene panel. Say if you meant
  something else.

  Chapter 3: reading the shrinkage estimate

  The section is now Regularisation shrinks covariance terms without setting
  them to zero (docs/03-covariance-estimation.qmd#L849). It covers:

  - Status: the first-pass covariances use the closed-form λ, not the
    cross-validated one the chapter announces.
  - What the estimator returns: variances untouched, every correlation
    multiplied by 1 − λ, none set to zero, and λ as a floor on the eigenvalues.
  - A computed table of λ per type and time, read from
    output/grn/plugin/plugin_summary.csv. Abundant types get λ of 0.006 to 0.03;
    types with 29 to 66 cells get 0.60 to 0.74.
  - Why a dense estimate does not mean every pair interacts: your three bullets,
    reworded for the book and linked to the sparse arms.
  - The denoised arm: for the rarest types the matrix has rank at most cells − 1
    before shrinkage.
  - A TODO pointing to the regularisation sweep in chapter 6.

  That folder (output/grn/plugin/) is still untracked, so the table will not
  render from a fresh clone until its CSVs are committed.

  The eigenvalue statement is written inline and is not in the verification
  ledger.

  Chapter 4: two scenarios

  Statistical interpretation of the results
  (docs/04-deconvolution-pipeline.qmd#L97) now has two subsections.

  - Scenario 1, each method selects its own genes.
    - It quotes the omnideconv policy you cited (p. 5).
    - A table sets DWLS, SCDC and BayesPrism side by side: genes kept, gene
      weights, input scale. I read the selection rules from the installed
      packages, and the input scales from the paper.
    - Your figure stays, with a caption saying it is AI-generated and a
      hypothesis, not a result.
  - Scenario 2, one panel for every method. It describes the current runs and
    what a difference means there. Results are left as a TODO until the 96 h
    fits finish.

  The paper does not support the hypothesis directly, and the text says so.
  Dietrich et al. report DWLS among the most stable methods across annotation
  resolutions (p. 8) and the one with the least spillover (p. 13). Their one
  signature comparison finds the better result with the lower condition number
  (p. 16). Their references are immune and tumour atlases, not a developmental
  continuum.

  Extraction

  I extracted the PDF with local tools only (pdftotext, pdftoppm); nothing was
  uploaded. The folder is
  .firecrawl/extract/dietrichOmnideconvUnifyingFramework2026/, with the text,
  six figure-page renders, figure records, a content summary and metadata. I
  added .firecrawl/ to .gitignore.

  - Skipped: supplementary files, runtime results, pages 17 to 25 and the
    reference list.
  - Citation: the paper is not in Zotero or the book's .bib, so the chapter
    carries a TODO-cite instead of a key. The folder name is a provisional slug,
    not a verified key.

  Runs

  At 11:27 the 96 h Welch-gated DeCovarT fits were at 14 of 24 chains. The
  second wave is slower than the first, so the end is more likely 12:15 to
  12:45. Nothing is committed.

❯ are the two backgroudn agents: Wait for the large-set comparators to finish,
  ... clsoe to be finished? where are their logs stored?
  ⎿  Referenced file docs/03-covariance-estimation.qmd
  ⎿  Referenced file docs/06-sensitivity-misspecification.qmd
  ⎿  Read convergence_failures.md (192 lines)
  ⎿  Referenced file R/04_01_vis_deconvolution.R
  ⎿  Referenced file R/04_01_stats_deconvolution.R

● Task "Wait for the large-set comparators to finish" still running in
  background

● Task "Wait for the large-set DeCovarT fits at 72 h and 96 h" still running in
  background
  ⎿  Plan file referenced
     (~/.claude/plans/the-general-purpose-is-parsed-pond.md)

● Compacted using a summary written while you were away · turn off with
  "idleCompaction": false in settings.json

  Read 3 files, ran 5 shell commands (ctrl+o to expand)
  ⎿  Loaded .claude/rules/data-visualisation.md

● Both runs are finished; the two "waiters" are not. They are stuck in their
  loop because of a mistake in my wait command: pgrep -f matches the waiter's
  own command line, so they would never have notified. They cost nothing, and I
  left them running since you asked me not to stop background tasks; say the
  word and I stop them.

  Where the logs are

  All in logs/:

  ┌─────────────┬────────────────────────────────────────────────┬─────────┐
  │     Run     │                      Log                       │ Finishe │
  │             │                                                │    d    │
  ├─────────────┼────────────────────────────────────────────────┼─────────┤
  │ Comparators │                                                │         │
  │ ,           │ logs/04_01_comparators_welch_2026-11-11.log    │ 09:05   │
  │ Welch-gated │                                                │         │
  ├─────────────┼────────────────────────────────────────────────┼─────────┤
  │ DeCovarT, W │                                                │         │
  │ elch-gated, │ logs/04_01_decovart_welch_72h_96h_2026-11-11.l │ 11:59   │
  │  72 h and   │ og                                             │         │
  │ 96 h        │                                                │         │
  ├─────────────┼────────────────────────────────────────────────┼─────────┤
  │ DeCovarT, W │                                                │         │
  │ elch-gated, │ logs/04_01_decovart_welch_48h_2026-11-11.log   │ 00:25   │
  │  48 h       │                                                │         │
  ├─────────────┼────────────────────────────────────────────────┼─────────┤
  │ Small gene  │ logs/04_01_r2_05_2026-11-11.log                │ 00:36   │
  │ sets        │                                                │         │
  └─────────────┴────────────────────────────────────────────────┴─────────┘

  Per-chain checkpoints are in output/deconvolution/fits/ (one RDS per chain).
  The waiters' own output is under /tmp/claude-1007/…/tasks/ and is empty.

  The 96 h run completed all 24 chains with no error, at 68 to 97 min per chain.
  Every comparator reports ok in all 18 cells of the design.

  What I did next

  Since the fits were done, I rebuilt the tables (2,904 rows, 44 methods) and
  redrew the seven figures in docs/images/deconvolution/. Both scripts exited
  cleanly. I looked at three of the seven figures (Welch-gated scatter, RMSE dot
  plot, 96 h estimates); the other four are unchecked.

  First reading of the full results

  - Comparators win on pooled error. On the Welch-gated genes, DeMixSC wNNLS has
    RMSE 0.054 (r = 0.995) and RNA-Sieve 0.078 (r = 0.98), pooled over the
    three times. The best DeCovarT flavour is at 0.24 (raw covariance, L-BFGS-B,
    bulk divided by b_g).
  - 96 h is where DeCovarT fails. On the Welch-gated genes, the three denoised
    inputs put nearly all the mass on caudal mesoderm (observed 0.22) and almost
    none on pre-somitic mesoderm (observed 0.73), giving RMSE of about 0.47. On
    the small gene set the same inputs give about 0.2.
  - 48 h, Σ_bio alone, Welch-gated: the estimate goes to the wrong vertex (r =
    −0.49, RMSE 0.76). Adding the Poisson term repairs it (RMSE 0.02).
  - Newton failures over the full design: 162 of 288 fits, still 0 of 288 for
    L-BFGS-B.

  ┌──────┬───────────────┬───────┬─────────────────┐
  │ Time │   Gene set    │ Genes │ Newton failures │
  ├──────┼───────────────┼───────┼─────────────────┤
  │ 48 h │ deviance ≥ 5% │    17 │          0 / 48 │
  ├──────┼───────────────┼───────┼─────────────────┤
  │ 72 h │ deviance ≥ 5% │   140 │         22 / 48 │
  ├──────┼───────────────┼───────┼─────────────────┤
  │ 96 h │ deviance ≥ 5% │   635 │         36 / 48 │
  ├──────┼───────────────┼───────┼─────────────────┤
  │ 48 h │ Welch-gated   │ 1,417 │         42 / 48 │
  ├──────┼───────────────┼───────┼─────────────────┤
  │ 72 h │ Welch-gated   │ 3,420 │         17 / 48 │
  ├──────┼───────────────┼───────┼─────────────────┤
  │ 96 h │ Welch-gated   │ 5,821 │         45 / 48 │
  └──────┴───────────────┴───────┴─────────────────┘

  The 72 h Welch-gated row breaks the "more genes, more failures" trend I wrote
  in convergence_failures.md (section 2.1). I have not updated that note.

  One figure problem: at 96 h the five cell-type colours of the repository
  palette are all light greens and hard to tell apart in the scatter.

  Not done: the scenario 2 TODO in chapter 4 is still unfilled, the
  moments_*.rds files are not under DVC, and nothing is committed.

❯ if nothing is longer ran, stop the backgorund taks.

  then try to identify several theorical and biolgoical epxlantions epxlaining
  why DeMixSC wNNLS and RNA-sieve thend to outperform decovart, espcially at 96
  h. why most of the mass foes to caudal mesoderm and almost none on
  pre-someitic mesoderm?

  are these two cell types too storngly corrleated (but accoutnogn for the
  covariances shoudl have helped to set them aprt)? also, why so much
  discrpenacy with rna-sieve, since our theoreicla framewokr is highly simialr
  to it, excpet that we try to incprorate for the interactions ebtwene the
  genes (might be intepreredt as a multivairate marginal generalisation fot he
  model)? With mutlviariate models, on average, we're expecting the same poitn
  wie estimates (same marginals indeed), but only differences in the derivation
  of confidnce itnervals (univeriate aprpoahces tendign to be
  overconfident)Maybe also our approahc to enforce the unit simplex
  paramrmaistrisaiton, using ILR repremaisatrion while rna-sieve

  Once done, rerun the piepline by:

  1) adding to it the mead and redeconv deconvolution algoirthm:
  https://github.com/bastienchassagnol/DeCovarT/blob/e71b5d714d4725ac052c6d3128
  1eaa39580b6f12/scripts/fig05_single_cell_hierarchical.R (I've archived
  previosu deocnvolution outputs: /mnt/DATA_11TB/projects/dtoo_project/DeCovarT
  _reproducibility/output/archive_deconvolution, such that you do not overwtire
  results with the next pipeline 2) also use the compsotional aware aithcison
  distance to compelment @docs/images/deconvolution/benchmark_pearson.png  and
  @docs/images/deconvolution/benchmark_rmse.png  .
  3) write locally a rna-sieve multviariate inspired version (so the same
  algorithm, usign an alternated or block graident approach, marginal variance
  inferred during the optimisation framework, ineuqality cosntraint only for
  the cellualr ratios that are post hoc nromalsied to reahc the unti simplex
  using quadratic prgramming, and the averaged expression and cell-type
  specific covariances refined in the second stage of the deocnvolution
  algorithm, in a process seemingly simarl to declas deocnovlution: ) t in R
  instead of Python, dropping the assumption of indpendence between genes
  belgongin to the same cell type), see @decals.pdf  and @rna-sieve-paper.pdf .
  in breif, I have two hypothses to explain the discrpnayc: either our
  repramisatrion strategy is not that eprformant comapred agaisnt a quadratic
  programming even thogh theorically they shoudl covnerge to the same final mle
  esitmate, but nuemricla dicrpancies can exist, with ILR posisbly tending to
  project in toehr dimensions, or the covariance and averaged profile is
  strongly impact and different in bulk profiles, something that we'r ento
  cpaturing by using plug-in estimates inferred bhegorahd  int he single-cell
  ⎿  Referenced PDF decals.pdf (30 pages)
  ⎿  Referenced PDF rna-sieve-paper.pdf (36 pages)

● 3 background agents launched (↓ to manage)
   ├ MEAD/ReDeconv and DeCovarT optimiser
   ├ RNA-Sieve and DeCLAS algorithms
   └ Map current benchmark pipeline

  Read 1 file, ran 7 shell commands (ctrl+o to expand)

● Agent "MEAD/ReDeconv and DeCovarT optimiser" finished · 2m 48s

● Agent "Map current benchmark pipeline" finished · 2m 52s

● Task "Wait for the large-set DeCovarT fits at 72 h and 96 h" was stopped by
the user or by the app connected to this session

● Task "Wait for the large-set comparators to finish" was stopped by the user or
by the app connected to this session

  Read 2 files, ran 2 shell commands (ctrl+o to expand)

● Updated plan
  ⎿  /plan to preview

  Ran 1 shell command (ctrl+o to expand)

● Agent "RNA-Sieve and DeCLAS algorithms" finished · 7m 10s

  Ran 1 shell command (ctrl+o to expand)

● User answered Claude's questions:
  ⎿  · MEAD needs at least two individuals in the reference (it weights genes by
     between-individual variance). The scored cells come from two cell lines at
     48 h but one line only at 72 h and 96 h. How should MEAD get its
     individuals? → Pseudo-individuals (Recommended)
     · DECALS estimates each cell-type covariance by regressing residual
     cross-products on squared proportions across 541 bulks; here there are 6
     bulks per time and up to 5 types, so a full gene-by-gene covariance cannot
     be re-estimated from the bulk. What should the second stage of the
     multivariate RNA-Sieve refine? → Means + gene-wise variances (Recommended)
     · The Aitchison distance is infinite when an estimated proportion is zero,
     and many estimates sit on a face (2e-14). The smallest observed fraction is
     0.0033. Which zero replacement should the metric use? → Floor at 0.001
     (Recommended)

● Updated plan
  ⎿  /plan to preview

● The two stale waiters are closed and no model run is active. The optimiser
  isn't the cause, since the log-determinant term drives DeCovarT toward a
  vertex; a QP on the same plug-ins recovers the composition on that sample
  (caudal 0.22, PSM 0.78), and RNA-Sieve's proportion update silently fails
  under the pinned cvxpy. I'm reusing the 180 archived checkpoints by copy next,
  unless you want a full rerun instead, which would take about 30 CPU-hours.
  · summary
  ⎿  Error: The user doesn't want to proceed with this tool use. The tool use
     was rejected (eg. if it was a file edit, the new_string was NOT written to
     the file). The user provided the following reason for the rejection:  close
     to eprfect plannign. Report findings in document
     04-deconvolution-pieplien.qmd; including your first prelimianry asusmptions

     Comments on the plan:
     [Re: "H2 is supported for the means. Unweighted NNLS gives 100% caudal
     mesoderm. Mean profiles are near-collinear (PSM against caudal r = 0.98,
     PSM against somite front r = 0.99, small panel). Covariance weighting
     separates them on the large panel (Mahalanobis separation 3.5 times the
     diagonal one) but not on the small one (0.22 times)."] Add using geom_tile
     the Pearson correlation ebtween cell types (and also kenda) for the
     averaged single cell profiles, splitte dby the 3 time points.
     [Re: "RNA-Sieve as run is a weighted least squares. In rnasieve 0.1.4 the
     proportion QP passes its constraints to cp.Minimize; under the pinned cvxpy
     1.6.0 that raises a TypeError (checked), which the code swallows.
     Proportions then come only from the final weighted least-squares projection
     (clip, renormalise). Consequence read from the code, not traced at run
     time."] is DeCovarT adapatable or currently present int eh form of a robsut
     mutviariate wieghted lienar regresison model?
     [Re: "The difference is the log-determinant. Between those two compositions
     the quadratic term moves by 75 units and −½ log det S(p) by 1,075. Each
     Σ_bio,j is close to rank C_j − 1; mixing fills the null directions, so log
     det rises steeply away from a vertex. With a diagonal covariance the same
     term moves by about 20."] does htis invovle to use an even mroe regularised
     version fot ehc voariance amtrix to further pnealise non diaogonal terms?
     [Re: "y_b ~ N(M w_b, D^½ V(w_b, M) D^½), with V = Σ_k w_k (Σ_k + m_k m_kᵀ)
     − (M w)(M w)ᵀ / Σ_k w_k. Mean and variance share n, as in RNA-Sieve. The
     between-type term is the multivariate form of the paper's eq. 7; the paper
     does not state it (my derivation, flagged in the code, with a switch)."] so
     we're adding as well covariance between cell types (at the proprotion
     level)? in adidtion to the cell type specific variances?
  Ran 1 shell command (ctrl+o to expand)

● Updated plan
  ⎿  /plan to preview

● Updated plan
  ⎿  /plan to preview

● Updated plan
  ⎿  /plan to preview

● Updated plan
  ⎿  /plan to preview

● Updated plan
  ⎿  /plan to preview

● Updated plan
  ⎿  /plan to preview
  ⎿  Error: User chose to stay in plan mode and continue planning

     Comments on the plan:
     [Re: "Does V add a covariance between cell types? No. Cell types stay
     independent. The extra term is the law of total covariance for one cell
     drawn at random: Cov(x) = Σ_k p_k Σ_k (within types) + Σ_k p_k m_k m_kᵀ −
     (Mp)(Mp)ᵀ (which type the cell is). The second part has rank at most K − 1.
     DeCovarT's S_W(p) keeps only the first part, i.e. it treats the number of
     cells of each type as fixed. The switch is stored both ways for the qp
     stage so its effect is visible."] indeed, the covariance must account for
     the fact that the total number of cells in a given cell type is lelyto
     vary, imapcting the covariance. Add a version or flavour fo decovartwitht
     aht update.
     [Re: "Step 1 · MEAD and ReDeconv
     Files: R/04_01_stats_deconvolution.R (two functions),
     scripts/04_01_stats_deconvolution.R (runners l. 299, comparator_label l.
     392), scripts/python/run_redeconv.py (new), R/04_01_vis_deconvolution.R l.
     39 (comparator levels).

     MEAD 1.0.2 (DongyueXie/MEAD @ d66fa10, mengyin/vashr, Rfast;
     renv::install() then renv::snapshot()). Call copied from
     ../DeCovarT/scripts/benchmark_metrics/second_generation_methods.R l.
     356–424: get_R01(bulk, alpha = 0.1) (fallback NULL), then MEAD(bulk, ref,
     cell_types, preprocessing_control = list(filter.gene = FALSE), R01 = r01).
     Pseudo-individuals (your choice): cells of each type split at random into 4
     groups, fixed seed. Label: "MEAD (pseudo-individuals)".
     MEAD divides each cell by its sum over the panel, so p_hat is a transcript
     fraction: convert with p_j ∝ q_j / s_j, s_j = colSums(mu)[j]. Cells with no
     count on the panel are dropped and counted in the status message.
     … +2 lines (ctrl+o to expand)

● Updated plan
  ⎿  /plan to preview

● Updated plan
  ⎿  /plan to preview

● Updated plan
  ⎿  /plan to preview

● Updated plan
  ⎿  /plan to preview

● Updated plan
  ⎿  /plan to preview
  ⎿  Error: User chose to stay in plan mode and continue planning

     Comments on the plan:
     [Re: "Step 3 · Multivariate RNA-Sieve in R
     New R/04_02_stats_sieve_multivariate.R; driver stage sieve; quadprog is
     already in the renv library. Code comments use the symbols of
     vignettes/_notation.qmd.

     Model (per bulk b; unnormalised weights w = n p ≥ 0, p = w / Σw):

     y_b ~ N(M w_b, D^½ V(w_b, M) D^½), with V = Σ_k w_k (Σ_k + m_k m_kᵀ) − (M
     w)(M w)ᵀ / Σ_k w_k. Mean and variance share n, as in RNA-Sieve. The
     between-type term is the multivariate form of the paper's eq. 7; the paper
     does not state it (my derivation, flagged in the code, with a switch).
     Reference: μ̃_k ~ N(m_k, Σ_k / C_k). D = diag(d_g), start at 1.
     … +15 lines (ctrl+o to expand)

  Ran 2 shell commands (ctrl+o to expand)

● Updated plan
  ⎿  /plan to preview

● Updated plan
  ⎿  /plan to preview

● Updated plan
  ⎿  /plan to preview

● User approved Claude's plan
  ⎿  Plan saved to: ~/.claude/plans/the-general-purpose-is-parsed-pond.md ·
     /plan to edit
     Second pass of the deconvolution benchmark: why DeCovarT loses at 96 h, and
     three additions

     Context

     The first pass (now in output/archive_deconvolution/, 180 checkpoints)
     put DeMixSC wNNLS (pooled RMSE 0.054) and RNA-Sieve (0.078) ahead of
     every DeCovarT flavour (best 0.24). At 96 h on the Welch-gated panel the
     denoised DeCovarT fits give about 0.98 caudal mesoderm (observed 0.22)
     and 0.01 pre-somitic mesoderm, PSM (observed 0.73).

     You asked for explanations, then a rerun with MEAD and ReDeconv, the
     Aitchison distance, and a multivariate RNA-Sieve written in R, to
     separate two hypotheses: H1 ILR optimiser against quadratic
     programming (QP); H2 bulk means and covariances differ from the
     single-cell plug-ins.

     Findings so far (read-only checks, scripts in the session scratchpad)

     One sample, 96h_CTL_rep1, Welch-gated panel (5,821 genes):

     ┌──────────────────────────────┬────────┬───────┬───────┐
     │          Estimator           │ Caudal │  PSM  │  Gut  │
     ├──────────────────────────────┼────────┼───────┼───────┤
     │ Observed                     │  0.215 │ 0.727 │ 0.028 │
     ├──────────────────────────────┼────────┼───────┼───────┤
     │ Unweighted NNLS              │  1.000 │ 0.000 │ 0.000 │
     ├──────────────────────────────┼────────┼───────┼───────┤
     │ GLS by QP, diagonal of raw Σ │  0.408 │ 0.590 │ 0.002 │
     ├──────────────────────────────┼────────┼───────┼───────┤
     │ GLS by QP, full raw Σ        │  0.219 │ 0.781 │ 0.000 │
     ├──────────────────────────────┼────────┼───────┼───────┤
     │ GLS by QP, full Σ_bio        │  0.000 │ 0.740 │ 0.012 │
     ├──────────────────────────────┼────────┼───────┼───────┤
     │ DeCovarT L-BFGS-B, raw Σ     │  0.066 │ 0.224 │ 0.700 │
     ├──────────────────────────────┼────────┼───────┼───────┤
     │ DeCovarT L-BFGS-B, Σ_bio     │  0.983 │ 0.009 │ 0.002 │
     └──────────────────────────────┴────────┴───────┴───────┘

     1. H1 is not supported as an optimiser failure. DeCovarT's estimate
        has a higher profile log-likelihood than the observed composition
        (Σ_bio: −57,115 against −58,264) and than the QP estimate. L-BFGS-B
        does not use ILR either (it works on non-negative weights). The
        objective differs, not the solver.
     2. The difference is the log-determinant. Between those two
        compositions the quadratic term moves by 75 units and
        −½ log det S(p) by 1,075. Each Σ_bio,j is close to rank C_j − 1;
        mixing fills the null directions, so log det rises steeply away from
        a vertex. With a diagonal covariance the same term moves by about 20.
     3. Profiling φ lets that term win. The fitted variance scale is
        10⁷ to 10⁸ (n̂ ≈ 0.05 cells): the residual is platform mismatch,
        not cell sampling. Profiling turns the quadratic term into
        −(G/2) log q(p), which is flat. RNA-Sieve ties mean and variance to
        one n, so its quadratic term dominates.
     4. RNA-Sieve as run is a weighted least squares. In rnasieve
        0.1.4 the proportion QP passes its constraints to cp.Minimize;
        under the pinned cvxpy 1.6.0 that raises a TypeError (checked),
        which the code swallows. Proportions then come only from the final
        weighted least-squares projection (clip, renormalise). Consequence
        read from the code, not traced at run time.
     5. H2 is supported for the means. Unweighted NNLS gives 100% caudal
        mesoderm. Mean profiles are near-collinear (PSM against caudal
        r = 0.98, PSM against somite front r = 0.99, small panel). Covariance
        weighting separates them on the large panel (Mahalanobis separation
        3.5 times the diagonal one) but not on the small one (0.22 times).
     6. On the small panel (635 genes) the ranking reverses: GLS-QP is poor
        (caudal 0.89), DeCovarT better (PSM about 0.5). The rerun covers all
        18 samples.

     Answers to your three comments

     - Is DeCovarT a robust multivariate weighted regression, or can it
       become one? Not today. deconvolute_single_cell() always maximises
       a likelihood that contains log det S(p), with profile = FALSE too;
       the package has no QP and no robust weights (deconvolute_ratios_rlm
       is robust but ignores covariances). It is adaptable: freezing S(p)
       and solving the weighted least squares is exactly the qp stage of
       step 3 with the full Σ, i.e. a multivariate GLS. Huber weights on
       the whitened residuals would make it robust; that is one more line
       in the same loop, added as stage qp_huber (step 3).
     - Does the log-determinant finding call for more regularisation?
       Partly. The cause is the near-zero eigenvalues, not the off-diagonal
       terms as such. Three remedies exist: more shrinkage (λ → 1 removes
       the effect: 20 units against 1,075), an eigenvalue floor, or dropping
       the log-determinant (GLS). Shrinkage also throws away the weighting
       that separated caudal from PSM on the large panel, so the answer
       needs the sweep of step 3b.
     - Does V add a covariance between cell types? No. Cell types stay
       independent. The extra term is the law of total covariance for one
       cell drawn at random: Cov(x) = Σ_k p_k Σ_k (within types) +
       Σ_k p_k m_k m_kᵀ − (Mp)(Mp)ᵀ (which type the cell is). The second
       part has rank at most K − 1. DeCovarT's S_W(p) keeps only the first
       part, i.e. it treats the number of cells of each type as fixed.
       The switch is stored both ways for the qp stage so its effect is
       visible, and a DeCovarT flavour with the full term is added
       (step 3d).

     Step 0 · Working folder

     Recreate output/deconvolution/ and copy the 180 checkpoints of
     output/archive_deconvolution/fits/ into it (archive untouched). The
     DeCovarT chains and the six existing comparators are deterministic and
     their code is unchanged; refitting would cost about 30 CPU-hours for
     the same numbers. New methods add checkpoints beside them.

     Step 1 · MEAD and ReDeconv

     Files: R/04_01_stats_deconvolution.R (two functions),
     scripts/04_01_stats_deconvolution.R (runners l. 299,
     comparator_label l. 392), scripts/python/run_redeconv.py (new),
     R/04_01_vis_deconvolution.R l. 39 (comparator levels).

     - MEAD 1.0.2 (DongyueXie/MEAD @ d66fa10, mengyin/vashr, Rfast;
       renv::install() then renv::snapshot()). Call copied from
       ../DeCovarT/scripts/benchmark_metrics/second_generation_methods.R
       l. 356–424: get_R01(bulk, alpha = 0.1) (fallback NULL), then
       MEAD(bulk, ref, cell_types, preprocessing_control = list(filter.gene =
       FALSE), R01 = r01).
       - Pseudo-individuals (your choice): cells of each type split at
         random into 4 groups, fixed seed. Label: "MEAD (pseudo-individuals)".
       - MEAD divides each cell by its sum over the panel, so p_hat is a
         transcript fraction: convert with p_j ∝ q_j / s_j,
         s_j = colSums(mu)[j]. Cells with no count on the panel are dropped
         and counted in the status message.
     - ReDeconv (redeconv==1.2.1, PEP 723 script as
       run_rnasieve.py): signature of mu and sqrt(diag(sigma_raw)),
       then redeconv.__ReDeconv_P.ReDeconv(), as in
       ../DeCovarT/scripts/benchmark_metrics/python/run_redeconv.py. Its own
       gene selection and size normalisation are bypassed (shared panel).
     - Robust regression: DeCovarT::deconvolute_ratios_rlm(y, mu)
       (MASS::rlm, M-estimator, no covariance), per sample, label "RLM".
       The package function stops when a coefficient is negative; the
       wrapper then clips at zero, renormalises, and says so in
       comparator_status.csv. If the function is not exported, the same
       MASS::rlm call is written in the wrapper (no :::).

     Step 2 · Aitchison distance

     - benchmark_metrics() (R/04_01_stats_deconvolution.R l. 704): add
       aitchison = Euclidean distance between clr vectors, after
       multiplicative zero replacement at 0.001 (your choice), and
       aitchison_01 at 0.01 as a sensitivity column.
     - Pooled table (driver l. 544–570): mean of the per-sample distances.
     - plot_benchmark_metric() (R/04_01_vis_deconvolution.R l. 228):
       scale and axis label from a named list per metric; add "aitchison".
       New benchmark_aitchison.{pdf,png}. Caption: floor value, and that
       distances compare within a time only (3 parts at 48/72 h, 5 at 96 h).

     Step 3 · Multivariate RNA-Sieve in R

     New R/04_02_stats_sieve_multivariate.R; driver stage sieve;
     quadprog is already in the renv library. Code comments use the
     symbols of vignettes/_notation.qmd.

     Model (per bulk b; unnormalised weights w = n p ≥ 0, p = w / Σw):

     - y_b ~ N(M w_b, D^½ V(w_b, M) D^½), with
       V = Σ_k w_k (Σ_k + m_k m_kᵀ) − (M w)(M w)ᵀ / Σ_k w_k. Mean and variance
       share n, as in RNA-Sieve. The between-type term is the multivariate
       form of the paper's eq. 7; the paper does not state it (my derivation,
       flagged in the code, with a switch).
     - Reference: μ̃_k ~ N(m_k, Σ_k / C_k). D = diag(d_g), start at 1.

     Stages, each stored (one chain per time × gene set × covariance
     input, over the 6 bulks):

     ┌──────────────┬──────┬───────────────────────────────────────────────────┐
     │    Stage     │ Free │                      Update                       │
     ├──────────────┼──────┼───────────────────────────────────────────────────┤
     │              │      │ freeze V, solve the QP min (y − Mw)ᵀV⁻¹(y − Mw),  │
     │ qp           │ w    │ w ≥ 0 (quadprog::solve.QP), recompute V; stop at  │
     │              │      │ max change in p < 1e-4 or 50 passes               │
     ├──────────────┼──────┼───────────────────────────────────────────────────┤
     │              │      │ add the closed-form mean update per bulk, M ≥ 0:  │
     │ qp_means     │ w, M │ m_k = μ̃_k + (w_k/C_k) Σ_k (V + Σ_k (w_k²/C_k)     │
     │              │      │ Σ_k)⁻¹ (y − M̃w)                                   │
     ├──────────────┼──────┼───────────────────────────────────────────────────┤
     │              │ w,   │ add d_g = mean over the 6 bulks of r_bg² / V_b,gg │
     │ qp_means_var │ M, d │  (your choice: the diagonal of the DECALS         │
     │              │      │ regression; correlations kept)                    │
     ├──────────────┼──────┼───────────────────────────────────────────────────┤
     │              │      │ from the previous state, maximise the full        │
     │ full_lik     │ w    │ likelihood with its log-determinant over w ≥ 0    │
     │              │      │ (L-BFGS-B, analytic gradient)                     │
     ├──────────────┼──────┼───────────────────────────────────────────────────┤
     │ qp_huber     │ w    │ as qp, with Huber weights (k = 1.345) on the      │
     │              │      │ whitened residuals, scale by MAD; reference fixed │
     └──────────────┴──────┴───────────────────────────────────────────────────┘

     - Covariance inputs: diagonal of raw Σ (the paper's independence
       model, in R), full raw Σ, full Σ_bio.
     - What each contrast tests: qp against DeCovarT L-BFGS-B with the
       same Σ → H1 (same plug-ins, QP estimating equation against profile
       likelihood). qp_means, qp_means_var against qp → H2.
       full_lik against qp_means_var → whether the log-determinant still
       moves the estimate once the scale is tied (finding 3). qp_huber
       against qp → whether a few genes drive the estimate.
     - The qp stage is also stored without the between-type term of V.
     - Departures from the paper, stated in the file header: means
       refined per bulk, not shared across bulks (a shared update is a
       35,000-dimensional dense system at 96 h); no +1 variance floor; no
       cell trimming or 20% detection filter (the panel is already gated).
     - Cost: about 15 s per Cholesky at 5,821 genes on one thread; a 96 h
       Welch chain is a few hundred factorisations. Three chains run as
       separate processes with 8 BLAS threads each; expected 1 to 2 h wall
       for that slot, minutes elsewhere.

     Step 3b · Shrinkage sweep for DeCovarT

     Driver stage lambda. Raw covariance with λ ∈ {closed form, 0.25, 0.5,
     0.75, 1}: off-diagonal terms of sigma_raw[,,j] rescaled by
     (1 − λ)/(1 − λ*_j), variances untouched. DeCovarT L-BFGS-B, bulk as is,
     both panels, three times (the closed-form level is already on disk).
     At 96 h on the Welch panel: 24 fits of about 9 min, about 20 min wall.
     This is the first slice of the sweep announced in the docs/06 TODO.

     Step 3c · Correlation between mean profiles

     plot_profile_correlation() in R/04_01_vis_deconvolution.R:
     geom_tile(), columns = 48, 72, 96 h, rows = Pearson (raw mean counts,
     the scale the solvers see) and Kendall. Cell types ordered by
     hierarchical clustering within a time, the coefficient printed in each
     tile, one sequential scale shared by all panels with its limits stated.
     One file per gene set (profile_correlation_<gene_set>), values in
     output/deconvolution/profile_correlation.csv.

     Step 3d · DeCovarT with random cell-type counts

     Driver stage total; code in R/04_02_stats_sieve_multivariate.R
     (the package is not changed). Same profile likelihood as
     deconvolute_single_cell() (free mean scale c and variance scale φ,
     both in closed form), with the covariance of a bulk whose cell-type
     counts are multinomial:

     S_T(p) = Σ_k p_k (Σ_k + m_k m_kᵀ) − (Mp)(Mp)ᵀ.

     - L-BFGS-B on non-negative weights renormalised inside the objective
       (the package's own closure), analytic gradient. The core (mixture,
       factorisation, gradient) is shared with the full_lik stage of
       step 3, with two flags: between-type term on or off, scale tied or
       profiled.
     - Regression test: with the term off and the scales profiled, the
       function must reproduce the archived DeCovarT L-BFGS-B log-likelihood
       and estimate.
     - Inputs: raw Σ and Σ_bio; bulk as is and divided by b_g. Four fits
       per sample and panel; about 4 CPU-hours at 96 h on the Welch panel.
     - Label "DeCovarT total covariance"; one more row block in the figures.

     Step 3e · Firth-penalised DeCovarT

     Driver stage firth. DeCovarT::fit_firth_single_cell() maximises the
     profile log-likelihood plus ½ log det of the observed information of
     the K − 1 ILR coordinates (BFGS in ILR, penalty gradient by central
     differences).

     - What to expect. The penalty is the determinant of a (K − 1) ×
       (K − 1) matrix, a few units; the term found in finding 2 is the
       G × G log det S(p), about a thousand units. Firth should pull
       estimates off the faces (the penalty is −∞ where the information is
       singular) but is not expected to offset the covariance
       log-determinant. The run tests that expectation.
     - Start. The package function starts from its own Newton fit and
       has no start argument; a start on the floor makes the penalty −∞ and
       the fit cannot move. The wrapper therefore rebuilds the same
       objective from the exported profile_loglik_single_cell_ilr(),
       gradient_… and hessian_profile_loglik_single_cell_ilr(), with two
       interior starts (the L-BFGS-B estimate floored at 0.001, and the
       barycentre), and keeps the higher penalised value. On a toy it must
       agree with fit_firth_single_cell().
     - Cost gate. One gradient needs 2(K − 1) Hessians. Small panels at
       the three times and the Welch panel at 48 h run in full. At 72 h and
       96 h on the Welch panel one Hessian is timed first; the slot runs
       only if a fit is projected under 1 h, otherwise it is reported as not
       run with the timing.
     - Inputs: raw Σ and Σ_bio, bulk as is and divided by b_g. Label
       "DeCovarT Firth".

     Remark · quadratic (population) covariance, not implemented

     The vignette's three-layer law is Σ_j p_j² Σ_B,j + (1/n) Σ_j p_j Σ_W,j
     - Σ_obs. The fitted n̂ ≈ 0.05 says the linear layer cannot carry the
       residual: what is missing does not shrink with n, so it belongs to
       Σ_B (quadratic) or Σ_obs (diagonal). Σ_B needs R ≥ 2 replicates; the
       reference has one cell line at 72 h and 96 h, and putting Σ_W in the
       quadratic slot silences rare types (p² against p), as the vignette
       warns. The gene-wise d_g of step 3 estimates the Σ_obs layer instead.
       The chapter states this as the next model to try at 48 h (two cell
       lines, scalar separation only); no fit in this pass.

     Step 4 · Tables, figures, report

     - stages=comparators,sieve,lambda,total,firth,tables, then
       scripts/04_01_vis_deconvolution.R. Main figures gain rows for MEAD,
       ReDeconv and the three sieve inputs (stage qp_means_var). New
       sieve_stages_{rmse,aitchison}: rows = covariance input, columns =
       time, one line per gene set, colour = stage.
     - Review each PNG against the data-visualisation checklist.
     - Correct convergence_failures.md (72 h Welch has 17 / 48 Newton
       failures, which breaks its gene-count trend) and point it to the
       chapter.

     Step 5 · Report in docs/04-deconvolution-pipeline.qmd

     Under scenario 2 (#sec-deconv-fixed-panel), replacing its TODO:

     1. Results on a shared panel: the three metric figures and the
        scatter, numbers read from output/deconvolution/*.csv in chunks.
     2. Why the profile likelihood loses at 96 h, new subsection
        {#sec-deconv-96h-diagnosis}:
        - the preliminary assumptions, each with its verdict and evidence:
          yours (too many or too few genes; Gaussian law; covariance not
          regularised enough; ILR against QP; bulk moments differ from the
          plug-ins) and mine from convergence_failures.md (ILR chart
          trap), kept even where the data went against them;
        - the estimator table above over the 18 samples, the likelihood
          decomposition (quadratic against log-determinant), the sieve
          stages, the λ sweep, the profile-correlation heatmap;
        - the RNA-Sieve implementation finding, stated as read from the
          code;
        - limits: one reference per time, 3 replicates, scored-cell
          fractions as truth, pseudo-individuals for MEAD.
     3. Speculation is labelled. Every display equation carries a label and
        a ledger id from step 6; the chapter prints the ledger table from
        verification/ledger.json. British English, {.highlight} for open
        items, no rendering of the whole book unless you ask (I render this
        chapter alone to check the chunks).

     Step 6 · Verification ledger (SymPy)

     New tree, following .claude/rules/symbolic-verification.md:
     verification/sympy/check_<group>_*.py (PEP 723, sympy==1.14.0,
     uv run --script), verification/run_all.py → verification/ledger.json
     and report/verification-ledger.tex; exit 0 only without fail. Each
     check derives the quantity from first principles on small symbolic
     cases (2 to 4 genes, 2 to 3 types) and compares it with the formula
     used in the code and the chapter. Symbols are compared with the shared
     notation file first (its path is looked up; the one named in the rule
     does not exist at that place).

     ┌─────┬────────────────────────────────────────────┬──────────────────────┐
     │ Id  │                   Claim                    │         Kind         │
     ├─────┼────────────────────────────────────────────┼──────────────────────┤
     │ C1  │ Bulk covariance with multinomial counts =  │ exact                │
     │     │ n S_T(p), from Cov(N) = n(diag p − ppᵀ)    │                      │
     ├─────┼────────────────────────────────────────────┼──────────────────────┤
     │ C2  │ S_T(p) − Σ_k p_k Σ_k = Σ_k p_k (m_k −      │ exact                │
     │     │ Mp)(m_k − Mp)ᵀ, rank ≤ K − 1               │                      │
     ├─────┼────────────────────────────────────────────┼──────────────────────┤
     │ C3  │ det and inverse of T − mmᵀ (rank-one       │ exact                │
     │     │ update)                                    │                      │
     ├─────┼────────────────────────────────────────────┼──────────────────────┤
     │     │ Closed-form ĉ and φ̂ are the stationary     │                      │
     │ P1  │ point; profile log-likelihood = −(G/2) log │ exact                │
     │     │  φ̂ − ½ log det S + const                   │                      │
     ├─────┼────────────────────────────────────────────┼──────────────────────┤
     │     │ Gradient of the profile log-likelihood in  │ exact, else          │
     │ P2  │ p_k, with and without the between-type     │ rational-sampling    │
     │     │ term                                       │                      │
     ├─────┼────────────────────────────────────────────┼──────────────────────┤
     │     │ Profiled against tied scale: the gradient  │                      │
     │ P3  │ of the quadratic term is multiplied by     │ exact                │
     │     │ φ̂⁻¹-free factor G/q against 1/n            │                      │
     ├─────┼────────────────────────────────────────────┼──────────────────────┤
     │ L1  │ log det(Σ_k p_k Σ_k) ≥ Σ_k p_k log det Σ_k │ numeric (mpmath),    │
     │     │                                            │ not a proof          │
     ├─────┼────────────────────────────────────────────┼──────────────────────┤
     │     │ Two rank-deficient matrices plus εI: the   │                      │
     │ L2  │ gap of L1 grows as −log ε per filled       │ series               │
     │     │ direction                                  │                      │
     ├─────┼────────────────────────────────────────────┼──────────────────────┤
     │ S1  │ Mean update of stage qp_means minimises    │ exact                │
     │     │ the stage objective (Woodbury form)        │                      │
     ├─────┼────────────────────────────────────────────┼──────────────────────┤
     │ S2  │ Rescaling off-diagonal terms by (1 − λ)/(1 │ exact                │
     │     │  − λ*) gives the λ-shrunk covariance       │                      │
     ├─────┼────────────────────────────────────────────┼──────────────────────┤
     │ M1  │ Aitchison distance = Euclidean distance of │ exact                │
     │     │  ILR (Helmert) coordinates                 │                      │
     ├─────┼────────────────────────────────────────────┼──────────────────────┤
     │ M2  │ Transcript to cell fractions: p_j ∝ q_j /  │ exact                │
     │     │ s_j inverts q_j ∝ p_j s_j                  │                      │
     ├─────┼────────────────────────────────────────────┼──────────────────────┤
     │     │ Firth penalty: ILR information = Jᵀ I_p J  │                      │
     │ F1  │ with J = (diag p − ppᵀ)V; it vanishes on a │ exact                │
     │     │  face                                      │                      │
     ├─────┼────────────────────────────────────────────┼──────────────────────┤
     │ Q1  │ Three-layer law: Cov = Σ p_j² Σ_B,j +      │ exact                │
     │     │ (1/n) Σ p_j Σ_W,j from the cell hierarchy  │                      │
     └─────┴────────────────────────────────────────────┴──────────────────────┘

     No Lean: this repository has no Lean build project, and the rule forbids
     reusing another project's. run_all.py records Lean as not configured.
     Probabilistic claims (C1) are verified conditional on the stated
     moments. A check that fails is reported with the value returned and the
     claim is rewritten.

     Verification

     1. Synthetic 3-type toy (DeCovarT::simulate_bulk_from_cells()): every
        sieve stage recovers p; the analytic gradient of full_lik matches
        finite differences; diagonal qp agrees with the Python RNA-Sieve
        projection on the same inputs to within its tolerance.
     2. Aitchison: zero for equal compositions; equals the distance between
        DeCovarT::isometric_log_ratio() vectors.
     3. MEAD: conversion to cell fractions checked on the toy.
     4. Smoke run: 48 h, small panel, all new methods; then the full run
        under nohup with a log in logs/; comparator_status.csv all
        ok, or the failure reported.

     Not in this pass

     No commit. No change to the DeCovarT package. No intervals for the new
     estimators. No eigenvalue-floor arm and no graphical-lasso arm (they
     stay in the docs/06 TODO).

  Read 3 files, ran 4 shell commands (ctrl+o to expand)

● Write(R/04_02_stats_sieve_multivariate.R)
  ⎿  Wrote 954 lines to R/04_02_stats_sieve_multivariate.R
       1 #' Estimators built on the Gaussian bulk law with a mixture cova
         riance.
       2 #'
       3 #' Notation (DeCovarT vignettes): `G` genes, `J` cell types, `M`
          the
       4 #' `G x J` matrix of mean counts per cell (`mu`), `Sigma_j` the
       5 #' within-type covariance, `C_j` the number of reference cells,
         `p` the
       6 #' cellular proportions, `w = n p` unnormalised weights with `n`
          the
     … +948 lines (ctrl+o to expand)

  Read 1 file, ran 9 shell commands (ctrl+o to expand)
● Update(R/04_01_stats_deconvolution.R)
  ⎿  Added 215 lines
      685    )
      686  }
      687  
      688 +#' MEAD with pseudo-individuals.                              
      689 +#'                                                            
      690 +#' `MEAD::MEAD()` on the cell-level raw counts of the panel, a
          +s in the                                                      
      691 +#' DeCovarT competitor runner (`get_R01(alpha = 0.1)` for the 
          +gene-gene                                                     
      692 +#' dependence set, falling back to `R01 = NULL`; gene filter o
          +ff, so                                                        
      693 +#' every method sees the same panel; other settings at their d
          +efaults).                                                     
      694 +#'                                                            
      695 +#' MEAD weights genes by their variance between reference indi
          +viduals.                                                      
      696 +#' The Suppinger reference has one cell line at 72 h and 96 h,
          + so the                                                       
      697 +#' cells of each type are split at random into `n_pseudo` grou
          +ps. Those                                                     
      698 +#' groups differ by sampling noise only: the weights measure t
          +he                                                            
      699 +#' precision of the mean, not donor variation.                
      700 +#'                                                            
      701 +#' MEAD divides each cell by its total over the supplied genes
          +, so its                                                      
      702 +#' estimate is the share of panel transcripts, `q`. It is conv
          +erted to                                                      
      703 +#' cell proportions with the mean panel total per cell `s_j`: 
      704 +#' `p_j = (q_j / s_j) / sum_k (q_k / s_k)`. MEAD itself drops 
          +the cells                                                     
      705 +#' with no count on the panel.                                
      706 +#'                                                            
      707 +#' @param bulk Bulk raw counts, genes by samples.             
      708 +#' @param counts Reference raw counts, cells by genes.        
      709 +#' @param cell_type Cell-type label of each reference cell.   
      710 +#' @param mu Mean signature, genes by types (gives `s_j`).    
      711 +#' @param n_pseudo Number of pseudo-individuals. Default `4`. 
      712 +#' @param seed Seed of the random split.                      
      713 +#' @return Long table, see [proportions_long()], with the attr
          +ibute                                                         
      714 +#'   `note` (route taken for `R01`, cells without counts).    
      715 +deconvolute_mead <- function(                                 
      716 +  bulk,                                                       
      717 +  counts,                                                     
      718 +  cell_type,                                                  
      719 +  mu,                                                         
      720 +  n_pseudo = 4L,                                              
      721 +  seed = 20261110L                                            
      722 +) {                                                           
      723 +  keep <- !is.na(cell_type)                                   
      724 +  counts <- counts[keep, , drop = FALSE]                      
      725 +  cell_type <- cell_type[keep]                                
      726 +  genes <- colnames(counts)                                   
      727 +  types <- colnames(mu)                                       
      728 +  bulk <- as.matrix(bulk)[genes, , drop = FALSE]              
      729 +  empty <- sum(rowSums(counts) == 0)                          
      730 +  individual <- withr::with_seed(seed, {                      
      731 +    out <- integer(length(cell_type))                         
      732 +    for (ty in types) {                                       
      733 +      idx <- which(cell_type == ty)                           
      734 +      out[idx] <- sample(rep_len(seq_len(n_pseudo), length(idx
          +)))                                                           
      735 +    }                                                         
      736 +    out                                                       
      737 +  })                                                          
      738 +  ref_counts <- t(counts)                                     
      739 +  colnames(ref_counts) <- paste0("cell_", seq_len(ncol(ref_cou
          +nts)))                                                        
      740 +  ref <- SingleCellExperiment::SingleCellExperiment(          
      741 +    assays = list(counts = ref_counts),                       
      742 +    colData = S4Vectors::DataFrame(                           
      743 +      cell_type = factor(cell_type, levels = types),          
      744 +      individual = factor(paste0("pseudo_", individual))      
      745 +    )                                                         
      746 +  )                                                           
      747 +  r01 <- tryCatch(MEAD::get_R01(bulk, alpha = 0.1), error = fu
          +nction(e) e)                                                  
      748 +  route <- "R01 from get_R01(alpha = 0.1)"                    
      749 +  if (inherits(r01, "error")) {                               
      750 +    route <- paste0("R01 = NULL (get_R01: ", conditionMessage(
          +r01), ")")                                                    
      751 +    r01 <- NULL                                               
      752 +  }                                                           
      753 +  fit <- suppressMessages(MEAD::MEAD(                         
      754 +    bulk,                                                     
      755 +    ref,                                                      
      756 +    cell_types = types,                                       
      757 +    preprocessing_control = list(filter.gene = FALSE),        
      758 +    R01 = r01                                                 
      759 +  ))                                                          
      760 +  transcript <- t(fit$p_hat[types, colnames(bulk), drop = FALS
          +E])                                                           
      761 +  size <- colSums(mu)[types]                                  
      762 +  cellular <- transcript / rep(size, each = nrow(transcript)) 
      763 +  cellular <- cellular / rowSums(cellular)                    
      764 +  out <- proportions_long(cellular)                           
      765 +  attr(out, "note") <- paste0(                                
      766 +    route,                                                    
      767 +    "; ",                                                     
      768 +    empty,                                                    
      769 +    " reference cells without a count on the panel"           
      770 +  )                                                           
      771 +  out                                                         
      772 +}                                                             
      773 +                                                              
      774 +#' ReDeconv through its Python implementation.                
      775 +#'                                                            
      776 +#' Writes the bulk and the reference moments (mean and standar
          +d                                                             
      777 +#' deviation of the raw counts per cell, per type), runs      
      778 +#' `scripts/python/run_redeconv.py` with `uv run --script`, an
          +d reads the                                                   
      779 +#' percentages. Supplying the moments bypasses ReDeconv's own 
      780 +#' signature-gene selection and its transcriptome-size normali
          +sation of                                                     
      781 +#' the reference: every panel gene is used, as for the other m
          +ethods.                                                       
      782 +#' ReDeconv clips negative coefficients and rescales the posit
          +ive ones;                                                     
      783 +#' the estimates are rescaled here to sum to one.             
      784 +#'                                                            
      785 +#' @param bulk Bulk raw counts, genes by samples.             
      786 +#' @param mu Mean signature, genes by types.                  
      787 +#' @param sd Standard deviations, genes by types.             
      788 +#' @param work_dir Folder for the exchange files (created, kep
          +t).                                                           
      789 +#' @param script Path to the Python runner.                   
      790 +#' @return Long table, see [proportions_long()].              
      791 +deconvolute_redeconv <- function(                             
      792 +  bulk,                                                       
      793 +  mu,                                                         
      794 +  sd,                                                         
      795 +  work_dir,                                                   
      796 +  script = file.path("scripts", "python", "run_redeconv.py")  
      797 +) {                                                           
      798 +  genes <- rownames(mu)                                       
      799 +  types <- colnames(mu)                                       
      800 +  in_dir <- file.path(work_dir, "input")                      
      801 +  out_dir <- file.path(work_dir, "output")                    
      802 +  unlink(c(in_dir, out_dir), recursive = TRUE)                
      803 +  dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)  
      804 +  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE) 
      805 +  bulk <- as.matrix(bulk)[genes, , drop = FALSE]              
      806 +  # ReDeconv reads tab-separated files and splits type names o
          +n "_".                                                        
      807 +  codes <- paste0("type", seq_along(types))                   
      808 +  write_tsv <- function(m, path) {                            
      809 +    data.table::fwrite(                                       
      810 +      data.table::data.table(gene = rownames(m), m),          
      811 +      path,                                                   
      812 +      sep = "\t"                                              
      813 +    )                                                         
      814 +  }                                                           
      815 +  write_tsv(bulk, file.path(in_dir, "bulk.tsv"))              
      816 +  signature <- matrix(                                        
      817 +    0,                                                        
      818 +    length(genes),                                            
      819 +    2L * length(types),                                       
      820 +    dimnames = list(                                          
      821 +      genes,                                                  
      822 +      as.vector(rbind(paste0(codes, "_mean"), paste0(codes, "_
          +std")))                                                       
      823 +    )                                                         
      824 +  )                                                           
      825 +  signature[, paste0(codes, "_mean")] <- mu[genes, types]     
      826 +  signature[, paste0(codes, "_std")] <- sd[genes, types]      
      827 +  write_tsv(signature, file.path(in_dir, "signature.tsv"))    
      828 +  status <- system2(                                          
      829 +    "uv",                                                     
      830 +    c("run", "--script", shQuote(script), shQuote(in_dir), shQ
          +uote(out_dir))                                                
      831 +  )                                                           
      832 +  if (status != 0L) {                                         
      833 +    stop("ReDeconv runner failed with status ", status, ".")  
      834 +  }                                                           
      835 +  tab <- utils::read.csv(                                     
      836 +    file.path(out_dir, "estimates.csv"),                      
      837 +    check.names = FALSE,                                      
      838 +    stringsAsFactors = FALSE,                                 
      839 +    fileEncoding = "UTF-8"                                    
      840 +  )                                                           
      841 +  est <- as.matrix(tab[, -1L, drop = FALSE])                  
      842 +  rownames(est) <- tab[[1L]]                                  
      843 +  est <- t(est)[colnames(bulk), codes, drop = FALSE]          
      844 +  colnames(est) <- types                                      
      845 +  if (any(rowSums(est) <= 0)) {                               
      846 +    stop("ReDeconv returned no positive coefficient for a samp
          +le.")                                                         
      847 +  }                                                           
      848 +  proportions_long(est / rowSums(est))                        
      849 +}                                                             
      850 +                                                              
      851 +#' Robust linear regression of the bulk on the signature.     
      852 +#'                                                            
      853 +#' `DeCovarT::deconvolute_ratios_rlm()`: Huber M-estimation   
      854 +#' (`MASS::rlm()`), no intercept, no covariance, no sign const
          +raint;                                                        
      855 +#' the coefficients are then divided by their sum. That functi
          +on stops                                                      
      856 +#' when a coefficient is negative. In that case the same `MASS
          +::rlm()`                                                      
      857 +#' fit is repeated here, negative coefficients are set to zero
          + and the                                                      
      858 +#' rest rescaled; the attribute `note` counts those samples.  
      859 +#'                                                            
      860 +#' @param bulk Bulk raw counts, genes by samples.             
      861 +#' @param mu Mean signature, genes by types.                  
      862 +#' @return Long table, see [proportions_long()], with the attr
          +ibute                                                         
      863 +#'   `note`.                                                  
      864 +deconvolute_rlm <- function(bulk, mu) {                       
      865 +  genes <- rownames(mu)                                       
      866 +  bulk <- as.matrix(bulk)[genes, , drop = FALSE]              
      867 +  clipped <- 0L                                               
      868 +  est <- vapply(                                              
      869 +    colnames(bulk),                                           
      870 +    function(s) {                                             
      871 +      y <- bulk[, s]                                          
      872 +      p <- tryCatch(                                          
      873 +        DeCovarT::deconvolute_ratios_rlm(y, mu),              
      874 +        error = function(e) NULL                              
      875 +      )                                                       
      876 +      if (is.null(p)) {                                       
      877 +        clipped <<- clipped + 1L                              
      878 +        p <- pmax(MASS::rlm(y ~ mu + 0, method = "M")$coeffici
          +ents, 0)                                                      
      879 +        if (sum(p) <= 0) {                                    
      880 +          stop("Robust regression gave no positive coefficient
          + for ", s, ".")                                               
      881 +        }                                                     
      882 +        p <- p / sum(p)                                       
      883 +      }                                                       
      884 +      as.numeric(p)                                           
      885 +    },                                                        
      886 +    numeric(ncol(mu))                                         
      887 +  )                                                           
      888 +  rownames(est) <- colnames(mu)                               
      889 +  out <- proportions_long(t(est))                             
      890 +  attr(out, "note") <- if (clipped > 0L) {                    
      891 +    paste0(                                                   
      892 +      clipped,                                                
      893 +      " of ",                                                 
      894 +      ncol(bulk),                                             
      895 +      " samples had a negative coefficient, set to zero"      
      896 +    )                                                         
      897 +  } else {                                                    
      898 +    ""                                                        
      899 +  }                                                           
      900 +  out                                                         
      901 +}                                                             
      902 +                                                              
      903  # ============================================================
           ==============
      904  # Benchmark metrics ----
      905  # ============================================================
           ==============

  Ran 2 shell commands (ctrl+o to expand)

● Write(scripts/python/run_redeconv.py)
  ⎿  Wrote 72 lines to scripts/python/run_redeconv.py
       1 # /// script
       2 # requires-python = ">=3.11,<3.12"
       3 # dependencies = [
       4 #     "redeconv==1.2.1",
       5 #     "pandas==2.0.3",
       6 # ]
       7 # ///
       8 """Run the official ReDeconv implementation on one signature and
          its bulks.
       9
     … +63 lines (ctrl+o to expand)

● Update(R/04_01_stats_deconvolution.R)
  ⎿  Added 2 lines, removed 1 line
      803    dir.create(in_dir, recursive = TRUE, showWarnings = FALSE)
      804    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      805    bulk <- as.matrix(bulk)[genes, , drop = FALSE]
      806 -  # ReDeconv reads tab-separated files and splits type names o
          -n "_".                                                        
      806 +  # ReDeconv reads tab-separated files and takes the type name
          + from the                                                     
      807 +  # column header; plain codes avoid spaces and accents there.
      808    codes <- paste0("type", seq_along(types))
      809    write_tsv <- function(m, path) {
      810      data.table::fwrite(

● Update(R/04_01_stats_deconvolution.R)
  ⎿  Added 18 lines, removed 21 lines
      865  deconvolute_rlm <- function(bulk, mu) {
      866    genes <- rownames(mu)
      867    bulk <- as.matrix(bulk)[genes, , drop = FALSE]
      868 -  clipped <- 0L                                               
      869 -  est <- vapply(                                              
      870 -    colnames(bulk),                                           
      871 -    function(s) {                                             
      872 -      y <- bulk[, s]                                          
      873 -      p <- tryCatch(                                          
      874 -        DeCovarT::deconvolute_ratios_rlm(y, mu),              
      875 -        error = function(e) NULL                              
      876 -      )                                                       
      877 -      if (is.null(p)) {                                       
      878 -        clipped <<- clipped + 1L                              
      879 -        p <- pmax(MASS::rlm(y ~ mu + 0, method = "M")$coeffici
          -ents, 0)                                                      
      880 -        if (sum(p) <= 0) {                                    
      881 -          stop("Robust regression gave no positive coefficient
          - for ", s, ".")                                               
      882 -        }                                                     
      883 -        p <- p / sum(p)                                       
      868 +  fits <- lapply(colnames(bulk), function(s) {                
      869 +    y <- bulk[, s]                                            
      870 +    p <- tryCatch(                                            
      871 +      DeCovarT::deconvolute_ratios_rlm(y, mu),                
      872 +      error = function(e) NULL                                
      873 +    )                                                         
      874 +    was_clipped <- is.null(p)                                 
      875 +    if (was_clipped) {                                        
      876 +      p <- pmax(MASS::rlm(y ~ mu + 0, method = "M")$coefficien
          +ts, 0)                                                        
      877 +      if (sum(p) <= 0) {                                      
      878 +        stop("Robust regression gave no positive coefficient f
          +or ", s, ".")                                                 
      879        }
      885 -      as.numeric(p)                                           
      886 -    },                                                        
      887 -    numeric(ncol(mu))                                         
      888 -  )                                                           
      889 -  rownames(est) <- colnames(mu)                               
      880 +      p <- p / sum(p)                                         
      881 +    }                                                         
      882 +    list(p = as.numeric(p), clipped = was_clipped)            
      883 +  })                                                          
      884 +  est <- vapply(fits, `[[`, numeric(ncol(mu)), "p")           
      885 +  dimnames(est) <- list(colnames(mu), colnames(bulk))         
      886 +  clipped <- sum(vapply(fits, `[[`, logical(1), "clipped"))   
      887    out <- proportions_long(t(est))
      888    attr(out, "note") <- if (clipped > 0L) {
      889      paste0(

  Read 1 file, ran 13 shell commands (ctrl+o to expand)

● Now the figure code: labels for the new methods, the generalised metric plot,
  and three new plots.

  Read 3 files, ran 13 shell commands (ctrl+o to expand)

● The second-pass runs are in progress; I will continue as soon as they finish.

  - Done so far: MEAD, ReDeconv and robust regression ran at the three times on
    both panels (all ok); every new stage ran on the small panels and on the
    Welch panel at 48 h; the code for the Aitchison distance, the new figures
    and the SymPy ledger is written, and the ledger passes (14 rows, 0 failed).
  - Still running: the Welch panel at 72 h and 96 h for the multivariate
    RNA-Sieve, the shrinkage sweep, the total-covariance DeCovarT and the Firth
    fits. A background waiter will notify me when the five processes exit.
  - Then: rebuild tables and figures, compute the likelihood split over the
    standard-Wnt samples, review the figures, and write the report into
    chapter 4.

❯ when do you plan stopping? what's remainign to eprform for the background
  agent?

  Ran 3 shell commands (ctrl+o to expand)

● I stop once the report is in chapter 4, which I expect around 18:30 to 19:00.
  The model runs are the slow part; they should end around 17:30 to 18:00.

  What is still running (16:10)

  ┌───────────────────────┬───────────────────────┬─────────────────────────┐
  │          Run          │         State         │          Left           │
  ├───────────────────────┼───────────────────────┼─────────────────────────┤
  │ Shrinkage sweep,      │ level 0.25 done in 28 │ levels 0.5, 0.75 and 1, │
  │ Welch panel, 96 h     │  min                  │  about 1 h 15           │
  ├───────────────────────┼───────────────────────┼─────────────────────────┤
  │ DeCovarT with total   │ just started, in      │ 24 fits, about 1 h to 1 │
  │ covariance, Welch     │ parallel with the     │  h 30                   │
  │ panel, 96 h           │ sweep                 │                         │
  ├───────────────────────┼───────────────────────┼─────────────────────────┤
  │ Firth, Welch panel,   │ running (one Hessian  │ 24 fits on 16 workers,  │
  │ 72 h                  │ 11 s, one fit         │ up to 1 h 30            │
  │                       │ projected at 45 min)  │                         │
  ├───────────────────────┼───────────────────────┼─────────────────────────┤
  │ Firth, Welch panel,   │                       │ most likely skipped by  │
  │ 96 h                  │ not started           │ the 1 h cost gate, with │
  │                       │                       │  its timing recorded    │
  └───────────────────────┴───────────────────────┴─────────────────────────┘

  Everything else has finished:

  - Comparators: MEAD, ReDeconv and robust regression, at the three times on
    both panels.
  - Multivariate RNA-Sieve: all 18 chains, including 96 h on the Welch panel (11
    to 20 min each).
  - Shrinkage sweep, total covariance and Firth: on the small panels at the
    three times, and on the Welch panel at 48 h; the sweep and total covariance
    also at 72 h.

  What I do after the runs

  1. Rebuild the tables and the 16 figures, and check each figure.
  2. Compute the likelihood split (quadratic term against log-determinant) on
     the nine standard-Wnt samples.
  3. Write the results and the diagnosis into
     docs/04-deconvolution-pipeline.qmd, render that chapter alone, and correct
     convergence_failures.md.

  Code, the Aitchison metric and the SymPy ledger (14 checks, all pass) are
  already written. Nothing is committed.

  One early result to expect in the report: refining the mean profiles from a
  single bulk makes the multivariate RNA-Sieve worse, not better (RMSE rises to
  about 0.75 on several panels). I implemented the refinement per bulk, not
  shared across bulks, so this may reflect that choice.

● Background command "Wait for the five large-panel runs to exit, then list what
they logged" completed (exit code 0)

  Ran 2 shell commands (ctrl+o to expand)

● Background command "Wait for the three remaining runs to exit" completed (exit
code 0)

  Read 3 files, ran 10 shell commands (ctrl+o to expand)

● I found a problem and am fixing it before reporting anything: some checkpoints
  were written while I was editing the driver script, and R reads a script file
  as it runs, so those blocks may have executed altered code. Eight shrinkage
  checkpoints were visibly malformed.

  - Set aside (not deleted): 129 checkpoints from the affected blocks, all on
    the Welch panel: the shrinkage sweep at 48 h and 72 h, the total-covariance
    fits at the three times, and the Firth fits at 48 h and 72 h.
  - Not affected: the comparators, the 18 multivariate RNA-Sieve chains, the 96
    h shrinkage sweep and everything on the small panels. Those blocks were read
    before the edit or ran in separate processes.
  - Refit in progress: from a frozen copy of the script, so further edits cannot
    touch it. I expect it to finish around 18:15; the Firth fits at 72 h are
    the slow part. I will compare the refits with the set-aside files to see
    whether the edit changed any number.

  Results that are already final, as mean RMSE over the three standard-Wnt
  replicates at 96 h on the Welch panel:

  ┌────────────────────────────────────────────────────────────────┬───────┐
  │                             Method                             │ RMSE  │
  ├────────────────────────────────────────────────────────────────┼───────┤
  │ MEAD (pseudo-individuals)                                      │ 0.048 │
  ├────────────────────────────────────────────────────────────────┼───────┤
  │ RNA-Sieve (Python)                                             │ 0.056 │
  ├────────────────────────────────────────────────────────────────┼───────┤
  │ Multivariate RNA-Sieve in R, QP with the full raw covariance   │ 0.060 │
  ├────────────────────────────────────────────────────────────────┼───────┤
  │ DeMixSC wNNLS                                                  │ 0.071 │
  ├────────────────────────────────────────────────────────────────┼───────┤
  │ ReDeconv                                                       │ 0.103 │
  ├────────────────────────────────────────────────────────────────┼───────┤
  │ DWLS                                                           │ 0.135 │
  ├────────────────────────────────────────────────────────────────┼───────┤
  │ BayesPrism                                                     │ 0.162 │
  ├────────────────────────────────────────────────────────────────┼───────┤
  │ Robust regression                                              │ 0.220 │
  ├────────────────────────────────────────────────────────────────┼───────┤
  │ SCDC                                                           │ 0.307 │
  ├────────────────────────────────────────────────────────────────┼───────┤
  │ Multivariate RNA-Sieve, same covariance, after mean refinement │ 0.503 │
  └────────────────────────────────────────────────────────────────┴───────┘

  - The QP result holds on all three replicates. With the same raw covariance
    plug-ins as DeCovarT, the quadratic programme converges and gives 0.060,
    where the first pass had about 0.36 for the DeCovarT L-BFGS-B fit.
  - It reverses on the small panel. At 96 h with 635 genes the same QP gives
    0.36 and does not converge.
  - Mean refinement hurts. Refining the mean profiles per bulk raises the error
    to about 0.5.
  - The denoised covariance alone is unstable. The QP iteration on it does not
    converge at 96 h.

  One change from the plan: the shared figures now show the multivariate
  RNA-Sieve at its QP stage, not at the stage with refined means and gene
  variances, because the refined stages are the degraded ones. All six stages
  remain in their own figure.

  Next, once the refit ends: final tables and figures, the likelihood split,
  then the report in chapter 4.

❯ state of the background jobs?

  Ran 1 shell command (ctrl+o to expand)

● One job is still running and it is nearly done: 22 of the 24 Firth refits at
  72 h on the Welch panel are on disk (17:41). The last two should finish within
  about 15 minutes.

  ┌─────────────────────────────────┬───────────────────────────────────────┐
  │               Job               │                 State                 │
  ├─────────────────────────────────┼───────────────────────────────────────┤
  │ Refit, Welch panel, 96 h: total │ finished, 24 of 24 fits               │
  │  covariance                     │                                       │
  ├─────────────────────────────────┼───────────────────────────────────────┤
  │                                 │ not run, by the cost gate: one        │
  │ Refit, Welch panel, 96 h: Firth │ Hessian takes 73 s, so one fit is     │
  │                                 │ projected at about 8 h                │
  ├─────────────────────────────────┼───────────────────────────────────────┤
  │ Refit, Welch panel, 48 h and 72 │                                       │
  │  h: shrinkage sweep and total   │ finished                              │
  │ covariance                      │                                       │
  ├─────────────────────────────────┼───────────────────────────────────────┤
  │ Refit, Welch panel, 48 h: Firth │ finished, 24 of 24                    │
  ├─────────────────────────────────┼───────────────────────────────────────┤
  │ Refit, Welch panel, 72 h: Firth │ 22 of 24, running                     │
  └─────────────────────────────────┴───────────────────────────────────────┘

  No refit has logged a failure. A waiter will notify me when the last process
  exits. Then I compare the refits with the set-aside checkpoints, rebuild
  tables and figures, compute the likelihood split, and write the report into
  chapter 4.



The figures are final and the report is in chapter 4, with one large surprise: DeCovarT with the raw covariance shrunk at λ = 0.5 is the best method of the benchmark at 96 h, with RMSE 0.034 against 0.36 at the closed-form shrinkage. Your hypothesis that the covariance was not regularised enough is supported; I had ranked it last, and that was wrong.

That 0.034 is optimistic: λ = 0.5 is the best of four levels, picked on the same nine samples that score it.

## State of the jobs

- **Refits:** finished. All 128 refitted estimates are identical to the set-aside ones, so editing the script during the run changed no number.
- **Figures and tables:** rebuilt, 15 figures in `docs/images/deconvolution/`. I viewed the RMSE, shrinkage, sieve-stage, correlation and Welch scatter figures in their final form; the others were redrawn but not inspected.
- **Firth, Welch panel, 96 h:** started in the background, as you asked, restricted to the three standard-Wnt samples (raw and denoised covariance, bulk as is). It is six fits at about 8 h each, so expect results tomorrow morning. Firth equals the unpenalised fit to the third decimal at 72 h on the Welch panel, so I expect no change at 96 h.

## Main results (96 h, Welch panel, mean RMSE over three replicates)

| Method | RMSE |
| --- | ---: |
| DeCovarT, raw covariance, λ = 0.5 | 0.034 |
| MEAD (pseudo-individuals) | 0.048 |
| RNA-Sieve (Python) | 0.056 |
| Multivariate RNA-Sieve in R, QP stage, full raw covariance | 0.060 |
| DeMixSC wNNLS | 0.071 |
| DeCovarT, raw covariance, closed-form λ | 0.36 |
| DeCovarT, denoised covariance | 0.47 |

- **Your two hypotheses:**
  - The optimiser is not the cause: DeCovarT's estimate has a higher likelihood than the observed composition in all 36 cases.
  - The estimating equation is: a quadratic programme on the same plug-ins gives 0.060.
  - Bulk and reference do differ (log-scale correlation 0.48 to 0.83), but refining the means from one bulk made things worse (0.50).
- **My log-determinant explanation** holds only for the denoised covariance on the Welch panel; with the raw covariance the gain is in the quadratic term.
- **The error at 96 h is a split inside one lineage.** After summing caudal mesoderm, pre-somitic mesoderm and the somite front, twelve of thirteen methods are at 0.04 or below, including DeCovarT with the denoised covariance.
- **No consistent effect:** the between-type covariance term and the Firth penalty.
- **Small panels:** none of the gains above appear there.

## Five directions, ranked by this evidence

1. **Choose the covariance regularisation for the deconvolution**, without using the observed composition (pseudo-bulks of held-out cells); test one common level against one per type.
2. **Estimate the composition from the mean equation** (quasi-likelihood, GEE, Fisher scoring), with a jackknife over gene clusters for the variance. HC3 leaving one gene out would be too narrow, because genes are correlated.
3. **Treat close types as a lineage**: fit the lineage first, then the split. This helps every method.
4. **Add a measurement layer between bulk and reference**: a gene-wise mean–variance law across bulks, which is the marginal of a Poisson log-normal model, or a mixed model with a gene effect. My single variance scale per gene did not work.
5. **Carry the precision of the reference** (replicates or a delete-a-group jackknife) instead of refitting the reference from the bulk.

The chapter also has a table of which direction should matter most by regime: cosine similarity of the profiles, bulk-to-reference shift, cells per gene, and covariance changes between conditions. Only the first three were measured; the last cannot be, because the reference follows one protocol.
