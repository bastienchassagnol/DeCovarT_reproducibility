# Targeted recommendation for benchmarking single-cell RNA-seq gene selection for cell-type prediction

## Executive recommendation

Your proposed benchmark is well posed **provided the target is ordinary multiclass cell-type annotation**, i.e. every cell has exactly one mutually exclusive cell-type label. In that setting, I would benchmark your **two-stage supervised marker procedure**—gene-wise omnibus screening followed by cell-type-versus-grand-mean contrasts—against exactly three deliberately different competitors: **multinomial elastic-net / structured multinomial regression**, **random forest**, and **sparse PLS-DA**. The DeCovarT vignette independently converges on essentially this strategy: it proposes multinomial elastic net as a joint multiclass candidate generator, penalised discriminant analysis, group/tree-guided penalties, and donor-level stability selection with tuning repeated inside resamples. citeturn613510view0turn613510view1

For family A, my ranking is **(1) `glmnet` multinomial elastic net, (2) `msgl` multinomial sparse-group lasso, (3) `grpnet` multinomial group elastic net**. `glmnet` should be the primary comparator because it is the mature, canonical implementation, gives class-specific coefficient vectors, works very well when \(p\gg n\), and directly supports lasso, ridge and elastic net. `msgl` is the best fit when your co-expression modules should enter the penalty explicitly while retaining within-module sparsity. `grpnet` is an especially interesting newer alternative when you want whole-module group selection combined with a group-ridge component, although its shorter track record makes it a better structured sensitivity analysis than the principal comparator. citeturn236094search0turn236094search3turn368561search4

One terminology correction matters: by “ridge and L2 combined”, you almost certainly mean **L1 + L2**, i.e. elastic net. Ridge *is* the L2 penalty; elastic net combines L1 sparsity with L2 shrinkage. `glmnet` implements this as
\[
\lambda\{\alpha\|\beta\|_1+(1-\alpha)\|\beta\|_2^2/2\},
\]
with `alpha=1` lasso, `alpha=0` ridge and intermediate values elastic net. citeturn249502search0

A second conceptual distinction is even more important. A multinomial model assumes **one categorical \(Y\) per cell**. If “multi-labelled” literally means that one cell can simultaneously possess several non-mutually-exclusive labels, multinomial logistic regression and PLS-DA are not the correct likelihood formulation. A simple proper baseline is then one penalised binary logistic model per label (“binary relevance”); if label dependencies matter, classifier chains extend that construction by feeding earlier label predictions into later classifiers. citeturn521237search0 The remainder of this report therefore assumes the biologically much more usual **single-label, multiclass cell-type prediction** unless stated otherwise.

The following is the benchmark I would actually run:

| Family | Recommended implementation | What it tests against your marker method | Main gene-selection object | Principal limitation |
|---|---|---|---|---|
| **Your two-stage markers** | omnibus screen → type-vs-grand-mean contrast | Classical supervised marginal marker selection | gene/contrast significance or score | ignores multigene predictive redundancy unless explicitly corrected |
| **Penalised multinomial** | `glmnet`, plus `msgl` structured sensitivity analysis | joint multivariate discriminative selection | non-zero class coefficients / module groups | coefficients are regularised; naive SEs/p-values are invalid |
| **Tree ensemble** | `tidymodels` + `ranger`; XGBoost optional | nonlinearities and interactions without linear-score assumptions | permutation importance / selection frequency | importance is unstable under correlated genes |
| **sPLS-DA** | `mixOmics::splsda`; `mint.splsda` when a genuine study factor exists | low-dimensional supervised latent signature | sparse loading support (`keepX`) | sparsity is on latent loadings, not an elastic-net logistic likelihood |

The comparison should be about **predictive marker selection**, not about which procedure produces the most DE discoveries. A gene can be biologically differential but redundant for prediction, while a jointly useful multigene predictor need not have the smallest marginal DE \(p\)-value. Penalised regression and sPLS-DA are explicitly joint supervised models, whereas your omnibus/contrast screen is fundamentally a marginal gene-wise procedure.

## Problem formulation and benchmark design

### Count predictors are not count responses

Your negative-binomial/quasi-Poisson requirement needs to be reformulated slightly. In a cell-type classifier, the response \(Y\) is the **categorical cell type** and the RNA counts are predictors \(X_g\). A negative-binomial likelihood would be appropriate if the **response itself were a gene count**, but it cannot replace the multinomial likelihood for categorical \(Y\). Your teaching resource makes exactly the underlying GLM distinction: Poisson and negative-binomial families describe count-valued *responses*, whereas Bernoulli/binomial families describe categorical responses. citeturn613510view2

Consequently, I would not search for a “negative-binomial multinomial classifier”. Instead, separate the two layers:

\[
\text{raw UMI counts}
\;\xrightarrow{\text{training-fold normalisation}}\;
X^*
\;\xrightarrow{\text{multinomial / RF / sPLS-DA}}\;
P(Y=k\mid X^*).
\]

`sctransform` is defensible for the first arrow: its model uses regularised negative-binomial regression to account for sequencing depth and returns variance-stabilised Pearson residuals intended for downstream multivariate analysis. citeturn529180search8 A simpler log-normalised expression matrix is worth retaining as a sensitivity analysis. The critical benchmark requirement is that all learned preprocessing parameters are estimated **inside the training portion of each outer split**, rather than once on all cells.

### The split unit matters more than the classifier

Do **not** make random cells the primary cross-validation unit when multiple cells originate from the same donor, biological replicate, culture or experimental unit. Recent large-scale single-cell benchmarks use donor-wise splitting precisely to reduce train–test leakage; scTab, for example, explicitly uses unseen-donor evaluation and macro-F1. citeturn902516search1turn902516search8 Recent direct investigations also show that random cell-level splitting can materially overestimate cross-donor cell-type performance relative to donor-held-out evaluation. citeturn902516search0

I would therefore define one outer split object and hand **exactly the same folds** to every family. Ideally, groups are donors/biological replicates. If cell line and time can themselves identify batches, the grouping constraint should also prevent the same experimental unit being represented on both sides of a split. The panel universe, preprocessing, feature selection, hyperparameter tuning and class-weight estimation all need to be redone within each training fold.

That last point includes your own method. The omnibus test and cell-type-specific contrasts use \(Y\), so they are **supervised feature selection**. Running them once before cross-validation would leak test labels into the gene panel.

I would make the primary endpoint **macro-F1 or balanced accuracy**, accompanied by per-class sensitivity/recall, confusion matrices, predicted-probability log loss or Brier-type calibration summaries, number of selected genes, and selection stability across outer resamples. Macro-F1 is especially appropriate when cell-type abundances differ substantially; it is also the metric used for donor-wise single-cell annotation benchmarking in scTab. citeturn902516search1

Finally, compare methods at **matched panel size** as well as at their native optimum. Otherwise one method selecting 30 genes and another retaining 800 genes is partly answering a different engineering problem. A useful presentation is performance versus panel size, complemented by Jaccard/selection-frequency stability and enrichment/coherence of selected modules.

## Penalised multinomial logistic regression

### Recommended first choice: `glmnet`

The strongest default comparator is **`glmnet` multinomial elastic net**. Friedman, Hastie and Tibshirani's original coordinate-descent paper covers multinomial regression with L1, L2 and their elastic-net mixture and was specifically designed for large sparse/high-dimensional problems. citeturn236094search0turn236094search1 As of 20 November 2026, CRAN `glmnet` is version **5.2**, released 28 October 2026. citeturn281258search3turn281258search19

The core model is

\[
P(Y_i=k\mid x_i)
=
\frac{\exp(\beta_{0k}+x_i^\top\beta_k)}
{\sum_{\ell=1}^K\exp(\beta_{0\ell}+x_i^\top\beta_\ell)} ,
\]

with a penalty controlled by `alpha` and `lambda`. This gives you a coefficient vector for every class, predicted probabilities, multinomial deviance and regularisation paths. citeturn249502search0turn236094search5

There are two useful multinomial penalty modes. `type.multinomial="ungrouped"` applies ordinary coefficient-wise sparsity; `type.multinomial="grouped"` groups the **coefficients of the same gene across all outcome classes**, so that a gene tends to enter or leave all class equations together. citeturn249502search8turn249502search9 This is an important distinction: `grouped` **does not mean biological gene modules**. It groups a gene's \(K\) class coefficients, not a set of co-expressed genes.

That is nevertheless attractive for your problem. With `ungrouped`, you obtain genuinely class-specific marker sparsity. With `grouped`, you ask whether a gene is informative somewhere in the multiclass system. I would benchmark both, but nominate **ungrouped elastic net as the main model if cell-type-specific interpretability is your first priority**.

`penalty.factor` also allows gene-specific shrinkage, including zero penalty for variables that must be retained. citeturn249502search9 It can therefore encode *soft* prior module information—for example, down-weighting genes in highly credible modules—but it cannot impose the all-in/all-out structure of arbitrary modules.

One subtle interpretability issue is worth documenting. `glmnet` fits a **symmetric multinomial parameterisation**, not a classical “class \(K\) is the reference” parameterisation; the redundancy is handled by the penalty. citeturn236094search5turn236094search7 Consequently, a single \(\beta_{gk}\) should be described as the contribution of gene \(g\) to class \(k\)'s linear score. The most invariant log-odds interpretation is through a class contrast,
\[
\log\frac{P(Y=k)}{P(Y=\ell)}
=
(\beta_{0k}-\beta_{0\ell})
+x^\top(\beta_k-\beta_\ell).
\]
That maps particularly nicely onto your “cell type versus overall background” motivation.

The main caveat is inference. `glmnet` gives penalised estimates, not ordinary maximum-likelihood standard errors and \(p\)-values. Naively treating the selected active set as if it had been specified in advance gives anti-conservative inference; specialised post-selection procedures exist for L1-penalised likelihood models, but they are a separate inferential layer. citeturn236094search6 For this benchmark I would therefore use coefficients, paths, stability frequencies, held-out calibration and prediction diagnostics as the principal diagnostics. If formal biological inference is required, keep your omnibus/contrast DE analysis as the inferential analysis rather than trying to manufacture conventional \(p\)-values from the classifier.

### Recommended structured option: `msgl`

My second choice is **`msgl`**, based on Vincent and Hansen's multinomial sparse-group lasso. Their 2014 paper explicitly formulates sparse group lasso for high-dimensional multiclass multinomial likelihood and reports an implementation in the R package `msgl`; their examples include problems with 10,000 predictors and many classes. citeturn236094search3 Current package documentation describes `msgl` 2.3.9 as multinomial logistic regression with sparse-group lasso and parallel CV/subsampling. citeturn355835search0

This is the most natural implementation of your module prior. With modules \(G_m\), it minimises a multinomial loss plus

\[
\lambda\left[
(1-\alpha)\sum_m \gamma_m\|\beta_{G_m}\|_2
+
\alpha\sum_j \xi_j|\beta_j|
\right].
\]

The first term encourages module/group selection; the second allows individual coefficients to disappear **inside an included module**. The software accepts an explicit feature `grouping` vector and group/parameter weights. citeturn355835search1turn355835search10

For your biology, I would define the groups from a **training-only fixed module map** if modules are data-derived; if modules come from external biology or a frozen prior network, that map may be fixed globally. One candidate grouping is co-expression modules irrespective of cell type; another sensitivity analysis is a biologically curated pathway/module map.

The principal limitation is overlapping membership. The standard `grouping` representation assigns a feature to one group; genes belonging to several modules require duplication/latent-group reformulation or custom code. `grpregOverlap` does implement overlapping grouped predictors, but its documented response families are linear, binary logistic, Poisson and Cox—not multinomial—so I would **not** substitute it for a joint multiclass model. citeturn355835search3turn355835search4

### Recommended newer module-centric option: `grpnet`

My third choice is **`grpnet`**, based on Helwig's group-elastic-net algorithms. The 2025 JCGS paper develops an adaptively bounded gradient method and explicitly includes high-dimensional multinomial examples with grouped predictors. citeturn368561search4turn368561search5 CRAN version **1.2**, published 29 April 2026, supports multinomial regression, group LASSO/MCP/SCAD, and an optional group-ridge penalty. citeturn368561search1turn368561search6

Its appeal for your priorities is unusually direct: `group` is a length-\(p\) vector defining predictor modules, and `alpha` interpolates between **group lasso and group ridge**. citeturn368561search0 It also returns class-specific multinomial coefficients and groupwise prediction/importance summaries. citeturn368561search2turn368561search3

I rank it below `msgl` for your primary analysis because `grpnet`'s native sparsity is fundamentally **between groups**: if a biological module is selected, it does not give the same explicit L1 within-module sparsity that `msgl` provides. Its main advantages are clean group-elastic-net semantics, support for group ridge, good coefficient access and a modern implementation. Its main disadvantages are its much newer ecosystem and less accumulated usage in omics.

The shortlist therefore looks like this:

| Rank | Method | Why it ranks here | Best role |
|---|---|---|---|
| **1** | `glmnet::glmnet(family="multinomial")` | most mature; elastic net; class coefficients; sparse matrices; fast paths; excellent baseline | principal joint classifier |
| **2** | `msgl` | true multinomial **sparse-group** lasso; module + within-module sparsity | biologically structured comparator |
| **3** | `grpnet` | native multinomial **group elastic net**, including group ridge | module-centric sensitivity analysis |

I would not prioritise `gglasso`, `SGL`, `sparsegl`, `grpreg` or `grpregOverlap` for the main joint multiclass benchmark. `gglasso` documents least-squares/binary-logistic/SVM group lasso rather than multinomial classification; `SGL` documents linear, binary logistic and Cox; standard `grpreg` documents grouped Gaussian/binomial/Poisson rather than multinomial; and `grpregOverlap` does not supply joint multinomial fitting. citeturn249502search3turn355835search6turn249502search2turn355835search3 `ordinalNet`-type models would likewise only be scientifically appropriate if your cell types had a genuine ordered outcome scale.

## Tree-based comparator

For this family I would deliberately keep things simple and use **`tidymodels` + `ranger`**, with XGBoost only as an optional secondary engine. The current `tidymodels` release is 1.5.0; `parsnip` directly supports `rand_forest(engine="ranger")` and `boost_tree(engine="xgboost")`, including class probabilities and tuning workflows. citeturn281258search12turn479698search0turn479698search3

This maps closely onto the supplied CardamomOT implementation. At the pinned commit, CardamomOT constructs `RandomForestClassifier(n_estimators=..., random_state=seed, n_jobs=-1)`; its default documentation/code uses 100 trees, and it assigns each training cell a weight proportional to \(n_c^{-\gamma}\), normalised to mean one, with \(\gamma=0.5\) as the documented default. Its optional logistic branch performs `log1p`, scaling and L2 logistic CV using balanced accuracy. citeturn613510view3

A deliberately boring R analogue is therefore:

```r
library(tidymodels)

rf_spec <-
  rand_forest(
    trees = 100,
    mtry  = tune(),
    min_n = tune()
  ) |>
  set_mode("classification") |>
  set_engine(
    "ranger",
    importance = "permutation",
    class.weights = class_weights
  )
```

`ranger`'s `class.weights` enter its cost-sensitive classification splitting and terminal voting; `case.weights`, by contrast, affect the probability that observations enter bootstrap/subsample samples. citeturn479698search10 Thus `class.weights ∝ n_c^-0.5` is the closer uncomplicated analogue of the CardamomOT intention. Because class frequencies change between outer training folds, calculate those weights from each training fold rather than from the complete dataset.

I prefer `tidymodels` here over `caret` because the former gives a current unified workflow for recipes, tuning, probability predictions and resampling; `caret` remains supported but is now the older orchestration layer. `mlr3` is also technically excellent and exposes `ranger`'s `class.weights` and permutation importance directly, so it would be a perfectly reasonable choice if your existing codebase is already `mlr3`-based. citeturn479698search4turn479698search7 There is no scientific gain in benchmarking all three orchestration APIs.

For XGBoost, simply substitute a `boost_tree()` specification and tune tree depth, number of trees and learning rate; `parsnip` currently provides a direct `xgboost` classification engine. citeturn479698search3 I would nevertheless designate RF, not XGBoost, as the prespecified family-B comparator because it reproduces your CardamomOT reference most closely and keeps the benchmark focused on feature-selection families rather than hyperparameter engineering.

Treat forest gene importance cautiously. Classical RF importance measures can be biased, and predictor correlation changes permutation-based importance behaviour; this is particularly relevant when hundreds of genes occur in tightly correlated modules. citeturn336157search0turn336157search3 Therefore use the forest first as a **predictive comparator**. For marker ranking, prefer held-out permutation importance and report selection/rank stability across outer resamples rather than interpreting the top ten genes as if they were independent biological effects.

## Sparse PLS-DA, integration and time

### Baseline sparse PLS-DA

Use **`mixOmics::splsda`** as the primary family-C model. The current Bioconductor release is `mixOmics` **6.36.0**. citeturn281258search0turn281258search2 `splsda()` accepts a categorical outcome, internally creates its indicator matrix, and selects variables by specifying the number of X variables retained on each component through `keepX`. citeturn529180search0 The mixOmics documentation explicitly describes sPLS-DA as adding lasso-type penalisation to the loading vectors. citeturn529180search2

That means sPLS-DA meets your **sparsity** requirement, but not literally your ridge requirement. `keepX`/sparse loadings are an L1-like selection mechanism; ordinary `mixOmics::splsda` does not expose a glmnet-style \(\alpha\) mixing L1 and ridge. I would not customise ridge into the first benchmark: doing so would make the sPLS-DA arm substantially less standard than the other methods.

Tune `ncomp`, `keepX` and the prediction distance inside the training data; mixOmics provides cross-validation tooling specifically for these parameters. citeturn529180search1turn529180search3

For RNA counts, my preferred input is **training-fold-fitted sctransform Pearson residuals**, followed by sPLS-DA, rather than pretending that PLS-DA itself has a negative-binomial observation model. `splsda` does offer a CLR option for compositional inputs, but the documentation notes that zeros must first be addressed. citeturn529180search0 SCT residuals give a cleaner separation between the count model and the discriminant model.

### P-integration and SBR/SB1

`mint.splsda()` is genuinely relevant, but its semantics should not be stretched. MINT was designed for **P-integration**: several studies measuring the same features, with systematic study-specific variation, while identifying a common discriminant signature. citeturn723925search4turn529180search5 The function takes a `study` factor for every observation and performs study-specific scaling; current documentation advises using studies with more than three observations and with outcome categories represented in each study. citeturn723925search6

Therefore, if SBR and SB1 are genuinely two cell lines/studies and **cell-line identity is known for every cell included in the analysis**, `mint.splsda(X, Y, study = cell_line, ...)` is a very sensible secondary model. It asks for a signature that generalises across those lines rather than letting line-specific shifts dominate.

If SBR/SB1 status is absent at some time points, however, MINT cannot magically treat “unknown study” as an ordinary missing predictor. I would either restrict the MINT sensitivity analysis to observations for which the study label is defined, or use the fully observed experimental unit/donor/batch as the study variable if that better represents the design. Do **not** impute SBR/SB1 merely to satisfy the algorithm.

Similarly, mixOmics' `multilevel` option models within-individual repeated-measure structure. citeturn529180search0turn529180search7 Destructive single-cell sampling means the cells themselves are not longitudinal repeated measures. If the same donor/culture is sampled at several times, then the repeated unit is the **donor/culture**, not a cell lineage inferred retrospectively from expression.

### Time and splines

I would not add a spline-enhanced sPLS-DA to the primary four-way benchmark. The Straube et al. paper you linked is a **linear mixed-model spline framework for repeated time-course omics profiles**, designed for modelling trajectories, clustering and differential expression; it is not a native PLS-DA extension. citeturn253719search0turn253719search3

It is nevertheless useful inspiration for a *secondary* analysis. There are three defensible interpretations of time:

1. **Time is nuisance/confounding.** Within each training fold, regress each gene on a spline basis of time (and possibly donor/cell line), freeze that model, and feed residuals to the classifier. This asks for cell-type information beyond predictable temporal shifts.

2. **Time is useful predictive context.** Add a low-dimensional spline basis \(B_1(t),\ldots,B_q(t)\) as extra predictors. For multinomial models this is easier to control than for mixOmics because those terms can be left unpenalised or weakly penalised.

3. **Generalisation across time is the scientific question.** Use leave-timepoint-out or forward-time evaluation rather than putting increasingly elaborate time terms into the classifier.

The third is often the most informative. Because scRNA-seq sampling is destructive, a cell observed at \(t_1\) and a cell at \(t_2\) are not the same longitudinal individual; claims should therefore be about **population-level temporal changes**, not observed single-cell trajectories. The LMMS paper itself assumes repeated molecular profiles over subjects/time and explicitly models subject-specific variability, illustrating why its interpretation does not transfer literally to unrelated cells. citeturn253719search0

Missing cell types at some time points create a more fundamental problem: if a held-out time contains a type that never occurred in training, a closed-set multinomial/RF/sPLS-DA classifier cannot learn that class. For time-generalisation analyses, either evaluate only types represented in training and test or explicitly label this as an open-set/unseen-class problem rather than scoring the missing class as an ordinary classification failure.

## Implementation blueprint and decision criteria

I would implement the benchmark as a **nested, group-aware outer CV** in which one frozen split object is reused everywhere. Inside every outer training fold: derive the gene universe from training cells only; fit SCT or other normalisation on training data; run your omnibus plus contrasts; tune `glmnet` \((\alpha,\lambda)\), `msgl` \((\alpha,\lambda)\), RF hyperparameters and sPLS-DA \((ncomp,keepX)\); estimate class weights; select the gene panel; refit; then transform and predict the untouched outer test donor(s). Donor-wise evaluation is the appropriate design when your target claim is generalisation to unseen biological replicates. citeturn902516search0turn902516search1

For the main paper or benchmark figure, I would prespecify four rows: **two-stage omnibus/contrast**, **`glmnet` multinomial elastic net**, **`ranger` RF**, and **`mixOmics::splsda`**. Then make `msgl` the structured extension of the multinomial arm rather than a fifth unrelated family. This keeps the scientific comparison easy to understand while testing the main hypothesis you care about: whether your marginal-but-interpretable two-stage selection performs as well as joint sparse discrimination, nonlinear ensembles and latent sparse discrimination.

A compact R skeleton is:

```r
# Outer split: whole donors / biological replicates
outer_folds <- rsample::group_vfold_cv(dat, group = donor_id)

# A. Main multinomial comparator
fit_en <- glmnet::cv.glmnet(
  x = X_train,
  y = y_train,
  family = "multinomial",
  alpha = alpha,                       # tune, e.g. 0.05 ... 1
  type.multinomial = "ungrouped",
  weights = cell_weights,
  foldid = inner_group_fold
)

# A. Structured module sensitivity analysis
fit_sgl <- msgl::cv(
  X_train,
  y_train,
  grouping = gene_module,
  alpha = alpha_sgl
)

# B. Deliberately simple CardamomOT-like RF
rf_spec <- rand_forest(
  trees = 100,
  mtry  = tune(),
  min_n = tune()
) |>
  set_engine(
    "ranger",
    class.weights = class_weights,
    importance = "permutation"
  ) |>
  set_mode("classification")

# C. Sparse PLS-DA
fit_splsda <- mixOmics::splsda(
  X_train,
  y_train,
  ncomp = ncomp,
  keepX = keepX
)
```

The `glmnet` syntax and multinomial modes above are directly supported in the current documentation. citeturn249502search8turn281258search19 `msgl` accepts arbitrary non-overlapping feature grouping for multinomial sparse-group lasso. citeturn355835search1 `ranger` is exposed through `parsnip` and supports class probability prediction and weighting. citeturn479698search0turn479698search10 `splsda`'s `keepX` controls sparse feature retention on each component. citeturn529180search0

For reporting, I would separate **prediction**, **selection** and **interpretation** rather than collapse them into one winner. Prediction should report donor-held-out macro-F1/balanced accuracy, per-class recall and calibration. Selection should report gene count, stability frequency across outer folds and stability at matched panel sizes. Interpretation should compare whether genes recovered for each cell type have coherent coefficient signs/loadings and whether selected genes cluster into the prior modules. Your two-stage method has a substantial advantage in the last category because the type-versus-grand-mean contrast has a direct marker interpretation; penalised multinomial regression has the advantage that redundancy and combinations of predictors are optimised jointly.

The most informative primary comparison is therefore **your two-stage markers versus `glmnet` ungrouped multinomial elastic net**. If the two-stage method matches its predictive performance with substantially fewer and more stable cell-type-specific markers, that is strong evidence in favour of your procedure. If `glmnet` wins substantially, `msgl` then tells you whether the gain comes from exploiting correlated/module-level structure. RF tells you whether nonlinearities/interactions account for additional gains, while sPLS-DA asks whether a low-rank multivariate signature captures discrimination more efficiently than either marginal screening or direct logistic coefficients. This hierarchy closely matches the sparse-classifier rationale already outlined in the DeCovarT feature-screening proposal, while making the evaluation donor-aware and leakage-safe. citeturn613510view0turn613510view1