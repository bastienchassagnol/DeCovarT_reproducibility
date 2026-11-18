# denoised_covariance

This folder holds the first-pass plug-in script: the closed-form Schäfer-Strimmer shrinkage level, applied to the raw and to the denoised covariance. It was archived because the closed form minimises a Frobenius risk, not the Gaussian likelihood, and the denoised covariance gave the worst deconvolution (RMSE 0.47 at 96 h, Welch panel). Its functions are in `../sparse_precision_arms/R_03_01_stats_covariance_estimation.R`.

- `scripts/03_05_stats_plugin_covariance.R` (archived as `03_05_stats_plugin_covariance_closed_form.R`): plug-in within-type moments for the quick deconvolution benchmark. For each time, gene set (welch, r2_05) and scored type it writes a raw and a denoised covariance with the same closed-form intensity, to `output/grn/plugin/`. No cross-validation. The current script of the same name in scripts/ is a later version. archived at commit 64f1454 (working tree of 2026-11-16).
