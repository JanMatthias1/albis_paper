#!/usr/bin/env Rscript
# RCTD (Robust Cell Type Decomposition, https://github.com/dmcable/spacexr)
# spatial deconvolution of Figure 2's spot data.
#
# Rationale: a Visium spot mixes several cells, so asking a spot-level
# clustering method to assign one discrete cell type (as the Figure 3
# cell-typing panels do for the `cell` modality) is the wrong task for
# `spot` -- recovering per-spot cell-type MIXTURE PROPORTIONS is. This
# simulation is a natural fit for that: `cell` (Figure 2's ground-truth-
# labeled single cells, SAME underlying tissue/gene panel/cell-type profiles
# as `spot`) serves directly as RCTD's reference, and `spot`'s
# obsm['cell_type_frac_true'] gives an EXACT per-spot ground truth to score
# against -- no need to fabricate either side of the evaluation.
#
# Column order of cell_type_frac_true verified empirically before writing
# this script (2026-09-18): columns are type1..type8 in that exact natural
# order -- argmax(cell_type_frac_true) matches obs['cell_type_true'] for
# 10034/10034 spots, 0 mismatches.
#
# Mode: "full" (arbitrary per-spot mixtures), not RCTD's "doublet" mode
# (assumes ~1-2 cell types per spot) -- the domain structure here creates
# realistic multi-cell-type mixing per spot, so full mode is the right
# comparison. Change DOUBLET_MODE below to switch.
#
# Reference and query both use the exact same 556-gene ALBIS panel by
# construction, so no gene-panel intersection step is needed.
#
# Input: Figure 2's smaller_sphere, locked 2026-09-17 (log_mu=-2.5 cell retune):
#   reference (cell): data/figure_2/smaller_sphere/data/
#                      log_mu_-2.5_theta_0.40_jitter0.15_bsigma15/simulation_cell_z_qc.h5ad
#   query (spot):      data/figure_2/smaller_sphere/data/
#                      packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03/simulation_spot_z_qc.h5ad
#
# Output, under data/figure_3/spatial_deconvolution/RCTD/spot/:
#   rctd_results.rds        the full myRCTD object (weights, singlet scores, etc.)
#   per_spot_metrics.csv    per-spot Pearson correlation + RMSE (est vs. true fractions)
#   per_celltype_metrics.csv  per-cell-type correlation across all spots
#   metrics_summary.json    overall correlation/RMSE, printed to stdout too
#   est_vs_true_scatter.png  diagnostic scatter, estimated vs. true fraction, by cell type

suppressPackageStartupMessages({
  library(spacexr)
  library(reticulate)
  library(Matrix)
  library(ggplot2)
})

CELL_H5AD <- "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_2/smaller_sphere/data/log_mu_-2.5_theta_0.40_jitter0.15_bsigma15/simulation_cell_z_qc.h5ad"
SPOT_H5AD <- "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_2/smaller_sphere/data/packing_pf0p04_log_mu_-2.0_theta_0.25_jitter0.10_bsigma03/simulation_spot_z_qc.h5ad"
OUT_DIR <- "/dcs04/hicks/data/Jan/sim_project/sim_paper/data/figure_3/spatial_deconvolution/RCTD/spot"
DOUBLET_MODE <- "full"
N_CELL_TYPES <- 8
MAX_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# The R "anndata" package's default (convert=TRUE) reticulate conversion of
# a scipy sparse matrix produces an R sparseMatrix whose 'x' slot is still
# the h5ad's original int32 dtype -- Matrix's validObject() rejects that the
# moment ANYTHING touches the object ("invalid class dgRMatrix object: 'x'
# slot is not of type double"), even just accessing $X itself. Bypassing
# that entirely: import anndata/numpy with convert=FALSE so nothing
# auto-converts, force the sparse data to float64 on the PYTHON side (where
# it's cheap and unambiguous) before it ever crosses into R, then build the
# dgCMatrix by hand from its raw CSC components.
py_ad <- reticulate::import("anndata", convert = FALSE)
py_np <- reticulate::import("numpy", convert = FALSE)

load_h5ad_counts_genes_by_obs <- function(path) {
  a <- py_ad$read_h5ad(path)
  Xcsc <- a$X$tocsc()
  Xcsc$data <- Xcsc$data$astype(py_np$float64)
  dims <- unlist(reticulate::py_to_r(Xcsc$shape))  # (n_obs, n_vars)
  m_obs_by_gene <- new(
    "dgCMatrix",
    i = as.integer(reticulate::py_to_r(Xcsc$indices)),
    p = as.integer(reticulate::py_to_r(Xcsc$indptr)),
    x = as.double(reticulate::py_to_r(Xcsc$data)),
    Dim = as.integer(dims)
  )
  list(
    counts = Matrix::t(m_obs_by_gene),  # genes x obs, what RCTD wants
    var_names = reticulate::py_to_r(a$var_names$tolist()),
    n_obs = as.integer(dims[1]),
    py_obj = a  # kept for pulling obs/obsm fields below
  )
}

get_obs_column <- function(py_obj, col) {
  reticulate::py_to_r(py_obj$obs[col]$astype("str")$values$tolist())
}

get_obsm <- function(py_obj, key) {
  reticulate::py_to_r(py_obj$obsm[key])
}

cat("[load] reference (cell):", CELL_H5AD, "\n")
cell_h5 <- load_h5ad_counts_genes_by_obs(CELL_H5AD)
ref_counts <- cell_h5$counts
rownames(ref_counts) <- cell_h5$var_names
ref_ids <- paste0("ref_", make.unique(as.character(seq_len(ncol(ref_counts)))))
colnames(ref_counts) <- ref_ids

cell_type_levels <- paste0("type", seq_len(N_CELL_TYPES))
cell_types <- factor(get_obs_column(cell_h5$py_obj, "cell_type_true"), levels = cell_type_levels)
names(cell_types) <- ref_ids
stopifnot(!any(is.na(cell_types)))

nUMI_ref <- Matrix::colSums(ref_counts)
names(nUMI_ref) <- ref_ids

cat("[reference] n_cells=", ncol(ref_counts), " n_genes=", nrow(ref_counts),
    " per-type counts:\n", sep = "")
print(table(cell_types))

reference <- Reference(ref_counts, cell_types, nUMI_ref)

cat("\n[load] query (spot):", SPOT_H5AD, "\n")
spot_h5 <- load_h5ad_counts_genes_by_obs(SPOT_H5AD)
query_counts <- spot_h5$counts
rownames(query_counts) <- spot_h5$var_names
spot_ids <- paste0("spot_", make.unique(as.character(seq_len(ncol(query_counts)))))
colnames(query_counts) <- spot_ids

coords <- as.data.frame(get_obsm(spot_h5$py_obj, "spatial")[, 1:2])
colnames(coords) <- c("x", "y")
rownames(coords) <- spot_ids

nUMI_query <- Matrix::colSums(query_counts)
names(nUMI_query) <- spot_ids

cat("[query] n_spots=", ncol(query_counts), " n_genes=", nrow(query_counts), "\n", sep = "")

puck <- SpatialRNA(coords, query_counts, nUMI_query)

cat("\n[RCTD] create.RCTD (doublet_mode=", DOUBLET_MODE, ", max_cores=", MAX_CORES, ")\n", sep = "")
myRCTD <- create.RCTD(puck, reference, max_cores = MAX_CORES)
myRCTD <- run.RCTD(myRCTD, doublet_mode = DOUBLET_MODE)

saveRDS(myRCTD, file.path(OUT_DIR, "rctd_results.rds"))

# --- score against ground truth ------------------------------------------
weights <- as.matrix(normalize_weights(myRCTD@results$weights))
weights <- weights[, cell_type_levels, drop = FALSE]

true_frac <- as.matrix(get_obsm(spot_h5$py_obj, "cell_type_frac_true"))
colnames(true_frac) <- cell_type_levels
rownames(true_frac) <- spot_ids

common <- intersect(rownames(weights), rownames(true_frac))
cat("\n[score] ", length(common), " / ", ncol(query_counts),
    " spots have an RCTD estimate (rest failed QC/min-UMI inside RCTD)\n", sep = "")
est <- weights[common, , drop = FALSE]
truth <- true_frac[common, , drop = FALSE]

safe_cor <- function(a, b) if (sd(a) == 0 || sd(b) == 0) NA_real_ else cor(a, b)

per_spot_cor <- vapply(seq_len(nrow(est)), function(i) safe_cor(est[i, ], truth[i, ]), numeric(1))
per_spot_rmse <- sqrt(rowMeans((est - truth)^2))
per_spot_df <- data.frame(
  spot_id = common,
  pearson_r = per_spot_cor,
  rmse = per_spot_rmse
)
write.csv(per_spot_df, file.path(OUT_DIR, "per_spot_metrics.csv"), row.names = FALSE)

per_type_cor <- vapply(seq_len(ncol(est)), function(j) safe_cor(est[, j], truth[, j]), numeric(1))
per_type_rmse <- sqrt(colMeans((est - truth)^2))
per_type_df <- data.frame(
  cell_type = cell_type_levels,
  pearson_r = per_type_cor,
  rmse = per_type_rmse
)
write.csv(per_type_df, file.path(OUT_DIR, "per_celltype_metrics.csv"), row.names = FALSE)

overall_cor <- safe_cor(as.vector(est), as.vector(truth))
overall_rmse <- sqrt(mean((est - truth)^2))
summary_list <- list(
  n_spots_scored = length(common),
  n_spots_total = ncol(query_counts),
  doublet_mode = DOUBLET_MODE,
  overall_pearson_r = overall_cor,
  overall_rmse = overall_rmse,
  mean_per_spot_pearson_r = mean(per_spot_cor, na.rm = TRUE),
  mean_per_spot_rmse = mean(per_spot_rmse)
)
writeLines(jsonlite::toJSON(summary_list, auto_unbox = TRUE, pretty = TRUE),
           file.path(OUT_DIR, "metrics_summary.json"))

cat("\n=== RESULTS ===\n")
cat("spots scored:", length(common), "/", ncol(query_counts), "\n")
cat("overall Pearson r (all spots x all cell types, pooled):", round(overall_cor, 4), "\n")
cat("overall RMSE:", round(overall_rmse, 4), "\n")
cat("mean per-spot Pearson r:", round(mean(per_spot_cor, na.rm = TRUE), 4), "\n")
cat("per-cell-type Pearson r:\n")
print(per_type_df)

# --- diagnostic plot -------------------------------------------------------
plot_df <- data.frame(
  estimated = as.vector(est),
  true = as.vector(truth),
  cell_type = rep(cell_type_levels, each = nrow(est))
)
p <- ggplot(plot_df, aes(x = true, y = estimated)) +
  geom_point(alpha = 0.15, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~cell_type, nrow = 2) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "RCTD spot deconvolution: estimated vs. true cell-type fraction",
       x = "true fraction (cell_type_frac_true)", y = "RCTD estimated fraction") +
  theme_bw()
ggsave(file.path(OUT_DIR, "est_vs_true_scatter.png"), p, width = 9, height = 5, dpi = 150)

cat("\n[done] outputs written to", OUT_DIR, "\n")
