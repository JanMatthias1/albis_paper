#!/usr/bin/env Rscript
# Figure 4A -- run nnSVG on flat files written by 3D_nnsvg.py: logcounts.mtx
# (genes x cells, log-normalized), gene_names.txt, cell_names.txt,
# locations.csv (cell_id, x, y). Writes nnsvg_results.csv (per-gene BRISC fit
# + LR-test columns, incl. pval/padj) into the same directory.
#
# nnSVG's own BRISC::BRISC_order hard-requires exactly 2D coordinates (see
# project memory, 2026-09-05) -- there is no 3D mode, unlike scBSP/SPARK-X
# which were run on the pooled 3D z-stack. This is why 3D_nnsvg.py restricts
# to a single representative slice before calling this script.
#
# Expected environment: the R bundled in sim_paper/env/stagate-pyg (nnSVG
# installed there 2026-09-05 -- see project memory for the zlib.h/Magick++
# install workarounds if this R library ever needs rebuilding).
#
# Reproducibility + observability (2026-09-05): nnSVG/BRISC is stochastic
# (AMMD ordering + per-gene BRISC optimisation), so back-to-back runs drifted
# (spot AUROC 0.638 -> 0.615). We now set.seed() and pin the BiocParallel
# RNGseed. We also pass verbose = FALSE and an explicit MulticoreParam with
# progressbar = TRUE + one task per gene: nnSVG's verbose per-gene BRISC
# blocks are captured by BiocParallel (never reach the SLURM log), so the
# progress bar is the only usable signal that a long run is advancing.
#
# Usage: Rscript run_nnsvg.R <indir> [n_threads] [seed]

suppressMessages({
  library(Matrix)
  library(nnSVG)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
n_threads <- if (length(args) >= 2) as.integer(args[2]) else 4
seed <- if (length(args) >= 3) as.integer(args[3]) else 42L

y <- as.matrix(readMM(file.path(indir, "logcounts.mtx")))  # genes x cells, dense
gene_names <- readLines(file.path(indir, "gene_names.txt"))
cell_names <- readLines(file.path(indir, "cell_names.txt"))
rownames(y) <- gene_names
colnames(y) <- cell_names

loc <- read.csv(file.path(indir, "locations.csv"), row.names = 1)
coords <- as.matrix(loc[cell_names, c("x", "y")])

cat(sprintf("[data] y: %d genes x %d cells; coords: %d x %d\n",
            nrow(y), ncol(y), nrow(coords), ncol(coords)))

# tasks: fixed chunk count (NOT nrow(y)) -- per-gene dispatch doubled spot
# wall (295s vs ~130s, 6.5min of fork/serialise sys time). 64 chunks gives
# ~9 genes/chunk: 64 progress-bar ticks for ETA, dispatch overhead amortised.
# Fixed absolute value + RNGseed => L'Ecuyer stream per chunk, so results are
# reproducible regardless of `workers`.
n_tasks <- min(nrow(y), 64L)
set.seed(seed)
bpparam <- MulticoreParam(workers = n_threads, RNGseed = seed,
                          progressbar = TRUE, tasks = n_tasks)
cat(sprintf("[nnsvg] start: workers=%d seed=%d tasks=%d  %s\n",
            n_threads, seed, n_tasks, format(Sys.time())))

t0 <- Sys.time()
res <- nnSVG(input = y, spatial_coords = coords, n_threads = n_threads,
             BPPARAM = bpparam, verbose = FALSE)
cat(sprintf("[nnsvg] done in %.1fs  %s\n",
            as.numeric(Sys.time() - t0, units = "secs"), format(Sys.time())))

out <- as.data.frame(res)
out$gene_names <- rownames(out)
write.csv(out, file.path(indir, "nnsvg_results.csv"), row.names = FALSE)
cat(sprintf("[save] %s\n", file.path(indir, "nnsvg_results.csv")))
