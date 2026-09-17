#!/usr/bin/env Rscript
# Figure 4A follow-up -- run SPARK-X's sparkx() on flat files written by
# domain_mix_sweep_generate.py (or 3D_sparkx-style exports): counts.mtx
# (genes x cells), gene_names.txt, cell_names.txt, locations.csv (cell_id,
# x, y[, z]). Writes sparkx_pvalues.csv (combinedPval, adjustedPval,
# gene_names) into the same directory.
#
# Expected environment: the R bundled in sim_paper/env/stagate-pyg (SPARK
# installed there 2026-09-05 -- see project memory for the Matrix/C++14
# install workarounds if this R library ever needs rebuilding).
#
# Usage: Rscript run_sparkx.R <indir> [--2d]
#   --2d uses locations.csv's x,y columns only (drops z) -- for a single
#   2D-slice export where z is otherwise a constant column that would make
#   SPARK-X's projection-kernel design matrix singular.

suppressMessages({
  library(Matrix)
  library(SPARK)
})

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
use_2d <- "--2d" %in% args

counts <- as(readMM(file.path(indir, "counts.mtx")), "CsparseMatrix")
gene_names <- readLines(file.path(indir, "gene_names.txt"))
cell_names <- readLines(file.path(indir, "cell_names.txt"))
rownames(counts) <- gene_names
colnames(counts) <- cell_names

loc_cols <- if (use_2d) c("x", "y") else c("x", "y", "z")
loc <- read.csv(file.path(indir, "locations.csv"), row.names = 1)
loc <- as.matrix(loc[cell_names, loc_cols])

cat(sprintf("[data] counts: %d genes x %d cells; loc: %d x %d\n",
            nrow(counts), ncol(counts), nrow(loc), ncol(loc)))

t0 <- Sys.time()
res <- sparkx(counts, loc, numCores = 4, option = "mixture")
cat(sprintf("[sparkx] done in %.1fs\n", as.numeric(Sys.time() - t0, units = "secs")))

out <- res$res_mtest
out$gene_names <- rownames(out)
write.csv(out, file.path(indir, "sparkx_pvalues.csv"), row.names = FALSE)
cat(sprintf("[save] %s\n", file.path(indir, "sparkx_pvalues.csv")))
