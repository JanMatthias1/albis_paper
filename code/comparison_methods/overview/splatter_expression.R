#!/usr/bin/env Rscript
# Splatter expression pool for Figure 5 (5A, and the 5D SPIDER runs): the synthetic
# single-cell counts that SPIDER and scCube draw their expression from. A direct
# newSplatParams + splatSimulateGroups call with fixed parameters (no empirical data);
# counts are written unchanged as Matrix Market + TSV so the method environments can
# read them with numpy/scipy/pandas, without an R-Python bridge. Group1..GroupK labels
# are renamed positionally to type1..typeK (the cell-type names used in every output).
#
# Run with the `splatter` env's Rscript (see submit_figure5A.sh):
#   env/splatter/bin/Rscript splatter_expression.R \
#     --contract contract.json [--params params.json] --seed 20260921 --out-dir <dir>

suppressMessages({
  library(optparse)
  library(jsonlite)
  library(splatter)
  library(Matrix)
})

option_list <- list(
  make_option("--contract", type = "character"),
  make_option("--params", type = "character", default = NULL),
  make_option("--seed", type = "integer"),
  make_option("--out-dir", type = "character", dest = "out_dir")
)
opt <- parse_args(OptionParser(option_list = option_list))

contract <- fromJSON(opt$contract)
params <- if (!is.null(opt$params)) fromJSON(opt$params) else list()

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
manifest_path <- file.path(opt$out_dir, "manifest.json")

n_cells <- as.integer(contract$n_cells)
n_genes <- as.integer(contract$n_genes)
n_cell_types <- as.integer(contract$n_cell_types)
group_prob <- as.numeric(contract$cell_type_proportions)
group_prob <- group_prob / sum(group_prob)

de_prob <- if (!is.null(params$de_prob)) params$de_prob else 0.15
de_fac_loc <- if (!is.null(params$de_fac_loc)) params$de_fac_loc else 0.5
de_fac_scale <- if (!is.null(params$de_fac_scale)) params$de_fac_scale else 0.4
dropout_type <- if (!is.null(params$dropout_type)) params$dropout_type else "experiment"
dropout_mid <- if (!is.null(params$dropout_mid)) params$dropout_mid else 0.0

t0 <- Sys.time()
status <- "ok"
failure_reason <- NA

result <- tryCatch({
  splat_params <- newSplatParams(
    batchCells = n_cells,
    nGenes = n_genes,
    group.prob = group_prob,
    de.prob = de_prob,
    de.facLoc = de_fac_loc,
    de.facScale = de_fac_scale,
    dropout.type = dropout_type,
    dropout.mid = dropout_mid,
    seed = as.integer(opt$seed)
  )
  sim <- splatSimulateGroups(splat_params, verbose = FALSE)
  list(sim = sim, splat_params = splat_params)
}, error = function(e) {
  status <<- "failed"
  failure_reason <<- paste0(class(e)[1], ": ", conditionMessage(e))
  NULL
})

runtime_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
peak_memory_bytes <- as.numeric(sum(gc()[, 6])) * 1024 * 1024  # "max used" columns, Mb -> bytes

output_paths <- list(
  counts_mtx = file.path(opt$out_dir, "counts.mtx"),
  genes_tsv = file.path(opt$out_dir, "genes.tsv"),
  cells_tsv = file.path(opt$out_dir, "cells.tsv"),
  manifest = manifest_path
)

if (!is.null(result)) {
  sim <- result$sim
  counts <- as(counts(sim), "CsparseMatrix")  # genes x cells, as Splatter returns it
  writeMM(counts, output_paths$counts_mtx)
  writeLines(rownames(counts), output_paths$genes_tsv)

  group <- as.character(colData(sim)$Group)
  group_levels <- sort(unique(group))
  canonical <- paste0("type", seq_len(n_cell_types))
  # Splatter's Group1..GroupK order matches group_prob's input order, which
  # is the contract's cell_type_proportions order -- map positionally to the
  # shared "type1".."typeK" convention used across every workflow's output.
  remap <- setNames(canonical[seq_along(group_levels)], group_levels)
  cell_type_true <- unname(remap[group])

  cells_df <- data.frame(
    cell_id = colnames(counts),
    cell_type_true = cell_type_true,
    stringsAsFactors = FALSE
  )
  write.table(cells_df, output_paths$cells_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  file.create(output_paths$counts_mtx)
  file.create(output_paths$genes_tsv)
  file.create(output_paths$cells_tsv)
}

manifest <- list(
  workflow = "splatter_expression",
  component_versions = list(
    language = "R",
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    splatter = as.character(packageVersion("splatter"))
  ),
  params = list(
    n_cells = n_cells, n_genes = n_genes, n_cell_types = n_cell_types,
    group_prob = group_prob, de_prob = de_prob, de_fac_loc = de_fac_loc,
    de_fac_scale = de_fac_scale, dropout_type = dropout_type, dropout_mid = dropout_mid
  ),
  seed = opt$seed,
  input_provenance = list(empirical_data_required = FALSE, pretrained_weights_required = FALSE),
  output_paths = output_paths,
  requested_targets = NULL,
  achieved_targets = NULL,
  runtime_seconds = runtime_seconds,
  peak_memory_bytes = peak_memory_bytes,
  status = status,
  failure_reason = failure_reason,
  host = Sys.info()[["nodename"]]
)
write(toJSON(manifest, auto_unbox = TRUE, null = "null", digits = NA), manifest_path)

cat(toJSON(list(status = status, out_dir = opt$out_dir), auto_unbox = TRUE), "\n")
