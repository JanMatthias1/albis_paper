# Direct Splatter reference generation; output values are preserved.
library(splatter)
library(Matrix)
library(jsonlite)
settings <- fromJSON(Sys.getenv('SETTINGS'))
out <- Sys.getenv('REFERENCE')
if (dir.exists(out)) stop('Reference directory already exists')
dir.create(out, recursive=TRUE)
started <- proc.time()
params <- newSplatParams(batchCells=settings$reference_cells, nGenes=556,
    group.prob=rep(1/8,8), de.prob=0.15, de.facLoc=0.5, de.facScale=0.4,
    dropout.type='experiment', dropout.mid=0, seed=settings$seed)
reference <- splatSimulateGroups(params, verbose=FALSE)
generation_seconds <- unname((proc.time()-started)['elapsed'])
stopifnot(nrow(reference)==556, ncol(reference)==settings$reference_cells)
# Serialize native counts and labels. Column names are input-schema names, not new model outputs.
writeMM(as(counts(reference), 'CsparseMatrix'), file.path(out,'counts.mtx'))
writeLines(rownames(reference),file.path(out,'genes.tsv'))
write.table(data.frame(Cell=colnames(reference), Cell_type=as.character(colData(reference)$Group)),
    file.path(out,'cells.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write_json(list(status='ok',method='splatter',seed=settings$seed,n_genes=556,
    n_cells=settings$reference_cells,generation_seconds=generation_seconds,
    version=as.character(packageVersion('splatter'))),file.path(out,'measurement.json'),auto_unbox=TRUE,pretty=TRUE)
