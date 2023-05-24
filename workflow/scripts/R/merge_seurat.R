### Logs
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

### Libraries
library(Seurat)
library(tidyverse)

#-----------------------------------------------------------------------------------------------------------------------
# Main: # Create seurat object & calculate doublets
#-----------------------------------------------------------------------------------------------------------------------
seurat_objects <- map(snakemake@input, \(x) readRDS(x))

if (length(seurat_objects) > 1) {
    seurat <- merge(seurat_objects[[1]], seurat_objects[2:length(seurat_objects)])
} else {
    seurat <- unlist(seurat_objects)
}

saveRDS(seurat, snakemake@output[["seurat"]])