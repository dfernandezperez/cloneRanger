### Logs
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

### Libraries
library(Seurat)
library(tidyverse)

### Read and merge seurat objects if there are multiple samples
if (length(snakemake@input)  == 1) {
    seurat <- readRDS(snakemake@input[[1]])
} else {
    seurat_objects <- map(snakemake@input, \(x) readRDS(x))
    seurat         <- merge(seurat_objects[[1]], seurat_objects[2:length(seurat_objects)] %>% unlist())
    seurat         <- JoinLayers(seurat)
}

saveRDS(seurat, snakemake@output[["seurat"]])
