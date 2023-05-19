### Logs
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

### Libraries
library(Seurat)
library(tidyverse)
library(scDblFinder)
library(BiocParallel)
library(SingleCellExperiment)

#-----------------------------------------------------------------------------------------------------------------------
# Remove cell doublets
#-----------------------------------------------------------------------------------------------------------------------
remove_doublets <- function(seurat, cores = 1) {
    sce <- NormalizeData(seurat) %>% as.SingleCellExperiment()
    sce <- scDblFinder(sce, samples = "cellhashing", BPPARAM = MulticoreParam(cores))

    dblt_info <- sce$scDblFinder.class
    names(dblt_info) <- rownames(sce@colData)
    seurat$scDblFinder.class <- dblt_info

    print("Total number of singlets and doublets")
    print(table(seurat$scDblFinder.class))

    return(seurat)
}

#-----------------------------------------------------------------------------------------------------------------------
# Main: # Create seurat object & calculate doublets
#-----------------------------------------------------------------------------------------------------------------------
seurat_objects <- map(snakemake@input, \(x) readRDS(x))

if (length(seurat_objects) > 1) {
    seurat <- merge(seurat_objects[[1]], seurat_objects[2:length(seurat_objects)])
} else {
    seurat <- unlist(seurat_objects)
}

seurat <- remove_doublets(seurat, cores = snakemake@threads[[1]])

# Filter doublets and save object
seurat_clean <- subset(seurat, subset = scDblFinder.class == "singlet")

saveRDS(seurat_clean, snakemake@output[["no_doublets"]])
saveRDS(seurat, snakemake@output[["raw"]])