save.image()
### Logs
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

### Libraries
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(scDblFinder)
library(BiocParallel)

#-----------------------------------------------------------------------------------------------------------------------
# Create seurat object
#-----------------------------------------------------------------------------------------------------------------------
create_seurat <- function(input_files, sample_name, cellhash_names, min_cells_gene = 1, min_cells_larry = 1, is_larry = FALSE, is_cell_hashing = FALSE) {
    # Load files
    names(input_files) <- sample_name
    data <- Read10X(data.dir = input_files)

    # Create seurat objects
    seurat <- CreateSeuratObject(counts = data$`Gene Expression`, min.cells = min_cells_gene)

    if (is_cell_hashing) {
        seurat[["Cellhashing"]] <- CreateAssayObject(counts = data$`Antibody Capture`, min.cells = min_cells_larry)
    }

    if (is_larry) {
        seurat[["Larry"]] <- CreateAssayObject(counts = data$`Custom`, min.cells = min_cells_larry)
    }

    # Fix sample names. cellhashing column is added to be consistent with the structure of seurat objects
    # coming from libraries in which cellhashing has been performed.
    # Also a subsample entry is created to be consistent with samples with cellhashing.
    seurat$sample      <- sample_name
    seurat$subsample   <- sample_name
    Idents(seurat)     <- seurat$sample

    return(seurat)
}

#-----------------------------------------------------------------------------------------------------------------------
# Remove cell doublets
#-----------------------------------------------------------------------------------------------------------------------
remove_doublets <- function(seurat, cores = 1) {
    sce <- NormalizeData(seurat) %>% as.SingleCellExperiment()
    sce <- scDblFinder(sce, samples = "sample", BPPARAM = MulticoreParam(cores))

    dblt_info <- sce$scDblFinder.class
    names(dblt_info) <- rownames(sce@colData)
    seurat$scDblFinder.class <- dblt_info

    print("Total number of singlets and doublets")
    print(table(seurat$scDblFinder.class))

    return(seurat)
}

#-----------------------------------------------------------------------------------------------------------------------
# Main & save output
#-----------------------------------------------------------------------------------------------------------------------
# Create seurat object & calculate doublets
seurat <- create_seurat(
    input_files     = snakemake@input[[1]],
    sample_name     = snakemake@wildcards[["sample"]],
    min_cells_gene  = snakemake@params[["min_cells_gene"]],
    min_cells_larry = snakemake@params[["min_cells_larry"]],
    is_larry        = snakemake@params[["is_larry"]],
    is_cell_hashing = snakemake@params[["is_cell_hashing"]],
    cellhash_names  = snakemake@params[["cellhash_names"]]
)

seurat <- remove_doublets(seurat, cores = snakemake@threads[[1]])

# Add mitochondrial and ribosomal RNA metrics
mito_pattern             <- snakemake@params[["mito_pattern"]]
ribo_pattern             <- snakemake@params[["ribo_pattern"]]
seurat[["percent.mt"]]   <- PercentageFeatureSet(seurat, pattern = mito_pattern)
seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = ribo_pattern)

# Filter doublets and save object
seurat_clean <- subset(seurat, subset = scDblFinder.class == "singlet")

saveRDS(seurat_clean, snakemake@output[["no_doublets"]])
saveRDS(seurat, snakemake@output[["raw"]])