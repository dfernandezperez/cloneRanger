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
create_seurat_cellhashing <- function(input_files, cell_assign, sample_name, subsample_names, min_cells_gene = 1, min_cells_larry = 1, is_larry = FALSE) {
    # Load files
    cell_assign <- read_csv(cell_assign)
    names(input_files) <- sample_name
    data <- Read10X(data.dir = input_files)

    # Create seurat objects
    seurat <- CreateSeuratObject(counts = data$`Gene Expression`, min.cells = min_cells_gene)
    seurat[["Multiplexing Capture"]] <- CreateAssayObject(counts = data$`Multiplexing Capture`, min.cells = 1)

    if (is_larry) {
        seurat[["Larry"]] <- CreateAssayObject(counts = data$`Custom`, min.cells = min_cells_larry)
    }

    # Filter valid cell barcodes since we read a raw barcode matrix
    cell_assign <- cell_assign %>% # Remove multiplets
        filter(!Assignment == "Multiplet")

    seurat <- subset(seurat, cells = paste(sample_name, cell_assign$Barcode, sep = "_"))

    # Split the cells into cellhashing assigned cells with the proper subsample name.
    # Those unassigned will be named with the "general" sample name (from units.tsv) followed by _Unassigned.
    cell_assign <- cell_assign %>%
        mutate(sample = case_when(
            Assignment %in% names(subsample_names) ~ unlist(subsample_names)[Assignment],
            TRUE ~ paste(sample_name, "Unassigned", sep = "_")
        ))

    subsample_names <- cell_assign %>%
        dplyr::select(Barcode, sample) %>%
        mutate(Barcode = paste(sample_name, Barcode, sep = "_")) %>%
        deframe()

    seurat$sample  <- subsample_names
    Idents(seurat) <- seurat$sample
    # Add a column with the cellhashing sample name. Important for the next step in which doublets will be removed.
    seurat$cellhashing <- sample_name


    return(seurat)
}


create_seurat_single <- function(input_files, sample_name, min_cells_gene = 1, min_cells_larry = 1, is_larry = FALSE) {
    # Load files
    names(input_files) <- sample_name
    data <- Read10X(data.dir = input_files)

    # Create seurat objects
    seurat <- CreateSeuratObject(counts = data$`Gene Expression`, min.cells = min_cells_gene)

    if (is_larry) {
        seurat[["Larry"]] <- CreateAssayObject(counts = data$`Custom`, min.cells = min_cells_larry)
    }

    # Fix sample names. cellhashing column is added to be consistent with the structure of seurat objects
    # coming from libraries in which cellhashing has been performed.
    seurat$sample      <- sample_name
    seurat$cellhashing <- sample_name
    Idents(seurat)     <- seurat$sample

    return(seurat)
}

create_seurat <- function(input_files, cell_assign, sample_name, subsample_names, min_cells_gene = 1, min_cells_larry = 1, is_larry = FALSE, is_cell_hashing = FALSE) {
    if (is_cell_hashing) {
        seurat <- create_seurat_cellhashing(
            snakemake@input[["cellranger_mtx"]],
            snakemake@params[["cell_assign"]],
            snakemake@wildcards[["sample"]],
            snakemake@params[["subsample_names"]],
            snakemake@params[["min_cells_gene"]],
            snakemake@params[["min_cells_larry"]],
            snakemake@params[["is_larry"]]
        )
    } else {
        seurat <- create_seurat_single(
            snakemake@params[["sample_path"]],
            snakemake@wildcards[["sample"]],
            snakemake@params[["min_cells_gene"]],
            snakemake@params[["min_cells_larry"]],
            snakemake@params[["is_larry"]]
        )
    }

    return(seurat)
}

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
# Main & save output
#-----------------------------------------------------------------------------------------------------------------------
# Create seurat object & calculate doublets
seurat <- create_seurat(
    snakemake@input[["cellranger_mtx"]],
    snakemake@params[["cell_assign"]],
    snakemake@wildcards[["sample"]],
    snakemake@params[["subsample_names"]],
    snakemake@params[["min_cells_gene"]],
    snakemake@params[["min_cells_larry"]],
    snakemake@params[["is_larry"]],
    snakemake@params[["is_cell_hashing"]]
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