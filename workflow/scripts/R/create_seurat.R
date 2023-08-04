### Logs
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

### Libraries
library(Seurat)
library(Signac)
library(tidyverse)
library(SingleCellExperiment)
library(scDblFinder)
library(BiocParallel)
library(DropletUtils)

#-----------------------------------------------------------------------------------------------------------------------
# Create seurat object
#-----------------------------------------------------------------------------------------------------------------------
create_seurat <- function(input_files, sample_name, cellhash_names, min_cells_gene = 1,
                          is_larry = FALSE, is_cell_hashing = FALSE, UMI_cutoff = 0) {
  # Load files
  names(input_files) <- sample_name
  data <- Read10X(data.dir = input_files)
  
  # Create seurat objects
  if (is.list(data)) {

    seurat <- CreateSeuratObject(counts = data$`Gene Expression`, min.cells = min_cells_gene)
    
    if (is_cell_hashing) {
      seurat[["Cellhashing"]] <- CreateAssayObject(counts = data$`Antibody Capture`)
    }
    
    if (is_larry) {
      seurat[["Larry"]] <- CreateAssayObject(counts = data$`Custom`)
    }

  } else {

    seurat <- CreateSeuratObject(counts = data, min.cells = min_cells_gene)

  }
  
  # Fix sample names. cellhashing column is added to be consistent with the structure of seurat objects
  # coming from libraries in which cellhashing has been performed.
  # Also a subsample entry is created to be consistent with samples with cellhashing.
  seurat$sample      <- sample_name
  seurat$subsample   <- sample_name
  Idents(seurat)     <- seurat$sample
  
  # Remove cells with less than UMI threshold
  seurat <- subset(seurat, subset = nCount_RNA >= UMI_cutoff)
  
  return(seurat)

}


create_seurat_arc <- function(input_files, input_fragments, input_larry = NULL, sample_name, min_cells_gene = 1,
                          is_larry = FALSE, UMI_cutoff = 0) {
  # Load files
  names(input_files) <- sample_name
  arc <- Read10X(data.dir = input_files)
  
  if (is_larry) {

    names(input_larry) <- sample_name
    larry <- Read10X(data.dir = input_larry)

    common_cells <- intersect(
      colnames(larry$Custom),
      colnames(arc$`Gene Expression`)
    )

    seurat            <- CreateSeuratObject(counts = arc$`Gene Expression`[,common_cells], min.cells = min_cells_gene)
    seurat[["ATAC"]]  <- CreateChromatinAssay(counts = arc$`Peaks`[,common_cells], fragments = input_fragments, sep = c(":", "-"))
    seurat[["Larry"]] <- CreateAssayObject(counts = larry$Custom[,common_cells])

  } else {

    seurat           <- CreateSeuratObject(counts = arc$`Gene Expression`, min.cells = min_cells_gene)
    seurat[["ATAC"]] <- CreateChromatinAssay(counts = arc$`Peaks`, fragments = input_fragments, sep = c(":", "-"))

  }
  
  # Fix sample names. cellhashing column is added to be consistent with the structure of seurat objects
  # coming from libraries in which cellhashing has been performed.
  # Also a subsample entry is created to be consistent with samples with cellhashing.
  seurat$sample      <- sample_name
  seurat$subsample   <- sample_name
  Idents(seurat)     <- seurat$sample
  
  # Remove cells with less than UMI threshold
  seurat <- subset(seurat, subset = nCount_RNA >= UMI_cutoff)
  
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
# Cellhashing assignment functions
#-----------------------------------------------------------------------------------------------------------------------
cell_hashing_assignment <- function(seurat, cellhash_names, sample_name) {
  
  sce               <- as.SingleCellExperiment(seurat)
  hash.stats        <- hashedDrops( counts(altExp(sce, "Cellhashing")), constant.ambient=TRUE)
  cellhash_ab_names <- metadata(hash.stats)$ambient %>% names()
  
  cell_hash_assignment <- hash.stats %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("cell_id") %>% 
    dplyr::mutate(subsample = cellhash_names[cellhash_ab_names[Best]]) %>% 
    dplyr::mutate(subsample = ifelse(Confident, subsample, paste0("Unassigned_", sample_name))) %>% 
    dplyr::select(cell_id, subsample) %>% 
    tibble::deframe() %>% 
    unlist()
  
  seurat$subsample <- cell_hash_assignment
  
  return(seurat)
  
}

cell_hashing_assignment_seurat <- function(seurat, cellhash_names, sample_name) {
  
  seurat <- NormalizeData(seurat, assay = "Cellhashing", normalization.method = "CLR")
  seurat <- HTODemux(seurat, assay = "Cellhashing", positive.quantile = 0.99)
  
  return(seurat)
  
}

summary_cellhashing <- function(seurat, cellhash_names) {
  
  # Plots based on seurat and dropletutils cell hashing assignment
  Idents(seurat) <- "Cellhashing_classification"
  p1 <- RidgePlot(seurat, assay = "Cellhashing", features = rownames(seurat[["Cellhashing"]]), ncol = 1)
  p2 <- VlnPlot(seurat, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
  Idents(seurat) <- "subsample"
  p3 <- RidgePlot(seurat, assay = "Cellhashing", features = rownames(seurat[["Cellhashing"]]), ncol = 1)
  p4 <- VlnPlot(seurat, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
  
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  
  # scatter plot of all pairwise combinations of Abs used for cellhashing.
  # First with emprydrops annotation, then with seurat annotation
  ab_combs <- snakemake@params[["cellhash_names"]] %>% names() %>% combn(2)
  
  for (i in ncol(ab_combs)) {
    comb <- ab_combs[,i]
    p <- FeatureScatter(seurat, feature1 = comb[1], feature2 = comb[2], raster = TRUE)
    print(p)
  }
  
  Idents(seurat) <- "Cellhashing_classification"
  for (i in ncol(ab_combs)) {
    comb <- ab_combs[,i]
    p <- FeatureScatter(seurat, feature1 = comb[1], feature2 = comb[2], raster = TRUE)
    print(p)
  }
  
}

#-----------------------------------------------------------------------------------------------------------------------
# Main & save output
#-----------------------------------------------------------------------------------------------------------------------
# Create seurat object & calculate doublets
if (snakemake@params[["library_type"]] == "ARC") {

  seurat <- create_seurat_arc(
    input_files     = snakemake@input[["arc"]],
    input_fragments = snakemake@input[["fragments"]],
    input_larry     = snakemake@input[["counts"]],
    sample_name     = snakemake@wildcards[["sample"]],
    min_cells_gene  = snakemake@params[["min_cells_gene"]],
    is_larry        = snakemake@params[["is_larry"]],
    UMI_cutoff      = snakemake@params[["umi_cutoff"]]
  )

} else {

  seurat <- create_seurat(
    input_files     = snakemake@input[["counts"]],
    sample_name     = snakemake@wildcards[["sample"]],
    min_cells_gene  = snakemake@params[["min_cells_gene"]],
    is_larry        = snakemake@params[["is_larry"]],
    is_cell_hashing = snakemake@params[["is_cell_hashing"]],
    cellhash_names  = snakemake@params[["cellhash_names"]],
    UMI_cutoff      = snakemake@params[["umi_cutoff"]]
  )

}

# Remove doublets
seurat <- remove_doublets(seurat, cores = snakemake@threads[[1]])

# Add mitochondrial and ribosomal RNA metrics
mito_pattern             <- snakemake@params[["mito_pattern"]]
ribo_pattern             <- snakemake@params[["ribo_pattern"]]
seurat[["percent.mt"]]   <- PercentageFeatureSet(seurat, pattern = mito_pattern)
seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = ribo_pattern)

# Filter doublets and save object
seurat_clean <- subset(seurat, subset = scDblFinder.class == "singlet")

# If the library has been processed with cell hashing, assign cells to subsamples in seurat without doublets
if (snakemake@params[["is_cell_hashing"]]) {
  
  seurat_clean <- cell_hashing_assignment(
    seurat         = seurat_clean,
    cellhash_names = snakemake@params[["cellhash_names"]],
    sample_name    = snakemake@wildcards[["sample"]]
  )
  
  seurat_clean <- cell_hashing_assignment_seurat(
    seurat         = seurat_clean,
    cellhash_names = snakemake@params[["cellhash_names"]],
    sample_name    = snakemake@wildcards[["sample"]]
  )
  
  pdf(paste0("results/02_createSeurat/", snakemake@wildcards[["sample"]], "_cellhash.pdf"), width = 7.5, height = 5)
  summary_cellhashing(seurat_clean, snakemake@params[["cellhash_names"]])
  dev.off()
  
}

# Save to rds
saveRDS(seurat_clean, snakemake@output[["no_doublets"]])
saveRDS(seurat, snakemake@output[["raw"]])