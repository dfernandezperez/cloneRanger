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
# Create seurat object
#-----------------------------------------------------------------------------------------------------------------------
create_seurat <- function(input_files, sample_names, min.cells = 1) {

	# Load files
	names(input_files) <- sample_names
	data               <- Read10X(data.dir = input_files)

	# Create seurat objects
	seurat            <- CreateSeuratObject(counts = data$`Gene Expression`, min.cells = min.cells)
	seurat[['Larry']] <- CreateAssayObject(counts = data$`Custom`, min.cells = min.cells)

	# Fix sample names
	samples_names        <- seurat$orig.ident %>% names() %>% gsub("_[A,C,G,T]*-1", "", .)
	names(samples_names) <- seurat$orig.ident %>% names()
	seurat$sample        <- samples_names
	Idents(seurat)       <- seurat$sample

	return(seurat)

}

#-----------------------------------------------------------------------------------------------------------------------
# Remove cell doublets
#-----------------------------------------------------------------------------------------------------------------------
remove_doublets <- function(seurat, cores = 1) {

	sce <- NormalizeData(seurat) %>% as.SingleCellExperiment
	sce <- scDblFinder(sce, samples = "sample", BPPARAM = MulticoreParam(cores))

	dblt_info                <- sce$scDblFinder.class
	names(dblt_info)         <- rownames(sce@colData)
	seurat$scDblFinder.class <- dblt_info

	print("Total number of singlets and doublets")
	print(table(seurat$scDblFinder.class))

	return(seurat)

}

#-----------------------------------------------------------------------------------------------------------------------
# Main & save output
#-----------------------------------------------------------------------------------------------------------------------
# Create seurat object & calculate doublets
seurat <- create_seurat(unlist(snakemake@input), snakemake@params[["sample_names"]], snakemake@params[["min_cells"]])
seurat <- remove_doublets(seurat, cores = snakemake@threads[[1]])

# Filter doublets and save object
seurat_clean <- subset(seurat, subset = scDblFinder.class == "singlet")
saveRDS(seurat_clean, snakemake@output[["no_doublets"]])
saveRDS(seurat, snakemake@output[["raw"]])