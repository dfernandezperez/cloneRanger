
### Log
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

### Libraries
library(Seurat)
library(tidyverse)
library(DropletUtils)
source("workflow/scripts/R/functions.R")

#-----------------------------------------------------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------------------------------------------------
get_larry_molinfo <- function(molinfo) {
  
  larry_bc_pos <- which(molinfo$feature.type == "Custom")

  larry_molinfo <- molinfo$data %>% 
    data.frame() %>% 
    dplyr::filter(gene %in% larry_bc_pos) %>% 
    mutate(larry_bc = molinfo$genes[gene]) %>% # Create col with larry bc name instead of idx from molinfo$genes
    mutate(larry_color = str_replace(larry_bc, "_.*$", "")) %>% # Add color information to split barcodes
    dplyr::select(-gene)
  
    return(larry_molinfo)
  
}


filter_umi_reads <- function(df_molinfo, read_threshold) {
  
  df_molinfo <- df_molinfo %>% 
    dplyr::filter(reads >= read_threshold)
    
  return(df_molinfo)
  
}


filter_bc_umis <- function(df_molinfo, umi_threshold) {
  
  df_molinfo <- df_molinfo %>% 
    dplyr::group_by(cell, larry_bc, larry_color) %>% 
    dplyr::count() %>% 
    dplyr::rename(n_umi = n) %>% 
    dplyr::filter(n_umi >= umi_threshold) %>% 
    dplyr::ungroup()
  
  return(df_molinfo)
  
}


#-----------------------------------------------------------------------------------------------------------------------
# Data laoding
#-----------------------------------------------------------------------------------------------------------------------
seurat_rds    <- snakemake@input[[1]]
molec_info_h5 <- snakemake@params[["molecule_info"]]
output_rds    <- snakemake@output[[1]]
  
umi_cutoff  <- snakemake@params[["umi_cutoff"]]
read_cutoff <- snakemake@params[["reads_cutoff"]]
  
seurat      <- readRDS(seurat_rds)
molinfo     <- read10xMolInfo(molec_info_h5,
                              get.cell    = TRUE,
                              get.umi     = TRUE,
                              get.gem     = FALSE,
                              get.gene    = TRUE,
                              get.reads   = TRUE,
                              get.library = FALSE
                              )

# Modify cell names form molinfo file to match those from seurat (basically add whatever is)
# before and after the cellular barcode (which consists in random 16 nucleotides, by now)
cell_prefix <- Cells(seurat) %>% str_replace("[AGTC]{16}-.*", "")
cell_suffix <- Cells(seurat) %>% str_replace(".*_[AGTC]{16}", "")

molinfo$data <- molinfo$data %>% 
  data.frame() %>%
  dplyr::mutate(cell = paste0(cell_prefix, cell, cell_suffix)) %>% 
  filter(cell %in% Cells(seurat)) %>% 
  DataFrame()

molinfo_larry <- get_larry_molinfo(molinfo)


#-----------------------------------------------------------------------------------------------------------------------
# Barcode calling
#-----------------------------------------------------------------------------------------------------------------------
larry_filt <- molinfo_larry %>% 
  filter_umi_reads(read_cutoff) %>% 
  filter_bc_umis(umi_cutoff)
  
larry_bc_calls <- larry_filt %>% 
  dplyr::group_by(cell, larry_color) %>% 
  dplyr::filter(n_umi == max(n_umi)) %>% 
  dplyr::filter(n() == 1) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(cell, larry_bc) %>% 
  dplyr::group_by(cell) %>% 
  dplyr::arrange(larry_bc) %>% 
  dplyr::summarise(larry_bc = paste(larry_bc, collapse = "__")) %>% 
  tibble::deframe()

seurat$larry <- larry_bc_calls


#-----------------------------------------------------------------------------------------------------------------------
# Filter larry matrix and save to new rds
#-----------------------------------------------------------------------------------------------------------------------
# Add barcode matrix to seurat object. Since not all cells are present in the barcode matrix, we will have to manually
# add the missing cells and set 0 to all barcode UMIs.
larry_filt_mat <- larry_filt %>% 
  mutate(cell = factor(cell, levels = Cells(seurat))) %>% 
  dplyr::select(cell, larry_bc, n_umi) %>% 
  pivot_wider(values_from = n_umi, names_from = cell, values_fill = 0) %>% 
  column_to_rownames("larry_bc") %>% 
  as.matrix()

# Create empty matrix for cells without larry barcodes
barcode_names <- unique(larry_filt$larry_bc)
cell_names_noLarry <- Cells(seurat)[!Cells(seurat) %in% larry_filt$cell]

mtx_noLarry_cells <- matrix(0, 
       length(barcode_names),
       length(cell_names_noLarry),
       dimnames = list(
         barcode_names,
         cell_names_noLarry
       )
)

# Merge both matrices
final_larry_mtx <- cbind(larry_filt_mat, mtx_noLarry_cells) %>% 
  as("dgCMatrix")

# Add filtered matrix to seurat object
seurat[["Larry_filt"]] <- CreateAssayObject(counts = final_larry_mtx)

saveRDS(seurat, file = output_rds)