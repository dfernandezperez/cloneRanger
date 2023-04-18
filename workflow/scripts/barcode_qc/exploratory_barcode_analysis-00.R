
library(Seurat)
library(data.table)
library(BiocParallel)
library(stringr)
library(ggplot2)

"%ni%" <- Negate("%in%")
source("../src/clone_calling_functions.R")

# Import data -------------------------------------------------------------

OCIAML3 <- readRDS("../data/seurat_objects/01_OCIAML3.rds")


# Specify metadata --------------------------------------------------------

barcodes_names <- c("GFP", "Sapphire", "Scarlet")
names(barcodes_names) <- barcodes_names
color_barcodes <- c("green", "blue", "red")
samples_names <- c("Control_1", "Control_2", "Treatment_1", "Treatment_2")
names(samples_names) <- samples_names
samples_acronyms <- paste(substr(samples_names, 1, 1), 
                          substr(samples_names, nchar(samples_names), 
                                 nchar(samples_names)), 
                          sep = "")


# Extract barcode data (no filtering) -------------------------------------

OCIAML3_barcodes_df <- get_bc_df(split_matrices_barcodes(OCIAML3, barcodes_names))


# Save thresholds of counts for later filtering ---------------------------

q_counts_OCIAML3 <- 0.075
OCIAML3_barcodes_df_threshold <- rbindlist(lapply(OCIAML3_barcodes_df, 
                                                  function(barcode) {
  data.frame(cutoff = ceiling(quantile(as.numeric(names(table(barcode$counts))), 
                                       q_counts_OCIAML3)))                                                  
}), idcol = "Barcode")
thresholds <- as.character(OCIAML3_barcodes_df_threshold$cutoff)
thresholds <- paste(thresholds, collapse = ", ")


# Summary metrics ---------------------------------------------------------

OCIAML3_barcodes_df_summ <- rbindlist(lapply(OCIAML3_barcodes_df, function(df) {
  data.frame(n_barcodes = nrow(df), 
             q0.25 = quantile(df$counts[df$counts > 1], 1/4), 
             median = quantile(df$counts[df$counts > 1], 1/2), 
             q0.75 = quantile(df$counts[df$counts > 1], 3/4))
}))
OCIAML3_barcodes_df_summ


