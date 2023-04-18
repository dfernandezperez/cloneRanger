
split_matrices_barcodes <- function(seurat_object, barcodes_names) {
  
  sapply(barcodes_names, function(barcode) {
    
    bc <- str_detect(seurat_object@assays[["Barcodes"]]@data@Dimnames[[1]], 
                     barcode)
    
    matrix <- seurat_object@assays[["Barcodes"]]@data[bc, ]
    
    # Barcodes with positive counts in at least one (>= 1) cell
    
    matrix[Matrix::rowSums(matrix > 0) >= 1, ] 
    
  }, simplify = FALSE)
}


get_bc_df <- function(matrices_list) {
  
  lapply(matrices_list, function(M){
    
    M <- as(M, "TsparseMatrix")
    
    # M@i indicates the barcode id associated to a specific cell
    # M@j indicates the cell id associated to a specific barcode
    # M@x indicates the counts of each barcode in each cell
    # Add +1 because the indexes of the dgTMatrix are 0-based
    
    barcodes_m <- M@Dimnames[[1]][M@i + 1] 
    cells_m    <- M@Dimnames[[2]][M@j + 1]
    counts_m   <- M@x
    
    df <- setNames(as.data.frame(cbind(cells_m, barcodes_m, counts_m)),
                   nm = c("cells", "barcodes", "counts"))
    df$counts <- as.numeric(df$counts)
    
    df
  })
}


utilsMultipleIntegration <- function(cells_with_barcode, reduce = TRUE) {
  
  ## If multiple barcodes as label:
  ##  retain them if found as labels of other cells, or
  ##  remove from the label if not found in other cells
  
  ## Labels with more than one barcode - "multiple"
  multiple_bc <- cells_with_barcode$barcode[str_detect(
    cells_with_barcode$barcode, "_")]
  names(multiple_bc) <- cells_with_barcode$cellID[
    str_detect(cells_with_barcode$barcode, "_")]
  
  multiple_bc_alone <- unlist(str_split(multiple_bc, "_"))
  
  ## Labels with just one barcode - "unique"
  unique_bc <- cells_with_barcode$barcode[!
                                            str_detect(cells_with_barcode$barcode, "_")]
  names(unique_bc) <- cells_with_barcode$cellID[!
                                                  str_detect(cells_with_barcode$barcode, "_")]
  
  ## Remove barcodes only found once in unique
  if (reduce) {
    freq_bc_in_unique <- table(unique_bc)
    unique_bc_once <- names(freq_bc_in_unique)[freq_bc_in_unique == 1]
    unique_bc_once <- unique_bc_once[unique_bc_once %ni% multiple_bc_alone]
    unique_bc <- unique_bc[unique_bc %ni% unique_bc_once]
  }
  
  ## Extract barcodes found only once in "multiple" and
  ##  not being labels in "unique"
  freq_bc_in_multiple <- table(multiple_bc_alone)
  multiple_bc_once <- names(freq_bc_in_multiple)[freq_bc_in_multiple == 1]
  multiple_bc_once <- multiple_bc_once[multiple_bc_once %ni% unique_bc]
  
  ## Remove those barcodes from the labels in "multiple"
  corrected_labels <- str_remove_all(multiple_bc, 
                                     paste(c(paste("^", multiple_bc_once, "_", 
                                                   sep = ""),
                                             paste("_", multiple_bc_once, "_",
                                                   sep = ""),
                                             paste("_", multiple_bc_once, "$",
                                                   sep = "")), 
                                           collapse = "|"))
  corrected_labels <- str_replace_all(corrected_labels, 
                                      "([:digit:])([:upper:])", "\\1_\\2")
  names(corrected_labels) <- names(multiple_bc)
  if (any(str_detect(corrected_labels, "_"))) {
    corrected_labels_1 <- corrected_labels[! str_detect(corrected_labels, "_")]
    
    ## Sort the labels to make them match as string if so
    corrected_labels_2 <- corrected_labels[str_detect(corrected_labels, "_")]
    names_corrected_labels_2 <- names(corrected_labels_2)
    corrected_labels_2 <- str_replace(str_replace_all(apply(
      apply(str_split(corrected_labels_2, "_", simplify = TRUE), 
            1, str_sort, numeric = TRUE), 
      2, paste, collapse = "_"), "__", ""), "^_", "")
    names(corrected_labels_2) <- names_corrected_labels_2
    
    df <- rbind(data.frame(cellID = names(unique_bc), barcode = unique_bc),
                data.frame(cellID = names(corrected_labels_1), 
                           barcode = corrected_labels_1),
                data.frame(cellID = names(corrected_labels_2), 
                           barcode = corrected_labels_2))
  } else {
    df <- rbind(data.frame(cellID = names(unique_bc), barcode = unique_bc),
                data.frame(cellID = names(corrected_labels), 
                           barcode = corrected_labels))
  }
  
  ## Return dataframe with cells and their labels
  if (reduce) {
    out1 <- length(unique_bc_once)
    out_multiple <- sum(nchar(corrected_labels) == 0)
    list(df = df, 
         cells_with_unique_labels = setNames(c(out1, out_multiple),
                                             c("unique", "multiple")))
  } else {
    df
  }
}



