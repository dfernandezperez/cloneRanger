
# Visualization 2: distribution of barcodes per cell ----------------------

## Raw --------------------------------------------------------------------

OCIAML3_barcodesXcell <- rbindlist(lapply(OCIAML3_barcodes_df, 
                                          function(barcode) {
  df <- aggregate(barcodes ~ cells, barcode, FUN = length)
  df <- df[order(df$barcodes, decreasing = TRUE), ]
  df$index <- 1:nrow(df)
  df
}), idcol = "Barcode")

ggplot(OCIAML3_barcodesXcell, aes(x = index, y = barcodes, color = Barcode)) + 
  geom_point(size = 0.1) + facet_wrap(~ Barcode, scales = "free") + 
  scale_color_manual(values = c("green", "blue", "red")) + 
  theme_bw() + scale_y_log10() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 60, vjust = 0.5)) +
  xlab("ranking of cells") + ylab("number of barcodes") + 
  ggtitle("Number of barcodes in each cell")


## counts > 1 -------------------------------------------------------------

OCIAML3_barcodesXcell <- rbindlist(lapply(OCIAML3_barcodes_df, 
                                          function(barcode) {
  barcode <- barcode[barcode$counts > 1, ]
  df <- aggregate(barcodes ~ cells, barcode, FUN = length)
  df <- df[order(df$barcodes, decreasing = TRUE), ]
  df$index <- 1:nrow(df)
  df
}), idcol = "Barcode")

ggplot(OCIAML3_barcodesXcell, aes(x = index, y = barcodes, color = Barcode)) + 
  geom_point(size = 0.1) + facet_wrap(~ Barcode, scales = "free") + 
  scale_color_manual(values = c("green", "blue", "red")) + 
  theme_bw() + scale_y_log10() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 60, vjust = 0.5)) +
  xlab("ranking of cells") + ylab("number of barcodes") + 
  ggtitle("Number of barcodes in each cell (counts > 1)")


## counts > threshold -----------------------------------------------------

OCIAML3_barcodesXcell <- rbindlist(lapply(OCIAML3_barcodes_df, 
                                          function(barcode) {
  threshold <- ceiling(quantile(as.numeric(names(table(barcode$counts))), 
                                q_counts_OCIAML3))
  barcode <- barcode[barcode$counts >= threshold, ]
  df <- aggregate(barcodes ~ cells, barcode, FUN = length)
  df <- df[order(df$barcodes, decreasing = TRUE), ]
  df$index <- 1:nrow(df)
  df
}), idcol = "Barcode")

ggplot(OCIAML3_barcodesXcell, aes(x = index, y = barcodes, color = Barcode)) + 
  geom_point(size = 0.1) + facet_wrap(~ Barcode, scales = "free") + 
  scale_color_manual(values = c("green", "blue", "red")) + 
  theme_bw() + scale_y_log10() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 60, vjust = 0.5)) +
  xlab("ranking of cells") + ylab("number of barcodes") + 
  ggtitle("Number of barcodes in each cell (counts >= individual threshold)",
          paste(thresholds))


# boxplot of counts peer number of integrations  --------------------------

OCIAML3_barcodesXcell <- rbindlist(lapply(OCIAML3_barcodes_df, 
                                          function(barcode) {
  threshold <- ceiling(quantile(as.numeric(names(table(barcode$counts))), 
                                q_counts_OCIAML3))
  barcode <- barcode[barcode$counts >= threshold, ]
  df <- aggregate(barcodes ~ cells, barcode, FUN = length)
  colnames(df)[2] <- "n_barcodes"
  df[, 2] <- as.character(df[, 2])
  df[, 2] <- ifelse(as.numeric(df[, 2]) > 5, ">5", df[, 2])
  df[, 2] <- factor(df[, 2], levels = c("1", "2", "3", "4", "5", ">5"))
  merge(barcode, df)
}), idcol = "Barcode")

ggplot(OCIAML3_barcodesXcell, 
       aes(x = n_barcodes, y = counts, fill = Barcode)) +
  geom_violin() + scale_fill_manual(values = c("green", "blue", "red")) + 
  theme_bw() + scale_y_log10() +
  theme(legend.position = "top") +
  xlab("cells by number of barcodes integrated") + ylab("counts") + 
  ggtitle("Distribution of counts per number of integrations (counts >= individual threshold)",
          paste(thresholds))
