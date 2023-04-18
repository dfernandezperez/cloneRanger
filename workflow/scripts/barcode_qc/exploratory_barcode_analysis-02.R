
# Visualization 2: distribution of barcodes per cell ----------------------

## Raw --------------------------------------------------------------------

OCIAML3_barcodes_df_distrib <- rbindlist(lapply(OCIAML3_barcodes_df, 
                                                function(barcode) {
  df <- aggregate(barcodes ~ cells, barcode, FUN = length)
  df <- df[order(df$barcodes, decreasing = TRUE), ]
  df$index <- 1:nrow(df)
  df
}), idcol = "Barcode")

ggplot(distrib_zeros, aes(x = index, y = counts, color = Barcode)) + 
  geom_point(size = 0.1) + facet_wrap(~ Barcode, scales = "free") + 
  scale_color_manual(values = c("green", "blue", "red")) + 
  theme_bw() + scale_y_log10() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 60, vjust = 0.5)) +
  xlab("rank of expression") +
  ggtitle("Expression of barcode/cell")


## counts > 1 -------------------------------------------------------------

OCIAML3_barcodes_df_distrib <- rbindlist(lapply(OCIAML3_barcodes_df, 
                                                function(barcode) {
  barcode <- barcode[barcode$counts > 1, ]
  df <- aggregate(barcodes ~ cells, barcode, FUN = length)
  df <- df[order(df$barcodes, decreasing = TRUE), ]
  df$index <- 1:nrow(df)
  df
}), idcol = "Barcode")

ggplot(distrib_no_ones, aes(x = index, y = counts, color = Barcode)) + 
  geom_point(size = 0.1) + facet_wrap(~ Barcode, scales = "free") + 
  scale_color_manual(values = c("green", "blue", "red")) + 
  theme_bw() + scale_y_log10() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 60, vjust = 0.5)) +
  xlab("rank of expression") +
  ggtitle("Expression of barcode/cell (counts > 1)")


## counts > threshold -----------------------------------------------------

OCIAML3_barcodes_df_distrib <- rbindlist(lapply(OCIAML3_barcodes_df, 
                                                function(barcode) {
  threshold <- ceiling(quantile(as.numeric(names(table(barcode$counts))), 
                                q_counts_OCIAML3))
  barcode <- barcode[barcode$counts > threshold, ]
  df <- aggregate(barcodes ~ cells, barcode, FUN = length)
  df <- df[order(df$barcodes, decreasing = TRUE), ]
  df$index <- 1:nrow(df)
  df
}), idcol = "Barcode")

ggplot(OCIAML3_barcodes_df_distrib, 
       aes(x = index, y = counts, color = Barcode)) + 
  geom_point(size = 0.1) + facet_wrap(~ Barcode, scales = "free") + 
  scale_color_manual(values = c("green", "blue", "red")) + 
  theme_bw() + scale_y_log10() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 60, vjust = 0.5)) +
  xlab("rank of expression barcode/cell") +
  ggtitle("Expression of barcode/cell (counts > individual threshold)",
          paste(thresholds))

