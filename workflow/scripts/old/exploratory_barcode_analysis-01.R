
# Visualization 1: count distribution of barcode/cell ---------------------

## Raw --------------------------------------------------------------------

OCIAML3_barcodes_df_distrib <- rbindlist(lapply(OCIAML3_barcodes_df, 
                                                function(barcode) {
  df <- data.frame(counts = sort(barcode[, 3],
  decreasing = TRUE))
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
  xlab("rank of expression") +
  ggtitle("Expression of barcode/cell")


## counts > 1 -------------------------------------------------------------

OCIAML3_barcodes_df_distrib <- rbindlist(lapply(OCIAML3_barcodes_df, 
                                                function(barcode) {
  df <- data.frame(counts = sort(barcode[barcode$counts > 1, 3],
                                 decreasing = TRUE))
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
  xlab("rank of expression") +
  ggtitle("Expression of barcode/cell (counts > 1)")


## distribution of counts (violin) - counts > 1 ---------------------------

ggplot(OCIAML3_barcodes_df_distrib, 
       aes(x = Barcode, y = counts, fill = Barcode)) + 
  geom_violin() + 
  scale_fill_manual(values = c("green", "blue", "red")) + 
  theme_bw() + scale_y_log10() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  xlab("rank of expression barcode/cell") +
  ggtitle("Distribution of counts in barcode/cell (counts > 1)")


## counts > threshold -----------------------------------------------------

OCIAML3_barcodes_df_distrib <- rbindlist(lapply(OCIAML3_barcodes_df, 
                                                function(barcode) {
  threshold <- ceiling(quantile(as.numeric(names(table(barcode$counts))), 
                                q_counts_OCIAML3))
  df <- data.frame(counts = sort(barcode[barcode$counts >= threshold, 3],
                                 decreasing = TRUE))
  
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
  ggtitle("Expression of barcode/cell (counts >= individual threshold)",
          paste(thresholds))



## distribution of counts (violin) - counts > threshold -------------------

ggplot(OCIAML3_barcodes_df_distrib, 
       aes(x = Barcode, y = counts, fill = Barcode)) + 
  geom_violin() + 
  scale_fill_manual(values = c("green", "blue", "red")) + 
  theme_bw() + scale_y_log10() + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  xlab("rank of expression barcode/cell") +
  ggtitle("Distribution of counts in barcode/cell (counts >= individual threshold)")


