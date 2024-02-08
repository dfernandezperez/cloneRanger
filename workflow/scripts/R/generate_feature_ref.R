log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(tidyverse)

#------------------------------------------------------------------------------------------
# Load feature ref from larry/cellhashing 
#------------------------------------------------------------------------------------------
if ("larry_ref" %in% names(snakemake@input) & "cell_hash_ref" %in% names(snakemake@input)) {

	# If sample has larry and cellhashing info, read both and merge them into 1 ref file
	larry_ref    <- read_csv(snakemake@input[["larry_ref"]])
	cellhash_ref <- read_csv(snakemake@input[["cell_hash_ref"]])

	combined_ref <- bind_rows(larry_ref, cellhash_ref)
	write_csv(combined_ref, snakemake@output[[1]])

} else { 
	
	# If just cellhash or larry are set, return the same csv
	feature_ref <- read_csv(snakemake@input[[1]])
	write_csv(feature_ref, snakemake@output[[1]])

}