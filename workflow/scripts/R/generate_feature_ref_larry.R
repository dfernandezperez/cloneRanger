log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(tidyverse)

#------------------------------------------------------------------------------------------
# Load feature ref from every sample and merge it into 1
#------------------------------------------------------------------------------------------
feature_ref <- snakemake@input %>%
  purrr::map(read_csv) %>% 
  bind_rows() %>% 
  mutate(
    name = id
  ) %>% 
  distinct() %>% 
  group_by(id) %>% 
  mutate(
    id   = paste0(id, "_", 1:n()),
    name = id
  ) %>%
  select(id, name, read, pattern, sequence, feature_type)

write_csv(feature_ref, snakemake@output[[1]])