---
title: "Plotting custom marker genes"
date: "`r format(Sys.time(), '%d %B, %Y')`"
output:
  html_document:
  highlight: tango
  number_sections: no
  theme: default
  toc: yes
  toc_depth: 3
  toc_float:
    collapsed: no
    smooth_scroll: yes
---

```{r logs}
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")
```

```{r libraries}
#| code-summary: Libraries
library(Seurat)
library(tidyverse)
library(viridis)
```


```{r input-files-and-params}
input_file   <- snakemake@input[[1]]
marker_genes <- snakemake@params[["marker_genes"]]
```