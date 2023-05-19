#-----------------------------------------------------------------------------------------------------------------------
# Put a list of objects as tabs in rmarkdown
#-----------------------------------------------------------------------------------------------------------------------
in_tabs <- function(l, labels = names(l), level, knit = TRUE, close_tabset = FALSE) {
  require(stringi)
  id <- stri_rand_strings(1, 15) # generate random id for chunk
  ## https://stackoverflow.com/questions/69353667/dynamic-creation-of-tabs-in-rmarkdown-does-not-work-for-ggplot-while-it-does-for
  if (is.null(labels)) {
    stop("labels are NULL, it is required not to be so that the tabs have proper names")
  }
  names(l) <- labels

  rmd_code <- lapply(seq_along(l), FUN = function(i) obj_to_rmd(l[[i]], name = names(l)[i], level = level + 1L, id = paste(id, i, sep = "-")))

  if (isTRUE(getOption("knitr.in.progress"))) {
    res <- knitr::knit(text = unlist(rmd_code), quiet = TRUE)
    cat(res)
  } else {
    if (!knit) {
      cat(unlist(rmd_code))
    } else {
      return(l)
    }
  }
  if (close_tabset) {
    cat(paste(get_section(level), "{.unlisted .unnumbered .toc-ignore .tabset}", "\n"))
  }
}

get_section <- function(level) {
  paste(rep("#", times = level), collapse = "")
}

get_tabset <- function(obj) {
  if (inherits(obj, "list")) "{.tabset}" else ""
}

obj_to_rmd <- function(obj, parent_name = "l", name, level, id) {
  section_code <- sprintf("%s %s %s\n", get_section(level), name, get_tabset(obj))
  if (!inherits(obj, "list")) {
    rmd_code <- c(
      sprintf("```{r plot-%s, echo = FALSE}\n", id),
      sprintf("%s$`%s`\n", parent_name, name),
      "```\n",
      "\n"
    )
  } else {
    rmd_code <- c(
      "\n",
      lapply(
        X = seq_along(obj),
        FUN = function(i) obj_to_rmd(obj[[i]], sprintf("%s$`%s`", parent_name, name), names(obj)[i], level + 1L)
      )
    )
  }
  return(c(section_code, rmd_code))
}

#-----------------------------------------------------------------------------------------------------------------------
# Plot genes with featureplot present in a given list of genes
#-----------------------------------------------------------------------------------------------------------------------
plot_detected_genes <- function(seurat, detected_genes, gene) {
  if (gene %in% detected_genes) {
    p <- FeaturePlot(seurat, features = gene, max.cutoff = "q99", raster = FALSE) + scale_color_viridis()
  } else {
    p <- "Not detected."
  }
  return(p)
}