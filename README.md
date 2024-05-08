# qrtpcr

Analyze qRT-PCR data using delta-delta Ct method in R

# Prerequisites

## Load packages

```r
packages <- c(
    "tidyverse",
    "readxl",
    "openxlsx"
)

tryCatch(
    {
    for (package in packages) {
      if (!(package %in% rownames(installed.packages()))) {
        message("Installing package ", package)
        install.packages(package, quiet = T)
      }
      library(package, character.only = T)
    }
  },
  error = function(e) {
    message("Ran into an error")
    message(e)
  }, 
  finally = {}
)
```

## Meta table

This is a tab-separated file containing a list of samples and conditions.
The samples column name (typically `Sample Name`) must match the column in the exported data file.
The file can have additional columns not used in the analysis but may contain extra information.

```r
meta <- tibble(
    `Sample Name` = paste("Sample", 1:12),
    `Growth media` = rep(c("Normal media", "New media"), each = 6),
    Treatment = rep(c("Vehicle", "1 uM drug", "10 uM drug"), each = 2)
)
write_tsv(meta, path = "meta.txt")
```

## Use

There are two main functions: `find_pcr_fc()` and `plot_pcr()`.
They depend on the following variables that match experiment details:

Variable            | Description
---                 | ---
metafile            | a path to the tsv meta file
datafile            | character string or vector containing a path to the data files
control             | character vector that matches the name of the control sample. This can be a single sample name (`Sample 1`) or two biological conditions (`c("Normal media", "Vehicle")`). Default value is `Sample 1`.
graph_group         | character string of a column name in `meta` that will be used for groups when plotting. Default value is `NULL`.
independent_x_calc  | character string of a column name in `meta` that will be used to calculate conditions independently. Default value is `NULL`.
housekeeping        | character string of the housekeeping gene used (case sensitive). Default value is `GAPDH`.
output_file         | character string for the output file

```r
datafiles <- c(
    "data plate 1.xlsx",
    "data plate 2.xlsx"
)
control <- c("Normal media", "Vehicle")
graph_group <- "Growth media"
source("qrtpcr.R")

meta <- read_tsv(path = metafile, col_types = cols()) %>%
    mutate(across(everything(), as_factor))

analyzed_data <- fund_pcr_fc(
    datafile = datafiles,
    meta = meta,
    control = control,
    graph_group = graph_group
)

pcr_plot <- plot_pcr(
    data = analyzed_data,
    meta = meta,
    control = control,
    graph_group = graph_group,
)
walk(pcr_plot, ~ print(.x))

write.xlsx(
    x = c(
        list(All = analyzed_data),
        prep_pcr_analysis_to_write(analyzed_data, housekeeping)
    ),
    file = output_file
)
```
