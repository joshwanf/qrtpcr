suppressPackageStartupMessages(require("tidyverse"))
suppressPackageStartupMessages(require("readxl"))


rebind <- function(name, value, env = .GlobalEnv) {
  if (identical(env, emptyenv())) {
    stop("Can't find ", name, call. = FALSE)
  } else if (!exists(name, envir = env, inherits = FALSE)) {
    assign(name, value, envir = env)
  }
}

read_pcr_data <- function(path, results_sheet = "Results") {
  data <- normalizePath(path, mustWork = T)
  sheets <- excel_sheets(data)
  if (results_sheet %in% sheets) {
    DataRowsIndex <- which(!is.na(suppressMessages(read_excel(path, 
                                                              col_names=F, 
                                                              sheet = results_sheet)[,3]
    )))
    data <- read_excel(
      path, 
      sheet = results_sheet,
      range=cell_rows(DataRowsIndex), 
      col_types = "text"
    )
  } else {
    stop_message <- paste(
      "The sheet name used to find Results:", results_sheet
    )
    message(stop_message)
    stop("Couldn't find the results sheet, check that the Results sheet name is correct.")
  }
  return(data)
}

check_sample_names <- function(data_samples, meta_samples) {
  checked <- unique(data_samples) %>% map_lgl(
    \(x) x %in% meta_samples
  ) %>%
    all()
  return(checked)
}

remove_outliers <- function(x) {
  sdx <- sd(x, na.rm = T)
  if (sdx < 0.5 | is.na(sdx)) {
    return(x)
  } else {
    # Remove outliers
    dobj <- dist(na.omit(x))
    dmat <- as.matrix(dobj)
    diag(dmat) <- NA
    dmin <- apply(dmat, 2, min, na.rm=TRUE) # Find min distance between objects
    if (var(dmin) %in% c(0, NA)) {
      # Leave unchanged if all points are equally far apart
      return(x)
    } else {
      # Find the point that is furthest away
      id <- which(dmin==max(dmin), arr.ind=T)
      if (sd(x[-id], na.rm = T) > 0.5) {
        return(x)
      } else {
        x[id] <- NA
        return(x)
      }
    }
  }
}

find_col_name <- function(data, value_in_column) {
  # Finds an exact match of character vector `value_in_column` in the columns of `data` and returns the column names
  col_names <- c()
  for (condition in value_in_column) {
    val_regex <- paste0("^", condition, "$")
    col_contains_val <- map_lgl(data, function(x) any(str_which(x, val_regex)))
    col_names <- c(col_names, names(data)[col_contains_val])
  }
  return(col_names)
}

flatten_pcr_table <- function(data) {
  # If a [r x c] tibble contains a column that is a [r x n] tibble (not a list-column),
  # then this flattens the column and brings n columns to the top level
  flattened_table <- tibble(.rows = dim(data)[1])
  for (col in names(data)) {
    if (typeof(data[[col]]) != "list") {
      flattened_table <- add_column(flattened_table, data[col])
    } else if (typeof(data[[col]]) == "list") {
      col_df <- pull(data, col)
      flattened_table <- add_column(flattened_table, col_df)
    }
  }
  return(flattened_table)
}

calc_dct <- function(data, target, housekeeping, grouping, n_samples, debug) {
  # Typically n_housekeeping should be 3
  n_housekeeping <- length(data[["Ct mean"]][data[[target]] == housekeeping])
  n_targets <-  n_samples / n_housekeeping
  housekeeping_ct <- na.omit(data[["Ct mean"]][data[[target]] == housekeeping])
  if (length(housekeeping_ct) > 1) {
    housekeeping_ct <- mean(housekeeping_ct)
  }
  dct <- tibble(dCt = data[["Ct mean"]] - housekeeping_ct)
  if (debug) {
    dct <- dct %>% mutate(
      `dCt group` = interaction(data[[grouping]], sep = " + "),
      `dCt housekeeping` = rep(c(housekeeping_ct, rep(NA, n_housekeeping - 1)), n_targets)
    )
    message(yellow("  > Finished finding dCt."))
  }
  return(dct)
}

calc_ddct <- function(data, meta, control, independent_x_calc, grouping, n_samples, debug) {
  tmp_data <- data
  for (condition in control) {
    condition_col <- find_col_name(meta, condition)
    if (is.null(independent_x_calc) || (condition_col != independent_x_calc)) {
      tmp_data <- tmp_data %>% filter(!!sym(condition_col) == condition)
    }
  }
  control_dct <- mean(tmp_data$dct$dCt, na.rm = T)
  
  ddct <- tibble(
    ddCt = data$dct$dCt - control_dct,
    `Fold change` = 2 ^ - ddCt
  )
  if (debug) {
    ddct <- ddct %>% mutate(
      `ddCt group` = interaction(data[[grouping]], sep = " + "),
      `ddCt control` = control_dct
    )
    message(yellow("  > Finished finding ddCt and Fold change."))
  }
  return(ddct)
}

find_pcr_fc <- function(datafile,
                        meta = Meta,
                        control = "Sample 1",
                        graph_group = NULL,
                        independent_x_calc = NULL,
                        housekeeping = "GAPDH",
                        results_sheet = "Results",
                        debug = FALSE,
                        ...) {
  if (debug) message("********** Processing data to find fold change.")
  if (debug) message(yellow("  Finding fold change."))
  fc_vars <- c("fc_vars")
  
  # If datafile contains more than one file, read them and bind_rows
  if (length(datafile) > 1) {
    names(datafile) <- basename(datafile)
    list_of_data <- datafile %>% map(read_pcr_data, results_sheet)
    raw_data <- bind_rows(list_of_data, .id = "Plate")
    fc_vars <- c(fc_vars, "list_of_data", "raw_data")
  } else {
    raw_data <- read_pcr_data(datafile, results_sheet)
    fc_vars <- c(fc_vars, "raw_data")
  }
  if (debug) message(yellow("  > Finished reading data."))
  
  
  # Find column names of important columns from `data` and `meta`
  col <- list()
  col$control <- find_col_name(meta, control)
  if (length(col$control) != length(control)) {
    stop("Something wrong with the controls.")
  }
  col$target <- find_col_name(raw_data, housekeeping)
  if (length(col$target) != 1) {
    stop(paste0("Couldn't find the Target column using housekeeping value `", housekeeping, "`."))
  }
  col$sample <- grep("sample", names(raw_data), ignore.case = T, value = T)
  col$ct <- grep("^c[t|т](\\.|$)", names(raw_data), ignore.case = T, value = T)
  fc_vars <- c(fc_vars, "col")
  if (debug) message(yellow("  > Finished finding column names."))
  
  # Check sample names match between data and meta tables
  if (!(col$sample %in% names(meta))) {
    stop(paste0("Could not find column `", col$sample, "` in the Meta table."))
  }
  sample_names_checked <- check_sample_names(raw_data[[col$sample]], meta[[col$sample]])
  if (!sample_names_checked) {
    stop("Sample names didn't match between datafile and Meta table")
  }
  
  # Clean data
  # # Convert Ct to numeric, convert "Undetermined" values to 40, create factor levels
  # # Remove Ct values if omit is selected on qRT-PCR machine
  cleaned_data <- raw_data %>%
    select(contains("Plate"), matches(c("well", "omit", "sample", "target", "^c[t|т](\\.|$)"), perl = T)) %>%
    filter(!is.na(!!sym(col$sample))) %>%
    mutate(`Ct removed` = as.numeric(!!sym(col$ct))) %>%
    mutate(`Ct removed` = replace_na(`Ct removed`, 40)) %>%
    mutate(across(where(is.character), as_factor))
  fc_vars <- c(fc_vars, "cleaned_data")
  if (debug) message(yellow("  > Finished cleanning data."))
  
  if ("Omit" %in% names(cleaned_data)) {
    cleaned_data <- cleaned_data %>%
      mutate(`Ct removed` = case_when(Omit == "true" ~ NA_real_,
                                      TRUE ~ `Ct removed`)
      )
    if (debug) message(yellow("  > Finished omitting data."))
  }
  
  # Mean data
  # # Remove Ct outliers if sd(Ct) > 0.5 and only one value is removed
  # # Average technical replicates and find dCt
  
  debug_data <- list()
  groupings <- list()
  exprs <- list()
  groupings$mean_ct <- c(col$target, col$sample)
  if (debug) {
    debug_data$meaned_data <- cleaned_data %>%
      group_by(!!!syms(groupings$mean_ct)) %>%
      mutate(`mean groups (Target, Sample)` = paste(cur_group(), collapse = ", "),
             `n per group/total groups` = paste0(n(), "/", n_groups(.)),
             `mean group id` = cur_group_id(),
             `mean wells` = paste(cur_data()[["Well Position"]], collapse = ", "),
             `mean formula` = paste(na.omit(cur_data()[["Ct removed"]]), collapse = ", ")
      ) %>%
      ungroup()
    message(yellow("  > Finished setting cleaned_data debug."))
  }
  meaned_data <- cleaned_data %>%
    group_by(!!!syms(groupings$mean_ct)) %>%
    mutate(
      `Ct removed` = remove_outliers(`Ct removed`),
      `Ct mean` = case_when(all(is.na(`Ct removed`)) ~ `Ct removed`,
                            TRUE ~ c(mean(`Ct removed`, na.rm = T), rep(NA, n()-1))
      )
    ) %>%
    ungroup()
  fc_vars <- c(fc_vars, "debug_data", "groupings", "meaned_data")
  if (debug) message(yellow("  > Finished meaning data."))
  
  # dCt data
  groupings$dct <- c(col$sample)
  dct_data <- meaned_data %>%
    group_by(!!!syms(groupings$dct)) %>%
    mutate(
      dct = calc_dct(cur_data_all(), col$target, housekeeping, groupings$dct, n(), debug)
    ) %>%
    ungroup()
  fc_vars <- c(fc_vars, "dct_data")
  
  
  joined_data <- left_join(dct_data, meta, by = col$sample[1]) # col$sample[1] in case there is more than one column
  fc_vars <- c(fc_vars, "joined_data")
  if (debug) message(yellow("  > Finished joining data."))
  
  # Final data
  # # Calculate the ddct, fold change, and sd
  groupings$ddct <- c(col$target, independent_x_calc)
  final_data <- joined_data %>%
    group_by(!!!syms(groupings$ddct)) %>%
    mutate(
      ddct = calc_ddct(cur_data_all(), meta, control, independent_x_calc, groupings$ddct, n(), debug)
    ) %>%
    ungroup()
  
  groupings$sd <- c(col$target, col$control)
  final_data <- final_data %>%
    group_by(!!!syms(groupings$sd)) %>%
    mutate(
      `Mean fold change` = c(mean(ddct$`Fold change`, na.rm = T), rep(NA, n() - 1)),
      sd = c(sd(ddct$`Fold change`, na.rm = T), rep(NA, n() - 1))
    ) %>%
    ungroup() %>%
    select(
      matches("well"), 
      colnames(meta),
      everything()
    ) %>%
    flatten_pcr_table()
  
  # Rearrange primers so the housekeeping rows are displayed first
  final_data[[col$target]] <- relevel(final_data[[col$target]], ref = housekeeping)
  final_data <- final_data %>% arrange(!!sym(col$target))
  fc_vars <- c(fc_vars, "final_data")
  
  if (debug) message("********** Successfully processed data to find fold change.")
  return(final_data)
}

plot_pcr <- function(data,
                     meta,
                     control = "Sample 1",
                     graph_group = NULL,
                     independent_x_calc = NULL,
                     housekeeping = "GAPDH",
                     graph_cols = 2, 
                     debug = F) {
  data <- data %>%
    select(-matches(c("dCt group", "dCt housekeeping", "ddCt group", "ddCt control")))
  if (debug) message("********** Processing data to plot graph.")
  plot_vars <- c("plot_vars")
  # Set graph variables
  # Find column names of important columns from `data` and `meta`
  col <- list()
  col$control <- find_col_name(meta, control)
  col$target <- find_col_name(data, housekeeping)
  col$sample <- grep("sample", names(data), ignore.case = T, value = T)
  col$ct <- grep("^c[t|т](\\.|$)", names(data), ignore.case = T, value = T)
  
  # Decide how to draw the graph
  
  # Set column to use for x-axis
  graph_x <- setdiff(col$control, graph_group)
  if (length(graph_x) != 1) {
    stop_message <- paste(
      c(
        paste("Check `control` and `graph_group` variables."),
        paste("`control` columns:", paste(col$control, collapse = " + ")),
        paste("`graph_group` value:", graph_group)
      ),
      collapse = "\r\n"
    )
    stop(stop_message)
  }
  
  if (!is.null(graph_group)) {
    
    # Set columns to use for legend groups
    if (length(graph_group) > 1) {
      combined_graph_group <- paste(graph_group, collapse = " - ")
      data <- data %>%
        mutate(!!combined_graph_group := interaction(data[graph_group], sep = " - ", lex.order = T))
      graph_group_title <- sym(combined_graph_group)
    } else if (length(graph_group) == 1) {
      graph_group_title <- sym(graph_group)
    }
  } else {
    graph_group_title <- graph_group
  }
  plot_vars <- c(plot_vars, "col", "graph_x", "graph_group_title")
  
  # Filter housekeeping data 
  # and any rows with NA in the "Fold change" column are gone
  na_removed_data <- data %>% 
    filter(!!sym(col$target) != housekeeping) %>%
    filter(!is.na(`Fold change`)) %>%
    mutate(plot_page = ceiling((as.integer(!!sym(col$target)) - 1) / 4))
  plot_vars <- c(plot_vars, "na_removed_data")
  
  # Plot
  pcr_plot <- list()
  for (page in unique(na_removed_data$plot_page)) {
    pcr_plot[[page]] <- na_removed_data %>%
      filter(plot_page == page) %>%
      ggplot(
        aes(
          x = !!sym(graph_x),
          y = `Mean fold change`,
          fill = !!graph_group_title,
          ymin = `Mean fold change` - sd,
          ymax = `Mean fold change` + sd
        )
      ) +
      facet_wrap(sym(col$target), ncol = graph_cols, scales = "free_y") +
      geom_col(position = "dodge", na.rm = T) +
      geom_errorbar(width = 0.5,
                    position = position_dodge(width=0.9)) +
      # geom_jitter(aes(y = `Fold change`), width = 0.3) +
      geom_point(
        aes(y = `Fold change`), 
        position = position_dodge(width = 0.9)) +
      # position = position_jitterdodge(jitter.width = 0.3)) +
      labs(y = "Fold change") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  
  if (debug) message("********** Sucessfully processed data to plot graph.")
  return(pcr_plot)
}

prep_pcr_analysis_to_write <- function(analyzed_pcr_data, housekeeping) {
  col <- list()
  col$target <- find_col_name(analyzed_pcr_data, housekeeping)
  if (length(col$target) != 1) {
    stop_message <- c(
      paste("The column name used to find the Targets:", col$target),
      paste("The housekeeping name used:", housekeeping)
    )
    message(stop_message)
    stop("Couldn't find the Target column, check the housekeeping name")
  }
  primers <- unique(analyzed_pcr_data[[col$target]]) %>% as.character()
  primers <- primers[primers != housekeeping]
  separated_data <- primers %>% 
    set_names() %>%
    map(
      \(x) AnalyzedPCRData %>% filter(!!sym(col$target) %in% c(housekeeping, x))
    )
  return(separated_data)
}

find_pcr_stats <- function(analyzed_pcr_data, Meta, control, housekeeping) {
  library(rstatix)
  col <- list()
  col$control <- find_col_name(Meta, control)
  if (length(col$control) != length(control)) {
    stop("Something wrong with the controls.")
  }
  col$target <- find_col_name(analyzed_pcr_data, housekeeping)
  if (length(col$target) != 1) {
    stop_message <- c(
      paste("The column name used to find the Targets:", col$target),
      paste("The housekeeping name used:", housekeeping)
    )
    message(stop_message)
    stop("Couldn't find the Target column, check the housekeeping name")
  }
  primers <- unique(analyzed_pcr_data[[col$target]]) %>% as.character()
  primers <- primers[primers != housekeeping]
  ttest <- analyzed_pcr_data %>%
    filter(!!sym(col$target) != housekeeping) %>%
    group_by(!!sym(col$target), !!!syms(col$control)) %>%
    mutate(ttest_group = interaction(!!!syms(col$control))) %>%
    group_by(!!(sym(col$target))) %>%
    t_test(`Fold change` ~ ttest_group, ref.group = paste(control, collapse = ".")) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()
  return(ttest)
}
