# MetaEngine reporting module.

pm_result_summary <- function(results,
                              p_thresholds = c(0.05, 0.01, 0.001),
                              fdr_thresholds = c(0.10, 0.05, 0.01),
                              fc_thresholds = c(1.5, 2, 3),
                              p_col = "P.Value",
                              fdr_col = "adj.P.Val",
                              logfc_col = "logFC") {
  
  res <- pm_read(results)
  
  required <- c(p_col, fdr_col, logfc_col)
  missing_cols <- setdiff(required, names(res))
  if (length(missing_cols)) {
    stop("Differential-results table is missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }
  
  p <- as.numeric(res[[p_col]])
  fdr <- as.numeric(res[[fdr_col]])
  logfc <- as.numeric(res[[logfc_col]])
  abs_fc <- 2^abs(logfc)
  
  rows <- list()
  
  add_row <- function(category, criterion, count, up = NA_integer_, down = NA_integer_) {
    rows[[length(rows) + 1L]] <<- data.frame(
      category = category,
      criterion = criterion,
      count = as.integer(count),
      up = as.integer(up),
      down = as.integer(down),
      stringsAsFactors = FALSE
    )
  }
  
  add_row("Features", "Total tested features", nrow(res))
  
  for (thr in p_thresholds) {
    hit <- is.finite(p) & p < thr
    add_row("Raw P",
            paste0("P < ", format(thr, scientific = FALSE)),
            sum(hit),
            sum(hit & logfc > 0, na.rm = TRUE),
            sum(hit & logfc < 0, na.rm = TRUE))
  }
  
  for (thr in fdr_thresholds) {
    hit <- is.finite(fdr) & fdr < thr
    add_row("FDR",
            paste0("FDR < ", format(thr, scientific = FALSE)),
            sum(hit),
            sum(hit & logfc > 0, na.rm = TRUE),
            sum(hit & logfc < 0, na.rm = TRUE))
  }
  
  for (pthr in p_thresholds) {
    for (fcthr in fc_thresholds) {
      hit <- is.finite(p) & p < pthr & is.finite(abs_fc) & abs_fc >= fcthr
      add_row("Raw P + FC",
              paste0("P < ", format(pthr, scientific = FALSE),
                     " & |FC| >= ", fcthr),
              sum(hit),
              sum(hit & logfc > 0, na.rm = TRUE),
              sum(hit & logfc < 0, na.rm = TRUE))
    }
  }
  
  for (qthr in fdr_thresholds) {
    for (fcthr in fc_thresholds) {
      hit <- is.finite(fdr) & fdr < qthr & is.finite(abs_fc) & abs_fc >= fcthr
      add_row("FDR + FC",
              paste0("FDR < ", format(qthr, scientific = FALSE),
                     " & |FC| >= ", fcthr),
              sum(hit),
              sum(hit & logfc > 0, na.rm = TRUE),
              sum(hit & logfc < 0, na.rm = TRUE))
    }
  }
  
  do.call(rbind, rows)
}

pm_print_summary <- function(summary_table,
                             n_samples_total = NA_integer_,
                             n_samples_used = NA_integer_,
                             n_samples_excluded = NA_integer_,
                             model = NULL,
                             contrast = NULL) {
  
  cat("\n============================================================\n")
  cat("                 MetaEngine Summary\n")
  cat("============================================================\n")
  
  if (!is.null(model)) {
    cat("\nModel:\n")
    cat(paste(deparse(model), collapse = " "), "\n")
  }
  
  if (!is.null(contrast)) {
    cat("\nContrast:\n")
    cat(as.character(contrast), "\n")
  }
  
  if (is.finite(n_samples_total)) {
    cat("\nSamples\n------------------------------------------------------------\n")
    cat(sprintf("%-34s %8d\n", "Total samples:", n_samples_total))
    if (is.finite(n_samples_used))
      cat(sprintf("%-34s %8d\n", "Samples analyzed:", n_samples_used))
    if (is.finite(n_samples_excluded))
      cat(sprintf("%-34s %8d\n", "Samples excluded:", n_samples_excluded))
  }
  
  cat("\nDifferential features\n------------------------------------------------------------\n")
  
  for (i in seq_len(nrow(summary_table))) {
    if (summary_table$category[i] == "Features") next
    
    cat(sprintf("%-40s %7d",
                summary_table$criterion[i],
                summary_table$count[i]))
    
    if (!is.na(summary_table$up[i]) && !is.na(summary_table$down[i])) {
      cat(sprintf("   (Up: %d, Down: %d)",
                  summary_table$up[i],
                  summary_table$down[i]))
    }
    
    cat("\n")
  }
  
  cat("============================================================\n\n")
  invisible(summary_table)
}

pm_write_results_summary <- function(results,
                                     output_dir,
                                     p_thresholds = c(0.05, 0.01, 0.001),
                                     fdr_thresholds = c(0.10, 0.05, 0.01),
                                     fc_thresholds = c(1.5, 2, 3),
                                     n_samples_total = NA_integer_,
                                     n_samples_used = NA_integer_,
                                     n_samples_excluded = NA_integer_,
                                     model = NULL,
                                     contrast = NULL) {
  
  summary_tab <- pm_result_summary(
    results,
    p_thresholds = p_thresholds,
    fdr_thresholds = fdr_thresholds,
    fc_thresholds = fc_thresholds
  )
  
  write.csv(summary_tab,
            file.path(output_dir, "differential_summary.csv"),
            row.names = FALSE)
  
  txt <- capture.output(
    pm_print_summary(
      summary_tab,
      n_samples_total = n_samples_total,
      n_samples_used = n_samples_used,
      n_samples_excluded = n_samples_excluded,
      model = model,
      contrast = contrast
    )
  )
  
  writeLines(txt, file.path(output_dir, "differential_summary.txt"))
  
  get_exact_count <- function(label) {
    z <- summary_tab$count[summary_tab$criterion == label]
    if (!length(z)) return(NA_integer_)
    z[1]
  }
  
  n_p05 <- get_exact_count("P < 0.05")
  n_fdr05 <- get_exact_count("FDR < 0.05")
  n_p05_fc2 <- get_exact_count("P < 0.05 & |FC| >= 2")
  n_fdr05_fc2 <- get_exact_count("FDR < 0.05 & |FC| >= 2")
  
  manuscript <- paste0(
    "Differential analysis identified ", n_p05,
    " features at nominal P < 0.05, including ", n_p05_fc2,
    " with absolute fold change >= 2. ",
    n_fdr05,
    " features remained significant at FDR < 0.05, including ",
    n_fdr05_fc2,
    " with absolute fold change >= 2."
  )
  
  writeLines(manuscript,
             file.path(output_dir, "results_summary_manuscript.txt"))
  
  invisible(list(summary = summary_tab, manuscript = manuscript))
}



