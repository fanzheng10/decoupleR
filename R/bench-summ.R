#' Function that provides summary and plots for the benchmark run
#'
#' @param .res_tibble formatted bench result tibble with added auroc column
#' @return AUROC summary with TF coverage, ROC, AUROC, PRROC, Run time,
#' ROC plots, and Heatmap plots
#' @import ggplot2
#' @import pheatmap
get_bench_summary <- function(.res_tibble) {
  # get roc results
  roc <- format_roc(.res_tibble, "roc")

  # get PR roc results
  pr_roc <- format_roc(.res_tibble, "prroc")

  # Plot ROC
  roc_plot <- ggplot(roc, aes(x = fpr, y = tpr, colour = run_key)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("FPR (1-specificity)") +
    ylab("TPR (sensitivity)")

  # Plot PR ROC
  pr_roc_plot <- ggplot(pr_roc, aes(x = recall, y = precision , colour = run_key)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("Recall/Sensitivity") +
    ylab("Precision")

  # Extract AUROC
  auroc_tibble <- .res_tibble %>%
    unnest(roc) %>%
    select(set_bench, lvls, statistic, auc) %>%
    distinct()

  # Plot AUROC
  auroc_plot <- auroc_tibble %>%
    unite("run_key", set_bench, statistic, lvls, remove = F) %>%
    ggplot(., aes(x = reorder(run_key, auc),
                  y = auc,
                  fill = run_key)) +
    geom_bar(stat = "identity") +
    xlab("networks") +
    ylab("AUROC") +
    coord_flip(ylim = c(0.5, 0.8)) +
    theme(legend.position = "none")

  # AUROC Heatmap
  auroc_heat <- auroc_tibble %>% get_auroc_heat()

  # Extract AU PRROC
  pr_auroc_tibble <- .res_tibble %>%
    unnest(prroc) %>%
    select(set_bench, lvls, statistic, auc) %>%
    distinct()

  # AU PRROC Heatmap
  prroc_heat <- pr_auroc_tibble %>% get_auroc_heat()

  # get computational time info
  comp_time <- .res_tibble %>%
    # get statistic time from activity
    mutate(statistic_time = activity %>%
             map(function(tib)
               tib %>%
                 select(statistic_time) %>%
                 unique)) %>%
    unnest(statistic_time) %>%
    # calculate regulon size
    group_by(set_bench, lvls) %>%
    mutate(regulon_time = sum(statistic_time)) %>%
    select(set_bench, lvls, statistic_time, regulon_time)


  # Join AUROC, PRROC, Coverage, and Comp time
  auroc_summary <- auroc_tibble %>%
    inner_join(pr_auroc_tibble %>%
                 select(set_bench, auc, statistic) %>%
                 rename(pr_auc = auc),
               by = c("set_bench", "statistic")) %>%
    inner_join(x=.,
               y=(roc %>%
                    group_by(name_lvl) %>%
                    summarise(source_coverage = "roc$coverage") %>%
                    distinct() %>%
                    ungroup() %>%
                    separate(col="name_lvl",
                             into=c("set_bench", "lvls"),
                             sep="\\.")),
               by = c("set_bench", "lvls")) %>%
    inner_join(x=.,
               y=comp_time,
               by = c("set_bench", "lvls")) %>%
    distinct()

  bench_summary <- list(auroc_summary, roc_plot, pr_roc_plot,
                        auroc_plot, auroc_heat, prroc_heat)

  names(bench_summary) <- c("auroc_summary", "roc_plot", "pr_roc_plot",
                            "auroc_plot", "auroc_heat", "prroc_heat")

  return(bench_summary)
}



#' Helper function to format (PR) Receiver Operator Curve results
#' @param .res_tibble formatted bench result tibble with added auroc column
#' @param roc_column (PR) ROC column to format
#' @return returns
format_roc <- function(.res_tibble, roc_column){
  apply(.res_tibble, 1, function(df) {
    df[roc_column] %>%
      as_tibble() %>%
      mutate(name = df$set_bench,
             lvls = df$lvls,
             statistic = df$statistic) %>%
      unite("name_lvl", name, lvls, remove = F, sep = ".") %>%
      unite("run_key", name, statistic, lvls, remove = F)
  }) %>%
    do.call(rbind, .)
}


#' Helper function to produce AUROC heatmap
#' @param auroc_tibble Tibble with calculated AUROC
#' @return returns an AUROC heatmap
#' @import ggplot2
get_auroc_heat <- function(auroc_tibble){
  auroc_tibble %>%
    select(statistic, auc, lvls, set_bench) %>%
    unite("name_lvl", set_bench, lvls) %>%
    pivot_wider(names_from = name_lvl, values_from = auc) %>%
    column_to_rownames(var = "statistic")  %>%
    pheatmap(.,
             cluster_rows = F,
             treeheight_col = 0,
             treeheight_row = 0,
             display_numbers = T,
             silent = T,
             cluster_cols=F)
}
