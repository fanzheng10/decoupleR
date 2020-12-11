# Helper functions for run_benchmark ----
#' Helper Function to to generate the booleans used to check if the current
#' locations/data objects are the same as the previous one
#' @param design location of a json file with run design specifications
#' @return a design tibble to be used in benchmarking
#' @export
format_design <- function(design){
  design %>%
    mutate(.net_bln = net_loc %>% check_prereq(),
           .expr_bln = bnch_expr %>% check_prereq(),
           .meta_bln = bench_meta %>% check_prereq())
}


#' Helper Function that checks if the  preceding vector element is the same
#' as the current element
#'
#' @param vector_loc char vector with directories
#' @return logical values describing whether the location of the loaded files
#' has changes
#' @export
check_prereq <- function(vector_loc){
  tib_loc <- tibble(current=vector_loc, behind=lag(vector_loc))

  pmap_lgl(tib_loc, function(behind, current){
    ifelse(is.na(behind) || behind!=current, FALSE, TRUE)
  })
}


#' Helper Function to filter and format the gene set resource
#'
#'
#'
#'
#'
filter_sets <- function(set_source, gene_source, .lvls, lvls, .minsize){

  n_duprows <- sum(duplicated(set_source))

  gs_filtered <- set_source %>%
    dplyr::filter(.data[[.lvls]] %in% lvls) %>%
    distinct_at(vars(-.data[[.lvls]]), .keep_all = F) %>%
    rename(.source = gene_source) %>% #*
    group_by(.source) %>%
    add_count() %>%
    filter(n >= .minsize) %>%
    ungroup() %>%
    rename(gene_source = .source) #* !!ensym(gene_source) not found

  if (n_duprows){
    warning(str_glue("{n_duprows} rows were duplicated in the set resource! ",
                     "{sum(duplicated(gs_filtered))} duplicated rows ",
                     "remain after filtering."))
  }
  return(gs_filtered)
}



# Downstream reformatting ----
#' Function to format benchmarking results
#'
#' @param bench_res benchmarking results
#' @returns formatted benchmarking results
#' @export
bench_format <- function(bench_res){
  res_format <- bench_res %>%
    unite("set_bench", set_name, bench_name) %>%
    unnest(activity) %>%
    # convert lvls from character to string
    rowwise() %>%
    mutate(lvls = paste0(unlist(lvls), collapse = "")) %>%
    ungroup() %>%
    # get statistic name
    mutate(statistic = activity %>%
             map(function(tib)
               unique(tib[["statistic"]]))) %>%
    unnest(statistic) %>%
    select(set_bench, lvls, statistic, activity)

  # Check and filter infinite values
  inf_sums <- lapply(res_format$activity,
                     function(x) sum(is.infinite(x$score))) %>%
    setNames(paste(res_format$set_bench, res_format$statistic, sep="_")) %>%
    enframe() %>% unnest(value)

  if(sum(inf_sums$value)){
    res_format <- res_format %>%
      mutate(activity = activity %>%
               map(function(tib) tib %>%
                     mutate_at(vars(score), ~replace(., is.infinite(.), 0))
               ))

    warning(inf_sums %>%
              filter(value > 0) %>%
              str_glue_data("{.$name} had {.$value} infinite values. \n "))

  }
  return(res_format)
}



#' Function that provides summary and plots for the benchmark run
#'
#' @param .res_tible formatted bench result tible with added roc column
#' @param title character string for title of plots
#' @return AUROC summary per row with TF coverage, ROC AUROC, Heatmap plots
#' @import ggplot2
#' @import pheatmap
bench_sumplot <- function(.res_tible, title = "") {

  # get roc results
  roc <- apply(.res_tible, 1, function(df) {
    df$roc %>%
      mutate(name = df$set_bench,
             lvls = df$lvls,
             statistic = df$statistic) %>%
      unite("name_lvl", name, lvls, remove = F, sep = ".") %>%
      unite("run_key", name, statistic, lvls, remove = F)
  }) %>%
    do.call(rbind, .)

  # get computational time info
  comp_time <- .res_tible %>%
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

  # Plot ROC
  roc_plot <- ggplot(roc, aes(x = fpr, y = tpr, colour = run_key)) +
    geom_line() +
    ggtitle(paste("ROC curve:", title)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("FPR (1-specificity)") +
    ylab("TPR (sensitivity)")

  # Extract AUROC
  auroc <- .res_tible %>%
    unnest(roc) %>%
    select(set_bench, lvls, statistic, auc) %>%
    distinct()

  # Plot AUROC
  auroc_plot <- auroc %>%
    unite("run_key", set_bench, statistic, lvls, remove = F) %>%
    ggplot(., aes(x = reorder(run_key, auc),
                  y = auc,
                  fill = run_key)) +
    geom_bar(stat = "identity") +
    ggtitle(paste("AUROC:", title)) +
    xlab("networks") +
    ylab("AUROC") +
    coord_flip(ylim = c(0.5, 0.8)) +
    theme(legend.position = "none")


  # Plot AUROC heat
  auroc_heat <- auroc %>%
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


  # Join coverage
  auroc_summary <- auroc %>%
  inner_join(x=.,
             y=(roc %>%
                  group_by(name_lvl) %>%
                  summarise(exp_coverage = coverage) %>%
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

  bench_summary <- list(auroc_summary, roc_plot,
                        auroc_plot, auroc_heat)
  names(bench_summary) <- c("auroc_summary", "roc_plot",
                            "auroc_plot", "auroc_heat")

  return(bench_summary)
}
