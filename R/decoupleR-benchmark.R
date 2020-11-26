#' Benchmark gene sets with decouple
#'
#' @inheritParams format_design
#' @param .minsize regulon/gene set minimum number of targets/members
#' @param .lvls column name in gene set/network resource that you wish to filter
#' to, in relation to the lvls in design
#' @param .form bool whether to format or not
#' @param .perform bool whether to calculate roc and performance summary
#' @return An S4 object of class BenchResult
run_benchmark <- function(design,
                          .minsize = 10,
                          .lvls = "confidence",
                          .form = T,
                          .perform = T
                          ){
  .lvls <- ensym(.lvls)

  res <- design %>%
    format_design() %>%
    mutate(activity = pmap(.,
                           .f=function(row_name, net_loc, lvls,
                                       gene_source, target, statistics,
                                       bnch_expr, bench_meta,
                                       .net_bln, .expr_bln, .meta_bln, opts){

      # Check conditions and load prerequisites
      if(!.net_bln){
        .GlobalEnv$network <- readRDS(net_loc)
      }
      if(!.expr_bln){
        .GlobalEnv$gene_expression <- readRDS(bnch_expr) %>% as.matrix()
      }
      if(!.meta_bln){
        .GlobalEnv$meta_data <- readRDS(bench_meta)
      }

      # filter network (to be changed and extended for additional filters)
      network_filtered <- network %>%
        dplyr::filter((!!.lvls) %in% lvls) %>%
        distinct_at(vars(-(!!.lvls)), .keep_all = T) %>%
        rename(.source = ensym(gene_source)) %>% #*
        group_by(.source) %>%
        add_count() %>%
        filter(n >= .minsize) %>%
        ungroup() %>%
        rename(gene_source = .source) #* !!ensym(gene_source) not found

      # Print to track libs
      print(paste(row_name, paste0(unlist(lvls), collapse=""), sep="_"))

      # Obtain Activity with decouple and format
      decouple(mat = gene_expression, network = network_filtered,
               .source = gene_source, .target = target,
               statistics = statistics,
               .options = opts)  %>%
        dplyr::rename(id=condition) %>%
        inner_join(meta_data, by="id")  %>%
        group_split(statistic, .keep=T)
    })) %>% {
      if(.form & !.perform) bench_format(.)
      else if(.form & .perform) bench_format(.) %>%
        mutate(roc = activity %>% map(calc_roc_curve))
      else .
    }

  # handle return
  if(.form & .perform){
    bench_result <-new("BenchResult",
                       bench_res=res,
                       summary=res %>% bench_sumplot(),
                       design=design_d)
  }
  else{
    bench_result <-new("BenchResult",
                       bench_res=res,
                       summary=list(NULL),
                       design=design_d)
  }
  return(bench_result)
}


#' Helper Function to to generate the booleans used to check if the current
#' locations/data objects are the same as the previous one
#' @param design location of a json file with run design specifications
#' @return a design tibble to be used in benchmarking
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


#' Helper function to format benchmarking results
#'
#' @param bench_res benchmarking results
#' @returns formatted benchmarking results
bench_format <- function(bench_res){
  res_format <- bench_res %>%
    unnest(activity) %>%
    # get statistic time from activity
    mutate(statistic_time = activity %>%
             map(function(tib)
               tib %>%
                 select(statistic_time) %>%
                 unique)) %>%
    unnest(statistic_time) %>%
    # calculate regulon size
    group_by(row_name, lvls) %>%
    mutate(regulon_time = sum(statistic_time)) %>%
    # convert lvls from character to string
    rowwise() %>%
    mutate(lvls = paste0(unlist(lvls), collapse = "")) %>%
    ungroup() %>%
    # get statistic row_name
    mutate(statistic = activity %>%
             map(function(tib)
               unique(tib[["statistic"]]))) %>%
    unnest(statistic) %>%
    select(row_name, lvls, statistic, statistic_time, regulon_time, activity)
  return(res_format)
}



#' Function that provides summary
#'
#' @param .res_tible formatted bench result tible with added roc column
#' @param title character string for title of plots
#' @return AUROC summary per row with TF coverage, ROC AUROC, Heatmap plots
#' @import ggplot2, pheatmap
bench_sumplot <- function(.res_tible, title = "") {

  .res_tible <- .res_tible %>%
    select(row_name, lvls, statistic, roc)

  roc <- apply(.res_tible, 1, function(df) {
    df$roc %>%
      mutate(name = df$row_name,
             lvls = df$lvls,
             statistic = df$statistic) %>%
      unite("name_lvl", name, lvls, remove = F, sep = ".") %>%
      unite("run_key", name, statistic, lvls, remove = F)
  }) %>%
    do.call(rbind, .)

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
    select(row_name, lvls, statistic, auc) %>%
    distinct() %>%
    # join coverage
    inner_join(x=.,
               y=(roc %>%
                    group_by(name_lvl) %>%
                    summarise(cov = coverage) %>%
                    distinct() %>%
                    ungroup() %>%
                    separate(col="name_lvl",
                             into=c("row_name", "lvls"),
                             sep="\\.")))


  # Plot AUROC
  auroc_plot <- auroc %>%
    unite("run_key", row_name, statistic, lvls, remove = F) %>%
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
    select(statistic, auc, lvls, row_name) %>%
    unite("name_lvl", row_name, lvls) %>%
    pivot_wider(names_from = name_lvl, values_from = auc) %>%
    column_to_rownames(var = "statistic") %>%
    pheatmap(.,
             cluster_rows = F,
             treeheight_col = 0,
             treeheight_row = 0,
             display_numbers = T,
             silent = T)


  bench_summary <- list(auroc, roc_plot, auroc_plot, auroc_heat)
  names(bench_summary) <- c("auroc_summary", "roc_plot", "auroc_plot", "auroc_heat")

  return(bench_summary)
}
