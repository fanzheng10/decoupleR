
###############################################################################.
# 2. Convert to Function -------------------------------------------------------

#' Benchmark decoupleR results Function
#' @param design_json json containing an experimental design matrix
#' @param
#'
#' @return A list with gene set coverage, roc plot, and auroc plot
debench <- function()

  design_tibble_test <- fromJSON(here("inst/testdata/inputs/",
                                      "design_tibble.json"))
  design_tibble_test


  bench_dt_test <- design_tibble %>%
    mutate(activity = pmap(., function(net_name, confidence_levels,
                                       gene_source, target, statistics){
      # filter network
      network <- dorothea_genesets %>%
        filter(confidence %in% confidence_levels)

      # calculate activity
      decouple(mat = knock_ex, network = network,
               .source = all_of(gene_source), .target = target,
               statistics = statistics,
               .options = list(
                 scira = list(),
                 pscira = list(),
                 mean = list(.likelihood = NULL),
                 viper = list(options = list(verbose = FALSE, minsize=0)),
                 gsva = list(options = list(verbose = FALSE)))) %>%
        dplyr::select(-c(.data$run_id, .data$p_value)) %>%
        dplyr::rename(id=condition)
    }))

  bench_dt_test <- bench_dt_test  %>%
    mutate(statistics = list(levels(as.factor(activity[[1]]$statistic))))  %>%
    mutate(activity = activity %>%
             map(function(tib) tib %>%
                   inner_join(meta_expr, by="id") %>%
                   group_split(statistic, .keep=FALSE) %>%
                   as.list()
             ))  %>%
    unnest(c(activity, statistics)) %>%
    rowwise() %>%
    mutate(name = paste(net_name, statistics,
                        paste0(confidence_levels, collapse = ""),
                        sep = "_")) %>%
    ungroup()

  bench_dt_test <- bench_dt_test %>%
    mutate(roc = activity %>% map(calc_roc_curve))


  bench_dt_plots <- plot_roc_auroc(bench_roc, title = "Knock_TF", coverage = TRUE)


