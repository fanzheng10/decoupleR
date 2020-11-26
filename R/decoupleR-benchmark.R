#' Benchmark gene sets with decouple
#'
#' @inheritParams format_design
#' @param .minsize regulon/gene set minimum number of targets/members
#' @param .lvls column name in gene set/network resource that you wish to filter
#' to, in relation to the lvls in design
#' @param .form bool whether to format or not
#' @param .perform bool whether to calculate roc and performance summary
#' @return An S4 object of class BenchResult
#' @import ggplot2, pheatmap
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
