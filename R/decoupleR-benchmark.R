#' Benchmark gene sets with decouple
#'
#' @inheritParams format_design
#' @param .minsize regulon/gene set minimum number of targets/members
#' @param .form bool whether to format or not
#' @param .perform bool whether to calculate roc and performance summary
#' @return An S4 object of class BenchResult
run_benchmark <- function(.design,
                          .minsize = 10,
                          .form = T,
                          .perform = T,
                          .silent = T,
                          .downsample_pr = F,
                          .downsample_roc = F,
                          .downsample_times = 100
                          ){
  res <- .design %>%
    format_design() %>%
    mutate(activity = pmap(.,
                           .f=function(set_name, bench_name,
                                       stats_list, opts_list,
                                       bexpr_loc, bmeta_loc, source_loc,
                                       source_col, target_col,
                                       filter_col, filter_crit,
                                       .source_bln, .expr_bln, .meta_bln){

      # Check_prereq
       if(!.expr_bln){
         .GlobalEnv$gene_expression <- readRDS(bexpr_loc) %>% as.matrix()
       }
       if(!.meta_bln){
         .GlobalEnv$meta_data <- readRDS(bmeta_loc)
       }
       if(!.source_bln){
         .GlobalEnv$set_source <- check_prereq(source_loc, target_col,
                                               source_col, filter_col)
       }

      # Filter set_source/network
      ss_filtered <- filter_sets(set_source, source_col,
                                 filter_col, filter_crit,
                                 .minsize, .silent)

      # Print Current Row/Run
      if(!.silent){
        .curr_row <- paste(set_name, bench_name,
                           paste0(unlist(filter_crit), collapse=""),
                           sep="_")
        message(str_glue("Currently Running: {.curr_row}"))
      }

      # Obtain Activity with decouple and format
      decouple(mat = gene_expression, network = ss_filtered,
               .source = source_col, .target = target_col,
               statistics = stats_list,
               .options = opts_list)  %>%
        dplyr::rename(id=condition) %>%
        inner_join(meta_data, by="id")  %>%
        group_split(statistic, .keep=T)
      })) %>% {
      if(.form & !.perform) bench_format(., silent=.silent)
      else if(.form & .perform) bench_format(., silent=.silent) %>%
        mutate(roc = activity %>%
                 map(~calc_roc_curve(df=.x,
                                    downsampling=.downsample_roc,
                                    times=.downsample_times)),
               prroc = activity %>%
                 map(~calc_pr_curve(df=.x,
                                    downsampling=.downsample_pr,
                                    times=.downsample_times)))
      else .
    }

  if(.form & .perform){
    bench_result <-new("BenchResult",
                       bench_res=res,
                       summary=res %>% get_bench_summary(),
                       design=.design)
  }
  else{
    bench_result <-new("BenchResult",
                       bench_res=res,
                       summary=list(NULL),
                       design=.design)
  }
  return(bench_result)
}
