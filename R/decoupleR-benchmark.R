#' Benchmark gene sets with decouple
#'
#' @inheritParams load_design
#' @param save_loc location in which each row is saved
#' @param opts options for each stats method. Note that this list should match
#' the statistics passed to .f decouple
#' @return A tibble with an appended activity column that corresponds
#' to the activities calculated for each row of the design tibble
run_benchmark <- function(design_loc, opts){

  # call time of pipeline (i.e. starting time point)
  .GlobalEnv$ctime_value <- Sys.time()

  res <- load_design(design_loc) %>%
    mutate(activity = pmap(., function(name, net_loc, regs,
                                       gene_source, target, statistics,
                                       bnch_expr, bench_meta,
                                       net_bln, expr_bln, meta_bln){

      # Check conditions and load prerequisites
      if(!net_bln){
        .GlobalEnv$network <- readRDS(net_loc)
      }
      if(!expr_bln){
        .GlobalEnv$gene_expression <- readRDS(bnch_expr) %>% as.matrix()
      }
      if(!meta_bln){
        .GlobalEnv$meta_data <- readRDS(bench_meta)
      }

      # filter network (to be changed and extended for additional filters)
      network_filtered <- network %>%
        dplyr::filter(confidence %in% regs) %>%
        distinct_at(vars(-confidence), .keep_all = T) %>%
        group_by(tf) %>%
        add_count() %>%
        filter(n >= 10) %>%
        ungroup()

      # Print to track libs
      print(paste(name, paste0(unlist(regs), collapse=""), sep="_"))

      # Obtain Activity with decouple and format
      row <- decouple(mat = gene_expression, network = network_filtered,
               .source = all_of(gene_source), .target = target,
               statistics = statistics,
               .options = opts) %>%
        dplyr::rename(id=condition) %>%
        inner_join(meta_data, by="id") %>%
        dplyr::select(-c(.data$run_id, .data$p_value)) %>%
        mutate(rtime=Sys.time()) %>%
        group_split(statistic, .keep=T) %>%
        as.list()

      return(row)
    })) # %>%
    # bench_format()

  return(res)
}


#' Load Design from JSON file and Format
#'
#' @param design_loc location of a json file with run design specifications
#' @return a design tibble to be used in benchmarking
load_design <- function(design_loc){
  # load RDS
  design <- readRDS(design_loc) %>%
    mutate(net_bln = net_loc %>% check_prereq(),
           expr_bln = bnch_expr %>% check_prereq(),
           meta_bln = bench_meta %>% check_prereq())
  return(design)
}


#' Function that checks if the  preceding vector element is the same
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


#' Function to format benchmarking results
#'
#' @param bench_res benchmarking results
#' @returns formatted benchmarking results
bench_format <- function(bench_res){
  res_format <- bench_res %>%
    rowwise() %>%
    dplyr::mutate(statistics =
                    list(flatten_chr(.$activity[[1]] %>%
                                       map(function(tib)
                                         unique((tib)$statistic))))) %>%
    unnest(c(activity, statistics)) %>%
    select(name, regs, statistics, activity)
  return(res_format)
}



#' Function to calculate runtime for each statistic and gene set
#'
#' @param format_tibble a tibble resulting from a run that has been formatted
#' with bench_format /codehere
#' @return returns a tibble with runtime for each statistic and gene set
#' @details rtime - time point at the end of benchmark run for each gene set
#'  stime - '' for each statistic
#'  Runtime formula:
#'  runtime_stat1 = stime1 (end time point for stat1) - ctime (call time)
#'  runtime_stat2 = stat2 - stat1
#'  runtime_stat3 = stat3 - stat1
#'  ...
bench_runtime <- function(format_tibble){

  run_tibble <- format_tibble %>%
    mutate(stime = activity %>% map(function(tib) unique(tib$stime))) %>%
    mutate(rtime = activity %>% map(function(tib) unique(tib$rtime))) %>%
    unnest(c(stime, rtime)) %>%
    arrange(stime) %>%
    mutate(stime_lag = lag(stime),
           rtime_lag = lag(rtime)) %>%
    mutate(stat_runtime = map2(stime_lag,
                               stime,
                               get_runtime)) %>%
    mutate(reg_runtime = map2(rtime_lag,
                              rtime,
                              get_runtime)) %>%
    unnest(c(stat_runtime, reg_runtime)) %>%
    select(-c(stime, rtime, stime_lag, rtime_lag))

  return(run_tibble)
}

#' Helper function to calculate runtime
#'
#' @param start_time lagged time point
#' @param end_time time point of the end of execution for a given statistic
#' or set of regulons
get_runtime = function(start_time, end_time){
  if(is.na(start_time)){
    difftime(end_time, ctime_value)
  }else{
    difftime(end_time, start_time)
  }
}
