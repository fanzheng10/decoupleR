#' Benchmark gene sets with decouple
#'
#' @inheritParams format_design
#' @param minsize regulon/gene set minimum number of targets/members
#' @param format bool whether to format or not
#' @return A tibble with an appended activity column that corresponds
#' to the activities calculated for each row of the design tibble
run_benchmark <- function(design,
                          .minsize = 10,
                          .format = T){

  res <- design %>%
    format_design() %>%
    mutate(activity = pmap(., function(name, net_loc, regs,
                                       gene_source, target, statistics,
                                       bnch_expr, bench_meta,
                                       net_bln, expr_bln, meta_bln, opts){

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
        filter(n >= .minsize) %>%
        ungroup()

      # Print to track libs
      print(paste(name, paste0(unlist(regs), collapse=""), sep="_"))

      # Obtain Activity with decouple and format
      decouple(mat = gene_expression, network = network_filtered,
               .source = all_of(gene_source), .target = target,
               statistics = statistics,
               .options = opts)  %>%
        dplyr::rename(id=condition) %>%
        inner_join(meta_data, by="id")  %>%
        group_split(statistic, .keep=T)
    })) %>% {
      if(.format) bench_format(.) else .
    }

  return(res)
}


#' Load Design from JSON file and Format
#'
#' @param design location of a json file with run design specifications
#' @return a design tibble to be used in benchmarking
format_design <- function(design){
  design %>%
    mutate(net_bln = net_loc %>% check_prereq(),
           expr_bln = bnch_expr %>% check_prereq(),
           meta_bln = bench_meta %>% check_prereq())
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
    mutate(statistic_time = activity %>% map(function(tib)
      tib %>%
        select(statistic_time) %>%
        unique)) %>%
    unnest(statistic_time) %>%
    rowwise() %>%
    mutate(regs = paste0(unlist(regs), collapse = "")) %>%
    group_by(name, regs) %>%
    mutate(regulon_time = sum(statistic_time)) %>%
    select(name, regs, statistics, activity, statistic_time, regulon_time)
  return(res_format)
}
