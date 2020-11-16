#' Benchmark gene sets with decouple
#'
#' @inheritParams load_design
#' @param save_loc location in which each row is saved
#' @param opts options for each stats method. Note that this list should match
#' the statistics passed to .f decouple
#' @return A tibble with an appended activity column that corresponds
#' to the activities calculated for each row of the design tibble
run_benchmark <- function(design_loc, opts){

  calltime <- Sys.time()

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

  res <- res %>% mutate(ctime = calltime)

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
