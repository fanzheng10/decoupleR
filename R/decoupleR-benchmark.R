#' Load Design from JSON file
#' @param design_loc
#' @return a design tibble to be used in benchmarking
load_design <- function(design_loc){
  library(jsonlite)
  # load JSON
  design <- fromJSON(design_loc) %>%
    as_tibble()

  # check prerequisites
  design <- design %>%
    mutate(net_bln = net_loc %>% check_prereq(),
           expr_bln = bnch_expr %>% check_prereq(),
           meta_bln = bench_meta %>% check_prereq())
  return(design)
}

#' Function that checks ifthe  preceding vector element is the same
#' as the current element
#' @param vector_loc char vector with directories
#' @return logical value
check_prereq <- function(vector_loc){
  tib_loc <- tibble(current=vector_loc, behind=lag(vector_loc))

  pmap_lgl(tib_loc, function(behind, current){
    ifelse(is.na(behind) || behind!=current, FALSE, TRUE)
  })
}


#' Benchmark gene sets with decouple
#' @param design_tibble tibble with design specifications
#' @param opts options for each stats method. Note that this list should match
#' the statistics passed to .f decouple
#' @retun A tibble with an appended activity column that corresponds
#' to the activities calculated for each row of the design tibble
bench_couple <- function(design_tibble, opts){
  res <- design_tibble %>%
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

      # filter network
      network_filtered <- network %>%
        dplyr::filter(confidence %in% regs) %>%
        distinct()

      # Obtain Activity with decouple and format
      decouple(mat = gene_expression, network = network_filtered,
               .source = all_of(gene_source), .target = target,
               statistics = statistics,
               .options = opts) %>%
        dplyr::rename(id=condition) %>%
        inner_join(meta_data, by="id") %>%
        dplyr::select(-c(.data$run_id, .data$p_value)) %>%
        group_split(statistic, .keep=T) %>%
        as.list()
    }))

  return(res)
}


#' Function to format benchmarking results
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
    # rowwise() %>%
    # unite(name, regs, sep="_", col="name") %>%
    # ungroup()  %>%
    select(name, regs, statistics, activity)
  return(res_format)
}
