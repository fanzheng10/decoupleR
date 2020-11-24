#' Benchmark gene sets with decouple
#'
#' @inheritParams format_design
#' @param minsize regulon/gene set minimum number of targets/members
#' @param format bool whether to format or not
#' @param .confidence confidence levels/regulon names,etc
#' @return A tibble with an appended activity column that corresponds
#' to the activities calculated for each row of the design tibble
run_benchmark <- function(design,
                          .minsize = 10,
                          .form = T,
                          .confidence = "confidence" # enable a way to be NA
                          ){

  .confidence <- ensym(.confidence)

  design %>%
    format_design() %>%
    mutate(activity = pmap(. , function(name, net_loc, regs, gene_source,
                                        target, statistics, bnch_expr,
                                        bench_meta, net_bln, expr_bln,
                                        meta_bln, opts){

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
        dplyr::filter((!!.confidence) %in% regs) %>%
        distinct_at(vars(-(!!.confidence)), .keep_all = T) %>%
        rename(.source = ensym(gene_source)) %>% #*
        group_by(.source) %>%
        add_count() %>%
        filter(n >= .minsize) %>%
        ungroup() %>%
        rename(gene_source = .source) #* !!ensym(gene_source) not found

      # Print to track libs
      print(paste(name, paste0(unlist(regs), collapse=""), sep="_"))

      # Obtain Activity with decouple and format
      decouple(mat = gene_expression, network = network_filtered,
               .source = gene_source, .target = target,
               statistics = statistics,
               .options = opts)  %>%
        dplyr::rename(id=condition) %>%
        inner_join(meta_data, by="id")  %>%
        group_split(statistic, .keep=T)
    })) %>% {
      if(.form) bench_format(.) else .
    }
}


#' Helper Function to to generate the booleans used to check if the current
#' locations/data objects are the same as the previous one
#' @param design location of a json file with run design specifications
#' @return a design tibble to be used in benchmarking
format_design <- function(design){
  design %>%
    mutate(net_bln = net_loc %>% check_prereq(),
           expr_bln = bnch_expr %>% check_prereq(),
           meta_bln = bench_meta %>% check_prereq())
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
    group_by(name, regs) %>%
    mutate(regulon_time = sum(statistic_time)) %>%
    # convert regs from character to string
    rowwise() %>%
    mutate(regs = paste0(unlist(regs), collapse = "")) %>%
    ungroup() %>%
    # get statistic name
    mutate(statistic = activity %>%
             map(function(tib)
               unique(tib[["statistic"]]))) %>%
    unnest(statistic) %>%
    select(name, regs, statistic, statistic_time, regulon_time, activity)
  return(res_format)
}
