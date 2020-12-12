# Helper functions for run_benchmark ----
#' Helper Function to to generate the booleans used to check if the current
#' locations/data objects are the same as the previous one
#' @param design location of a json file with run design specifications
#' @return a design tibble to be used in benchmarking
#' @export
format_design <- function(design){
  design %>%
    mutate(.source_bln = source_loc %>% check_prereq(),
           .expr_bln = bexpr_loc %>% check_prereq(),
           .meta_bln = bmeta_loc %>% check_prereq())
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
filter_sets <- function(set_source, source_col, .lvls, lvls, .minsize, silent){

  n_duprows <- sum(duplicated(set_source))

  gs_filtered <- set_source %>%
    dplyr::filter(.data[[.lvls]] %in% lvls) %>%
    distinct_at(vars(-.data[[.lvls]]), .keep_all = F) %>%
    rename(.source = source_col) %>% #*
    group_by(.source) %>%
    add_count() %>%
    filter(n >= .minsize) %>%
    ungroup() %>%
    rename(source_col = .source) #* !!ensym(source_col) not found

  if (n_duprows & !silent){
    warning(str_glue("{n_duprows} rows were duplicated in the set resource! ",
                     "{sum(duplicated(gs_filtered))} duplicated rows ",
                     "remain after filtering."))
  }
  return(gs_filtered)
}
