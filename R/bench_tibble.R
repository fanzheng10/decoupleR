#' Benchmark input tibble containing the (experimental) design for each of the
#' benchmark runs corresponding to rows.
#'
#' @description To be used as default input or as a starting point and
#'  modified according to the benchmark run in mind.
#'
#' @name bench_tibble
#'
#' @docType data
#'
#' @format A tibble with locations, options, and filter options for the desired benchmark setting
#'
#' \describe{
#' \item{set_name}{user-defined name of the set resource}
#' \item{bench_name}{user-defined name of the benchmark data}
#' \item{stats_list}{List of statistics to run}
#' \item{opts_list}{Named list containing the options for each stat. method}
#' \item{bexpr_loc}{benchmark expression data location (.rds format tibble)}
#' \item{bmeta_loc}{benchmark metadata location (.rds format tibble)}
#' \item{source_loc}{set source (e.g. network resource, gene ontology sets,
#'  kinase sets, etc.) location (.rds format tibble)}
#' \item{source_col}{name of the column with the source for the set source}
#' \item{target_col}{name of the column with the targets for the set source}
#' \item{filter_col}{name of the column by which we wish to filter}
#' \item{filter_crit}{criteria by which we wish to filter the filter column}
#'
#' @importFrom tibble tibble
#' @details Read and executed in a row-wise manner.
NULL
