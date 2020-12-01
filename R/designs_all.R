#' Design tibble with all networks, stats, and benchmarks.
#'
#' @description To be used as default input or as a starting point and
#'  modified according to the benchmark run in mind.
#'
#' @name design_all
#'
#' @docType data
#'
#' @format A tibble with locations, options, and filter options for the desired benchmark setting
#'
#' \describe{
#' \item{row_name}{User-defined name of row/run - e.g. network + benchmark data}
#' \item{net_loc}{Location of the network/gene set resource}
#' \item{lvls}{Column to filter by - e.g. confidence levels of the network}
#' \item{gene_source}{Source of nodes for the network AKA .source}
#' \item{target}{Network nodes AKA .target}
#' \item{statistics}{List of statistics to run}
#' \item{bnch_expr}{Location of the benchmark data set. The data set is in the form of a matrix
#'    with each experiment summarized in a column}
#' \item{bench_meta}{Location of benchmark metadata. The metadata contains information about the
#'    different experiments - e.g. type of method, cell type, etc.}
#' \item{opts}{Named list to be used to pass the options for each statistical method}
#' }
#'
#' @importFrom tibble tibble
#' @details Read and executed in a row-wise manner.
NULL
