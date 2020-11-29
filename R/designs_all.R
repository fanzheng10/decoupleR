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
#' @usage readRDS(system.file("benchdata/inputs", "design_all.rds" , package = "decoupleR"))
#'
#' \describe{
#' \item{row_name}{user-defined name of row/run - e.g. network + benchmark data}
#' \item{net_loc}{location of the network/gene set resource}
#' \item{lvls}{column to filter by - e.g. confidence levels of the network}
#' \item{gene_source}{source of nodes for the network AKA .source}
#' \item{target}{network nodes AKA .target}
#' \item{statistics}{list of statistics to run}
#' \item{bnch_expr}{location of the benchmark data set (in the form of a matrix
#'    with each experiment summarized in a column}
#' \item{bench_meta}{location of benchmark metadata (information about the
#'    different experiments - e.g. type of method, cell type, etc.)}
#' \item{opts}{named list to be used to pass the options for each statistic}
#' }
#'
#' @importFrom tibble tibble
#' @details Read and executed in a row-wise manner. Say that this can also
#' be used with different type of gene set resources (e.g. KEGG, GO, etc)
NULL
