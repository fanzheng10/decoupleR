#' S4 Class used to format benchmark wrapper results.
#'
#' @slot bench_res Formatted or non-formatted Benchmark results
#' @slot summary Summary return by the bench_sumplot function
#' @slot design The input design tibble used to generate the results
#'
#' @name BenchResult-class
#' @rdname BenchResult-class
#' @export
setClass("BenchResult",
         slots=list(bench_res="tbl_df",
                    summary="list",
                    design="tbl_df"))
