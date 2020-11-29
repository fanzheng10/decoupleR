#' S4 Class used to format benchmark output.
#'
#' @slot bench_res Formatted or non-formatted Benchmark output
#' @slot summary Summary returned by the bench_sumplot functions - it countains
#' auroc and coverage summary table, roc plot, and auroc heatmap and barplot
#' @slot design The input design tibble used to generate the benchmark results
#'
#' @name BenchResult-class
#' @rdname BenchResult-class
#' @export
setClass("BenchResult",
         slots=list(bench_res="tbl_df",
                    summary="list",
                    design="tbl_df"))
