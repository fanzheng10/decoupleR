#' eGSA
#'
#' This function is a convenient wrapper for using one of the eGSA statistics
#' in decoupleR.
#'
#' @param emat An expression matrix with genes (HGNC symbol) in rows and samples
#'  in columns.
#' @param genesets A data frame of gene sets. The structure is dependent on the
#' gene set resource.
#' @param statistic A character indicating which of the following statistics
#' should be used to calculate the activity: mean, normalized_mean, GSEA
#' @param options A list of named options to pass to
#' SAFE_mean() such as \code{minsize} or \code{mixed}.
#'
#' @return A matrix of normalized enrichment scores for each gene set across all
#'  samples.
#'

run_eGSA = function(emat, genesets, statistic, options = list()) {

  # gs_resource to list
  regulon_list = df_to_list(genesets)

  # get tf activity by using one of the eGSA statistics
  result_list = apply(emat, 2, function(sample){
    if (statistic == "mean"){
      ES = SAFE_mean(FLS_vect = sample, FSC = regulon_list)$directional %>%
        select(SetName, ES)
    } else if (statistic == "normalized_mean"){
      ES = SAFE_mean_normalized(FLS_vect = sample, FSC = regulon_list)$directional %>%
        select(SetName, ES)
    } else if (statistic == "GSEA"){
      ES = SAFE_GSEA_scores(FLS_vect = sample, FSC = regulon_list)$directional %>%
        select(SetName, ES)
    }
  })

  tf_activity <- bind_rows(result_list, .id = "sample") %>%
    data.frame() %>%
    spread(key = sample, value = ES) %>%
    rename(tf = SetName)

  print("finished eGSA")

  return(tf_activity)

}

#' Create regulon list
#'
#' This function converts the regulon data frame into a list needed for the mean
#' statistics.
#'
#' @param genesets A data frame of gene sets including the following columns:
#' tf, confidence, target, mor and likelihood.
#'
#' @return A list of the regulons that can be used as FSC in SAFE_mean_scores().
#'
df_to_list = function(genesets){
  tfs = as.list(unique(genesets$tf))
  regulon_list = sapply(tfs, function(x) NULL)
  names(regulon_list) = tfs

  for(regulon_name in tfs){
    regulon_list[[regulon_name]][["MEMBERS"]] = filter(genesets, tf == regulon_name)$target
    regulon_list[[regulon_name]][["CONTRIBUTION"]] = filter(genesets, tf == regulon_name)$likelihood
    regulon_list[[regulon_name]][["MOR"]] = filter(genesets, tf == regulon_name)$mor
    regulon_list[[regulon_name]][["CONFIDENCE"]] = filter(genesets, tf == regulon_name)$confidence
  }
  return(regulon_list)
}



