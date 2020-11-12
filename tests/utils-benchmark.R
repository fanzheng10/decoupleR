#' This function calculates precision recall curves
#'
#' @param df run_method_viper() output
#' @return tidy data frame containing recall, precision, auc, tp, tn and
#'   coverage
calc_pr_curve = function(df) {
  library(ggplot2)

  df = df %>% prepare_for_roc(., filter_tn = T)

  feature_coverage = length(unique(df$tf))

  tn = df %>% filter(response == 0)
  tp = df %>% filter(response == 1)

  r = pr.curve(scores.class0 = df$predictor,
               weights.class0 = df$response,
               curve=T)
  res = r$curve %>%
    as_tibble() %>%
    setNames(., c("recall", "precision", "th")) %>%
    mutate(auc = r$auc.davis.goadrich,
           type = r$type,
           n = sum(df$response),
           tp = nrow(tp),
           tn = nrow(tn),
           coverage = feature_coverage) %>%
    arrange(recall, desc(precision))
}


#' This function calculates receiver operating characteristic
#'
#' @param df run_method_viper() output
#' @param downsampling logical flag indicating if the number of TN should be
#'   downsampled to the number of TP
#' @param times integer showing the number of downsampling
#' @param ranked logical flag indicating if input is derived from composite
#'   ranking that already took up-/downregulation (sign) into account
#'
#' @return tidy data frame with tpr, fpr, auc, n, tp, tn and coverage
#'
#' @import yardstick
calc_roc_curve = function(df, downsampling = F, times = 1000, ranked = F) {

  library(yardstick)

  if (ranked == T) {
    df = df %>% prepare_for_roc(., filter_tn = T, ranked = T)
  } else {
    df = df %>%
      prepare_for_roc(., filter_tn = T, ranked = F)
  }

  if (length(which(df$response == 0)) == nrow(df)){
    return(as_tibble(NULL))
  }

  tn = df %>% filter(response == 0)
  tp = df %>% filter(response == 1)

  feature_coverage = length(unique(df$tf))

  if (downsampling == T) {
    # number of true positives
    num_tp = nrow(tp)

    res = map_df(seq(from=1, to=times, by=1), function(i) {
      df_sub = sample_n(tn, num_tp, replace=TRUE) %>%
        bind_rows(tp)

      r_sub = df_sub %>%
        roc_curve(response, predictor)
      auc = df_sub %>%
        roc_auc(response, predictor) %>%
        select(.estimate)

      res_sub = tibble(tpr = r_sub$sensitivities,
                       fpr = 1-r_sub$specificities,
                       th = r_sub$thresholds,
                       auc = r_sub$auc,
                       n = length(which(df$response == 1)),
                       tp = nrow(tp),
                       tn = nrow(tn),
                       coverage = feature_coverage) %>%
        mutate_("run" = i)

    })
  } else {
    r = df %>%
      roc_curve(response, predictor)
    auc = df %>%
      roc_auc(response, predictor)

    res = tibble(tpr = r$sensitivity,
                 fpr = 1-r$specificity,
                 th = r$.threshold,
                 auc = auc$.estimate,
                 n = length(which(df$response == 1)),
                 tp = nrow(tp),
                 tn = nrow(tn),
                 coverage = feature_coverage) %>%
      arrange(fpr, tpr)

  }
  return(res)
}


#' This function takes the result from viper
#' and prepares the data frame for roc/pr curve analysis. For each TF-experiment
#' combination a response and a predictor value is generated.
#'
#' @param df run_method_viper() output
#' @param filter_tn logical flag indicating if unnecessary true negatives should
#'   be filtered out (unnecessary means that there are no true positives for a
#'   given tf)
#' @param ranked logical flag indicating if input is derived from composite
#'   ranking that already took up-/downregulation (sign) into account
#'
#' @return tidy data frame with meta information for each experiment and the
#'   response and the predictor value which are required for roc curve analysis
prepare_for_roc = function(df, filter_tn = F, ranked = F) {
  res = df %>%
    dplyr::mutate(response = case_when(tf == target ~ 1,
                                       tf != target ~ 0),
                  predictor =  case_when(ranked == F ~ score*sign,
                                         ranked == T ~ score))
  res$response = factor(res$response, levels = c(1, 0))

  if (filter_tn == TRUE) {
    # Only TF which are perturbed and predicted are considered
    z = intersect(res$tf, res$target)
    res = res %>%
      filter(tf %in% z, target %in% z)
  }
  res %>%
    select(c(tf, id, response, predictor))
}




#' FUNCTION: Plot ROC curve and AUROC
#'
#' @param data
#' @param title character string for title of plots
#' @param coverage print TF coverage (TRUE/FALSE)
#' @return TF coverage (optional), ROC and AUROC plot
bench_sumplot <- function(data, title, coverage = FALSE) {

  roc <- apply(data, 1, function(df) {
    mutate(df$roc, name = df$name)
  }) %>%
    do.call(rbind, .)

  # TF Coverage
  if (coverage == TRUE) {
    cov <- roc %>%
      group_by(name) %>%
      summarise(coverage) %>%
      distinct() %>%
      arrange(-coverage)
  }

  # Plot ROC
  roc_plot <- ggplot(roc, aes(x = fpr, y = tpr, colour = name)) +
    geom_line() +
    ggtitle(paste("ROC curve:", title)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("FPR (1-specificity)") +
    ylab("TPR (sensitivity)")

  # Extract AUROC
  auroc <- data %>%
    unnest(roc) %>%
    select(name, confidence = regs, statistics, auroc = auc) %>%
    distinct()

  # Plot AUROC
  auroc_plot <- ggplot(auroc, aes(x = reorder(name, auroc), y = auroc, fill = name)) +
    geom_bar(stat = "identity") +
    ggtitle(paste("AUROC:", title)) +
    xlab("networks") +
    ylab("AUROC") +
    coord_flip(ylim = c(0.5, 0.8)) +
    theme(legend.position = "none")

  # Plot AUROC heat
  auroc_heat <- auroc %>%
    select(statistics, auroc, confidence) %>%
    pivot_wider(names_from = statistics, values_from = auroc) %>%
    column_to_rownames(var = "confidence") %>%
    pheatmap(.,
             cluster_rows = F,
             treeheight_col = 0,
             treeheight_row = 0)


  out <- list(auroc ,cov, roc_plot, auroc_plot, auroc_heat)
  names(out) <- c("auroc", "coverage", "roc_plot", "auroc_plot", "auroc_heat")

  return(out)
}

