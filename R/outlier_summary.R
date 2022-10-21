outlier_summary <- function(actual, flagged, roundto = 3, n_out = NULL){
  n_out <- length(actual)
  n_flagged <- length(flagged)
  n_correct <- length(intersect(actual, flagged))
  precision <- n_correct / length(flagged)# precision
  recall<- n_correct / n_out # recall
  if(roundto>0){
    c("n_flagged" = n_flagged, "n_correct" = n_correct, round(c("precision" = precision, "recall" = recall, "F-score" = 2*(precision*recall)/(precision + recall)),roundto))
  } else {
    c("n_flagged" = n_flagged, "n_correct" = n_correct, "precision" = precision, "recall" = recall, "F-score" = 2*(precision*recall)/(precision + recall))
  }
}
