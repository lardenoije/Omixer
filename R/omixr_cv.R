#' Cramer's V Test
#'
#' This function calculates the Cramer's V estimate of correlation
#' between two categorical variables with a value between 0 and 1.
#'
#' @param x Randomization variable (e.g. sex)
#' @param y Technical covariate (e.g. plate number)
#'
#' @return Cramer's V estimate of correlation
#'
#' @importFrom stats chisq.test
#'
#' @keywords internal
#' @noRd

omixr_cv <- function(x, y) {

  # Calculate Cramer's V
  cv_stat <- chisq.test(x, y, correct = FALSE)$statistics

  if(is.null(cv_stat)) { cv_stat <- 0 }

  cv <- sqrt(cv_stat /
               (length(x) *
                  min(length(unique(x)),
                      length(unique(y))) - 1))

  return(as.numeric(cv))
}
