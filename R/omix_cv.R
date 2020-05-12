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
#' @examples
#' sex <- sample(1:2, 20, replace = TRUE)
#' plate <- rep(1:4, each = 5)
#' omix_corr(sex, plate)
#'
#' @importFrom stats chisq.test
#' @export

omix_cv <- function(x, y) {

  # Calculate Cramer's V
  cv <- sqrt(chisq.test(x, y, correct=FALSE)$statistics /
               (length(x) *
                  min(length(unique(x)),
                      length(unique(y))) - 1))
  return(as.numeric(cv))
}
