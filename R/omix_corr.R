#' Correlation Tests
#'
#' This function uses appropriate tests of correlation between
#' two variables and stores the estimate and p-value in a list.
#'
#' For two categorical variables, the Cramer's V estimate is
#' stored alongside chi-square p-value. For all other combinations
#' of variables, Pearson's correlation coefficient and p-value
#' are stored.
#'
#' Please note: variables will be converted to numeric class
#' within this function.
#'
#' @param x Randomization variable (e.g. age)
#' @param y Technical covariate (e.g. plate number)
#'
#' @return List of correlation estimate and p-value
#'
#' @examples
#' age <- rnorm(n = 20, mean = 70, sd = 3)
#' plate <- rep(1:4, each = 5)
#' omix_corr(age, plate)
#'
#' @importFrom stats chisq.test
#' @importFrom stats cor.test
#' @export

omix_corr <- function(x, y) {

  # Convert variables to numeric
  if(class(x) %in% c("factor", "character", "Date")){
    x <- as.numeric(factor(x))
  }
  if(class(y) %in% c("factor", "character", "Date")){
    y <- as.numeric(factor(y))
  }

  # Save correlation estimates and p values
  if(length(unique(x)) < 5 & length(unique(y)) < 5){
    # For two categorical variables, save Cramer's V correlation estimate and chi2 p-value
    corr_val <- omix_cv(x, y)
    corr_p <- chisq.test(x, y)$p.value
  } else {
    # Otherwise, use Pearson's correlation coefficient and p-value
    corr_val <- cor.test(x, y)$estimate
    corr_p <- cor.test(x, y)$p.value
  }

  # Return as a list of estimate and p-value
  corr_list <- list("corr_val" = corr_val,
                    "corr_p" = corr_p)
  return(corr_list)
}
