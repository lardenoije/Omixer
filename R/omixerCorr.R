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
#' @importFrom stats chisq.test
#' @importFrom stats cor.test
#' @importFrom forcats as_factor
#' @importFrom tibble tibble
#' @export

omixerCorr <- function(x, y) {

    ## Convert variables to numeric
    if(class(x) %in% c("character")) x <- as.numeric(as_factor(x))
    if(class(y) %in% c("character")) y <- as.numeric(as_factor(y))
    if(class(x) %in% c("Date", "factor")) x <- as.numeric(x)
    if(class(y) %in% c("Date", "factor")) y <- as.numeric(y)

    ## Save correlation estimates and p values
    if(length(unique(x)) < 5 & length(unique(y)) < 5) {
        ## Two categorical variables
        cvStat <- chisq.test(x, y, correct=FALSE)$statistics
        if(is.null(cvStat)) cvStat <- 0
        cv <- sqrt(cvStat/(length(x)*min(length(unique(x)),
            length(unique(y)))-1))
        corVal <- as.numeric(cv)
        corP <- chisq.test(x, y)$p.value
    } else {
        ## Otherwise, use Kendall's correlation coefficient and p-value
        corVal <- cor.test(x, y, method="kendall")$estimate
        corP <- cor.test(x, y, method="kendall")$p.value
    }

    ## Return as a data frame of estimate and p-value
    corTb <- tibble("corVal"=corVal, "corP"=corP)

  return(corTb)
}
