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
#' @import tibble
#' @import forcats
#' @import stringr
#' @importFrom stats chisq.test
#' @importFrom stats cor.test
#' @importFrom forcats as_factor
#' @importFrom tibble tibble
#' @export
#' 
#' @examples 
#' library(tibble)
#' library(forcats)
#' library(stringr)
#' 
#' sampleList <- tibble(sampleId=str_pad(1:48, 4, pad="0"),
#' sex=as_factor(sample(c("m", "f"), 48, replace=TRUE)), 
#' age=round(rnorm(48, mean=30, sd=8), 0), 
#' smoke=as_factor(sample(c("yes", "ex", "never"), 48, replace=TRUE)),
#' date=sample(seq(as.Date('2008/01/01'), as.Date('2016/01/01'), 
#'     by="day"), 48))
#'                 
#' omixerCorr(sampleList$age, sampleList$sex)

omixerCorr <- function(x, y) {

    ## Convert variables to numeric
    if(any(class(x) %in% c("character"))) x <- as.numeric(as_factor(x))
    if(any(class(y) %in% c("character"))) y <- as.numeric(as_factor(y))
    if(any(class(x) %in% c("Date", "factor", "ordered"))) x <- as.numeric(x)
    if(any(class(y) %in% c("Date", "factor", "ordered"))) y <- as.numeric(y)

    ## Save correlation estimates and p values
    if(length(unique(x)) < 5 & length(unique(y)) < 5) {
        ## Two categorical variables
        cvStat <- chisq.test(x, y, correct=FALSE)$statistic
        if(is.null(cvStat)) cvStat <- 0
        cv <- sqrt(cvStat/(length(x)*min(length(unique(x)),
            length(unique(y)))-1))
        corVal <- as.numeric(cv)
        corP <- chisq.test(x, y)$p.value
    } else {
        ## Otherwise, use Kendall's correlation coefficient and p-value
        corVal <- cor.test(x, y, method="kendall", exact=FALSE)$estimate
        corP <- cor.test(x, y, method="kendall", exact=FALSE)$p.value
    }

    ## Return as a data frame of estimate and p-value
    corTb <- tibble("corVal"=corVal, "corP"=corP)

    return(corTb)
}
