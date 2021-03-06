% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/omixerRand.R
\name{omixerRand}
\alias{omixerRand}
\title{Multivariate Randomization}
\usage{
omixerRand(
  df,
  sampleId = "sampleId",
  block = "block",
  iterNum = 1000,
  wells,
  div = "none",
  positional = FALSE,
  plateNum = 1,
  layout,
  mask = 0,
  techVars,
  randVars
)
}
\arguments{
\item{df}{Sample list}

\item{sampleId}{String specifying sample ID variable}

\item{block}{Paired sample identifier}

\item{iterNum}{Number of layouts to generate}

\item{wells}{Number of wells on a plate}

\item{div}{Plate subdivisions}

\item{positional}{Logical indicator of positional batch effects}

\item{plateNum}{Number of plates}

\item{layout}{Custom plate layout as data frame}

\item{mask}{Wells to be left empty}

\item{techVars}{Technical covariates}

\item{randVars}{Randomization variables}
}
\value{
Selected randomized sample list as a data frame

Randomization environment of optimal list generation
}
\description{
As the main function of the Omixer package, this function
outputs a randomized sample list that minimizes correlations
between biological factors and technical covariates.
}
\examples{
library(tibble)
library(forcats)
library(stringr)

sampleList <- tibble(sampleId=str_pad(1:48, 4, pad="0"),
sex=as_factor(sample(c("m", "f"), 48, replace=TRUE)), 
age=round(rnorm(48, mean=30, sd=8), 0), 
smoke=as_factor(sample(c("yes", "ex", "never"), 48, replace=TRUE)),
date=sample(seq(as.Date('2008/01/01'), as.Date('2016/01/01'), 
               by="day"), 48))
               
randVars <- c("sex", "age", "smoke", "date")

omixerLayout <- omixerRand(sampleList, sampleId="sampleId", 
block="block", iterNum=10, wells=48, div="row", 
plateNum=1, randVars=randVars)
}
