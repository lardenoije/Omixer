#' Multivariate Randomization
#'
#' As the main function of the Omixer package, this function
#' outputs a randomized sample list that minimizes correlations
#' between biological factors and technical covariates.
#'
#' @param df Sample list
#' @param sampleId String specifying sample ID variable
#' @param block Paired sample identifier
#' @param iterNum Number of layouts to generate
#' @param wells Number of wells on a plate
#' @param div Plate subdivisions
#' @param positional Logical indicator of positional batch effects
#' @param plateNum Number of plates
#' @param layout Custom plate layout as data frame
#' @param mask Wells to be left empty
#' @param techVars Technical covariates
#' @param randVars Randomization variables
#'
#' @return Selected randomized sample list as a data frame
#'
#' @import dplyr
#' @import ggplot2
#' @import magrittr
#' @importFrom tibble tibble
#' @importFrom readr write_delim write_csv write_csv2
#' @importFrom tidyselect everything all_of
#' @importFrom grid unit
#' @export

omixerRand <- function(df, sampleId="sampleId", block="block", iterNum=1000,
    wells, div="none", positional=FALSE, plateNum=1, layout, mask=0, techVars,
    randVars) {

  ## Set up plate layout
    if (!missing(layout) & !missing(techVars)) {
        layout <- layout; techVars <- techVars
    } else if (!missing(layout) & missing(techVars)) {
        stop("For custom layouts you must supply technical covariates.")
    } else if (missing(layout) & !missing(wells)) {
        rowNum <- case_when(wells == 96 ~ 8, wells == 48 ~ 6, wells == 24 ~ 4)
        colNum <- case_when(wells == 96 ~ 12, wells == 48 ~ 8, wells == 24 ~ 6)

    if (!wells %in% c(96, 48, 24)) {
        stop("Automated layouts only support 96, 48, or 24 well plates.")
    }
    if(!div %in% c("none", "col", "row", "col-block", "row-block")) {
        stop("Please specify a valid div (see documentation for options).")
    }
    if (positional == TRUE && div == "none") {
        stop("Positional batches not allowed without plate subdivisions.")
    }
    well <- NULL
    plate <- NULL
    layout <- tibble(plate=rep(1:plateNum, each=wells),
        well=rep(1:wells, plateNum),
        row=factor(rep(rep(1:rowNum, each=colNum), plateNum),
            labels=toupper(letters[1:rowNum])),
        column=rep(1:colNum, rowNum*plateNum), mask=mask,
        chip=case_when(div == "col" ~ column, div == "row" ~ as.integer(row),
            div == "col-block" ~ as.integer(ceiling(column / 2)),
            div == "row-block" ~ as.integer(ceiling(as.numeric(row) / 2))),
        chipPos = case_when(div == "col" ~ as.numeric(row),
            div == "row" ~ as.numeric(column),
            div == "col-block" ~ ifelse(column %% 2 == 0,
                as.numeric(row) + rowNum, row),
            div == "row-block" ~ ifelse(as.numeric(row) %% 2 == 0,
                column + colNum, column)))
    } else {
        stop("You must either specify a custom layout or a number of wells.")
    }

    ## Set up technical covariates
    if (plateNum == 1){
        if (div == "none"){
            stop("With one plate and no subdivisions, there are no batches.")
        } else {
            if (positional == TRUE) {
                techVars <- c("chip", "chipPos")
            } else {
                techVars <- "chip"
            }
        }
    } else if (plateNum > 1) {
        if (div == "none") {
            techVars <- "plate"
        } else {
            if (positional == TRUE) {
                techVars <- c("plate", "chip", "chipPos")
            } else {
                techVars <- c("plate", "chip")
            }
        }
    } else {
        stop("Plate number must be a positive integer.")
    }

    ## Define sample ID, blocks, and permutation variables
    permVar <- NULL
    if(block %in% colnames(df)) {
        df <- df %>% select(sampleId=all_of(sampleId), block=all_of(block),
            permVar = all_of(block), everything())
    } else if(sampleId %in% colnames(df)) {
        df <- df %>% select(sampleId=all_of(sampleId),
            permVar = all_of(sampleId), everything())
    } else {
        stop("Sample ID not found in provided sample list.")
    }

    ## Create the specified number of permutations of block or ID
    permSet <- lapply(1:iterNum, function(x) {
        set.seed(x*999)
        sample(unique(df$permVar))
    })

    ## Create list of randomized sample layouts using permSet
    dfRandList <- lapply(1:iterNum, function(x) {
        set.seed(x*999)
        dfShuffle <- df %>% group_by(permVar) %>% slice(sample(1:n())) %>%
            ungroup()
        dfRandList <- lapply(1:length(permSet[[x]]), function(y){
          dfRand <- dfShuffle %>% filter(permVar == permSet[[x]][y])
          return(dfRand)
        })
        dfRand <- bind_rows(dfRandList)
        dfRand <- dfRand %>% mutate(seed=x*999)
        return(dfRand)
    })

    ## Filter masked wells from plate layout
    layoutMasked <- layout %>% filter(mask == 0)
    if(nrow(layoutMasked) != nrow(dfRandList[[1]])) {
        stop("Number of unmasked wells must equal number of samples.")
    }

    ## Combine randomized sample lists with plate layout
    sampleLayoutList <- lapply(1:length(dfRandList), function(x) {
        sampleLayout <- cbind(dfRandList[[x]], layoutMasked)
        return(sampleLayout)
    })

    ## If randomization variables are not specified, then use all except IDs
    if(missing(randVars)){
        randVars <- colnames(dfRandList[[1]])[!colnames(dfRandList[[1]]) %in%
            c("sampleId", "blockVar")]
    }

    ## Save correlation estimates and p-values
    corVal <- NULL
    corP <- NULL
    corTbList <- lapply(1:length(sampleLayoutList), function(x){
        sampleLayout <- sampleLayoutList[[x]]
        corTbList <- lapply(randVars, function(y){
            corTbList <- lapply(techVars, function(z){
                cor <- omixerCorr(sampleLayout[, y], sampleLayout[, z])
                corTb <- tibble(randVars=y, techVars=z, corVal=cor$corVal,
                    corP=cor$corP)
                return(corTb)
            })
            corTb <-bind_rows(corTbList)
            return(corTb)
        })
        corTb <- bind_rows(corTbList)
        return(corTb)
    })

    ## Create correlation table
    pTest <- NULL
    absSum <- NULL
    corSumList <- lapply(1:length(corTbList), function(x){
        corTb <- corTbList[[x]]
        corSum <- tibble(seed=x*999, absSum=sum(abs(corTb$corVal)),
            pTest=any(corTb$corP < 0.05))
    })
    corSum <- bind_rows(corSumList)

    #Find the optimal layout
    chosenSeed <- (corSum %>% filter(pTest == FALSE) %>%
        filter(absSum == min(absSum)))$seed

    ## Check number of optimized layouts
    if(length(chosenSeed)==0) {
        warning("All randomized layouts contained unwanted correlations.")
        warning("Returning best possible layout.")
        nonoptSeed <- (corSum %>% filter(absSum == min(absSum)))$seed
    }
    if(length(chosenSeed) > 1) {
        chosenSeed <- chosenSeed[1]
        message("Several layouts were equally optmized.")
    }

    ## Save correlations for chosen layout
    if(length(chosenSeed)>0){
        corSelect <- corTbList[[chosenSeed/999]]
        omixerLayout <- sampleLayoutList[[chosenSeed/999]]
    } else {
        corSelect <- corTbList[[nonoptSeed/999]]
        omixerLayout <- sampleLayoutList[[nonoptSeed/999]]
    }

    ## Rejoin layout with masked wells
    omixerLayout <- full_join(omixerLayout, layout,
        by=c("well", "plate", "row", "column", "mask", "chip", "chipPos"))
    omixerLayout <- omixerLayout %>% arrange(plate, well)

    ## Print information
    print(paste("Layout created using a seed of:", omixerLayout$seed[1]))

    ## Visualize correlations
    print(ggplot(corSelect, aes(x=randVars, y=techVars)) +
    geom_tile(aes(fill=corVal), size=3, colour="white", show.legend=FALSE) +
    geom_text(aes(label=round(corVal, 3)),
        colour=ifelse(corSelect$corVal < mean(corSelect$corVal), "white",
        "grey30"), fontface="bold", nudge_y=0.2, size=8) +
    geom_text(aes(label=paste("p =", round(corP, 3))),
        colour=ifelse(corSelect$corVal < mean(corSelect$corVal), "white",
        "grey30"), nudge_y=-0.2, size=6) +
    scale_fill_distiller(palette="YlGnBu") +
    scale_x_discrete(position="top", name="Randomization variables \n",
        label=function(x) abbreviate(x, minlength=6), expand=c(0,0)) +
    scale_y_discrete(name="Technical \n covariates",
        label=function(x) abbreviate(x, minlength=6), expand=c(0,0)) +
    ggtitle("Correlations present in the chosen layout") + coord_equal() +
    theme(plot.title=element_text(hjust=0.5, size=24),
        axis.title=element_text(face="bold", size=18),
        axis.ticks=element_blank(),
        axis.text.x=element_text(size=16),
        axis.text.y=element_text(angle=90, size=16, vjust=1)))

  return(omixerLayout)
}

