#' Sample List Regeneration
#'
#' Regenerate an Omixer-produced randomized sample list quickly,
#' by providing the seed returned from omixerRand().
#'
#' @param df Sample list
#' @param seed Seed
#' @param sampleId String specifying sample ID variable
#' @param block Paired sample identifier
#' @param wells Number of wells on a plate
#' @param div Plate subdivisions
#' @param positional Logical indicator of positional batch effects
#' @param plateNum Number of plates
#' @param layout Custom plate layout as data frame
#' @param mask Wells to be left empty
#' @param techVars Technical covariates
#' @param randVars Randomization variables
#'
#' @return Chosen layout as a data frame
#'
#' @import dplyr
#' @import ggplot2
#' @import magrittr
#' @importFrom tibble tibble
#' @importFrom readr write_delim write_csv write_csv2
#' @importFrom tidyselect everything all_of
#' @importFrom grid unit
#' @export

omixerSpecific <- function(df, seed=seed, sampleId="sampleId", block="block",
    wells, div="none", positional=FALSE, plateNum=1, layout, mask=0, techVars,
    randVars) {

    ## Create Plate Layout
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
    plate <- NULL
    well <- NULL
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
            permVar=all_of(block), everything())
    } else if(sampleId %in% colnames(df)) {
        df <- df %>% select(sampleId=all_of(sampleId),
            permVar=all_of(sampleId), everything())
    } else {
        stop("Sample ID not found in provided sample list.")
    }

    ## Use seed to generate permSet of the specified layout
    set.seed(seed)
    permSet <- sample(unique(df$permVar))

    ## Create randomized sample layout using permSet
    set.seed(seed)
    dfShuffle <- df %>% group_by(permVar) %>% slice(sample(1:n())) %>%
        ungroup()

    dfRandList <- lapply(1:length(permSet), function(x){
      dfRand <- dfShuffle %>% filter(permVar == permSet[[x]])
      return(dfRand)
    })
    dfRand <- bind_rows(dfRandList)

    ## Set up mask
    layout$mask <- mask
    layoutMasked <- layout %>% filter(mask == 0)

    ## Check if number of unmasked wells is equal to number of samples
    if(nrow(layoutMasked) != nrow(dfRand)) {
        stop("Number of unmasked wells must equal number of samples")
    }

    ## Combine randomized sample lists with plate layout
    omixerLayout <- cbind(dfRand, layoutMasked)

    ## If randomization variables are not specified, then use all except IDs
    if(missing(randVars)){
        randVars <- colnames(dfRand)[!colnames(dfRand) %in%
            c("sampleId", "block")]
    }

    ## Save correlation estimates and p-values
    corVal <- NULL
    corP <- NULL
    corTb <- NULL
    corTbList <- lapply(randVars, function(y){
        corTbList <- lapply(techVars, function(z){
            cor <- omixerCorr(omixerLayout[, y], omixerLayout[, z])
            corTb <- tibble(randVars=y, techVars=z, corVal=cor$corVal,
                corP=cor$corP)
            return(corTb)
        })
        corTb <-bind_rows(corTbList)
        return(corTb)
    })

    ## Create correlation table
    corSum <- tibble(seed=seed, absSum=sum(abs(corTb$corVal)),
        pTest=any(corTb$corP < 0.05))

    ## Rejoin layout with masked wells
    omixerLayout <- full_join(omixerLayout, layout,
        by=c("well", "plate", "row", "column", "mask", "chip", "chipPos"))
    omixerLayout <- omixerLayout %>% arrange(plate, well)

    ## Print information
    print(paste("This sample layout was regenerated using a seed of", seed))

    ## Visualize correlations
    print(ggplot(corTb, aes(x=randVars, y=techVars)) +
        geom_tile(aes(fill=corVal), size=3, colour="white",
            show.legend=FALSE) +
        geom_text(aes(label=round(corVal, 3)),
            colour=ifelse(corTb$corVal < mean(corTb$corVal), "white",
            "grey30"), fontface="bold", nudge_y=0.2, size=8) +
        geom_text(aes(label=paste("p =",round(corP, 3))),
            colour=ifelse(corTb$corVal < mean(corTb$corVal), "white",
            "grey30"), nudge_y=-0.2, size=6) +
        scale_fill_distiller(palette="YlGnBu") +
        scale_x_discrete(position="top", name="Randomization variables \n",
            label=function(x) abbreviate(x, minlength=6), expand=c(0,0)) +
        scale_y_discrete(name="Technical \n covariates",
            label=function(x) abbreviate(x, minlength=6), expand=c(0,0)) +
        ggtitle("Correlations present in the chosen layout") + coord_equal() +
        theme(plot.title=element_text(hjust=0.5, size=24),
            axis.title=element_text(face="bold", size=18),
            axis.ticks=element_blank(), axis.text.x=element_text(size=16),
            axis.text.y=element_text(angle=90, size=16, vjust=1)))

    return(omixerLayout)
}
