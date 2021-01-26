#' Sample Sheet Generation
#'
#' This function will generate visually intuitive plate layouts
#' for the wet lab, with the option to colour code different types
#' of samples (e.g. for studies investigating multiple tissues).
#'
#' @param omixerLayout Randomized sample list
#' @param group Colour-coding indicator
#'
#' @return PDF of sample layout in working directory
#'
#' @import ggplot2
#' @import magrittr
#' @import tibble
#' @import forcats
#' @import stringr
#' @importFrom dplyr select filter
#' @importFrom gridExtra marrangeGrob
#' @importFrom tidyselect everything all_of
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
#'                by="day"), 48))
#'                
#' randVars <- c("sex", "age", "smoke", "date")
#' 
#' omixerLayout <- omixerRand(sampleList, sampleId="sampleId", 
#' block="block", iterNum=10, wells=48, div="row", 
#' plateNum=1, randVars=randVars)
#' 
#' omixerSheet(omixerLayout)

omixerSheet <- function(omixerLayout=omixerLayout, group, group.text.size = 3.5, sample.text.size = 4) {

    ## Set labels
    sampleId <- NULL
    column <- NULL
    bottom <- NULL
    top <- NULL
    plate <- NULL

    omixerLayout <- omixerLayout %>% mutate(top=sampleId)
    if(!missing(group)) {
        omixerLayout <- omixerLayout %>% 
            select(bottom=all_of(group), everything())
    } else {
        omixerLayout$bottom <- NA
    }

    
    
    ## Create list of plate layouts
    ggPlateList <- lapply(seq_len(max(omixerLayout$plate)), function(x) {
        ggPlate <- omixerLayout %>% filter(plate == x) %>%
        ggplot(aes(x=column, y=row)) +
        geom_tile(aes(x=column, y=row, fill=factor(bottom)), colour="grey20", size=1.5,
            show.legend=FALSE) + coord_equal() +
        geom_text(aes(label=ifelse(is.na(bottom), "", as.character(bottom))),
            colour="grey30", size=group.text.size, nudge_y=0.2) +
        geom_text(aes(label=ifelse(is.na(top), "", as.character(top))),
            colour="grey30", fontface="bold", size=sample.text.size, nudge_y=-0.1) +
        scale_fill_brewer(palette="Set3") +
        scale_x_discrete(name="",
            limits=factor(c(min(as.numeric(omixerLayout$column)):
                         max(as.numeric(omixerLayout$column)))),
            position="top", expand=c(0,0)) +
        scale_y_discrete(name="", limits=(toupper(letters[
            c(max(as.numeric(omixerLayout$row)):
                min(as.numeric(omixerLayout$row)))])),
            expand=c(0,0)) +
        ggtitle(paste("Sample Overview for Plate", x)) +
        theme(plot.title=element_text(size=22, face="bold", hjust=0.5,
            colour="grey30"), axis.text=element_text(size=18, face="bold"),
            axis.ticks=element_blank(), panel.grid=element_blank(),
            panel.background=element_rect("grey80"),
            panel.border=element_rect("grey20", fill=NA, size=3))
        return(ggPlate)
    })

    ggList <- lapply(ggPlateList, ggplotGrob)

    ## Save layouts as a pdf
    ggsave("./omixer_sample_sheets.pdf", paper="a4", height=11, width=8,
        marrangeGrob(grobs=ggList, nrow=2, ncol=1, top=NULL))
}
