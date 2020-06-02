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
#' @importFrom dplyr select filter
#' @importFrom gridExtra marrangeGrob
#' @importFrom tidyselect everything all_of
#' @import magrittr
#' @export

omixerSheet <- function(omixerLayout=omixerLayout, group=NULL) {

    ## Set labels
    sampleId <- NULL
    column <- NULL
    bottom <- NULL
    top <- NULL
    plate <- NULL

    omixerLayout <- omixerLayout %>% select(top = sampleId,
        bottom = ifelse(is.null(group), NA, all_of(group)), everything())

    ## Create list of plate layouts
    ggPlateList <- lapply(1:max(omixerLayout$plate), function(x) {
        ggPlate <- omixerLayout %>% filter(plate == x) %>%
        ggplot(aes(x=column, y=row)) +
        geom_tile(aes(x=column, y=row, fill=bottom), colour="grey20", size=1.5,
            show.legend=FALSE) + coord_equal() +
        geom_text(aes(label=ifelse(is.na(bottom), "", as.character(bottom))),
            colour="grey30", size=3.5, nudge_y=0.2) +
        geom_text(aes(label=ifelse(is.na(top), "", as.character(top))),
            colour="grey30", fontface="bold", size=4, nudge_y=-0.1) +
        scale_fill_brewer(palette="Set3") +
        scale_x_discrete(name="",
            limits=c(min(omixerLayout$column):max(omixerLayout$column)),
            position="top", expand=c(0,0)) +
        scale_y_discrete(name="", limits=toupper(letters[max(as.numeric(
            omixerLayout$row)):min(as.numeric(omixerLayout$row))]),
            expand=c(0,0)) +
        ggtitle(paste("Sample Overview for Plate", x)) +
        theme(plot.title=element_text(size=22, face="bold", hjust=0.5,
            colour="grey30"), axis.text=element_text(size=18, face="bold"),
            axis.ticks=element_blank(), panel.grid=element_blank(),
            panel.background=element_rect("grey50"),
            panel.border=element_rect("grey20", fill=NA, size=3))
        return(ggPlate)
    })

    ggList <- lapply(ggPlateList, ggplotGrob)

    ## Save layouts as a pdf
    ggsave("./omixer_sample_sheets.pdf", paper="a4", height=11, width=8,
        marrangeGrob(grobs=ggList, nrow=2, ncol=1, top=NULL))
}
