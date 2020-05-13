#' Sample Sheet Generation
#'
#' This function takes a sample layout, either specified or
#' from the omix_rand() function and generates a PDF of the
#' layout to aid with manual pipetting.
#'
#' @param layout List of randomized samples
#' @param group (Optional) Colour-coding indicator
#'
#' @return PDF of sample layout in working directory
#'
#' @importFrom dplyr mutate filter
#' @importFrom ggplot2 ggplot geom_tile geom_text scale_fill_brewer scale_x_discrete ggsave scale_y_discrete ggtitle theme element_blank element_text element_rect
#' @importFrom gridExtra marrangeGrob
#' @import magrittr
#' @export

omix_sheet <- function(layout = final_layout, group) {

  # Initialize variables
  sample_id <- NULL
  plate <- NULL
  aes <- NULL
  columns <- NULL
  bottom <- NULL
  top <- NULL

  # Set labels
  final_layout <- final_layout %>%
    mutate(
      top = sample_id,
      bottom = group)

  # Initialize variables
  ggPlateList <- list()
  ggIndex <- 0

  # Loop over plates
  for(p in 1:7) {
    ggIndex <- ggIndex + 1


    if(!is.null(final_layout$bottom)) {
      # If a group is specified, use colour coding
      ggPlateList[[ggIndex]] <- final_layout %>%
        filter(plate == p) %>%
        ggplot(aes(x=columns, y=row)) +
        geom_tile(aes(x=columns, y=row, fill = bottom), colour = "grey30", size = 1.5, show.legend = FALSE) +
        geom_text(aes(label=bottom), size=3.5, nudge_y=0.2, colour="grey30") +
        geom_text(aes(label=top), size=4, nudge_y=-0.1, colour="grey30", fontface="bold") +
        scale_fill_brewer(palette = "Set3") +
        scale_x_discrete(name="", limits=c(1:12), expand=c(0.01,0.01), position="top") +
        scale_y_discrete(name="", limits=toupper(letters[8:1]),expand=c(0.08,0.08)) +
        ggtitle(paste("Sample Overview for Plate",p))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_rect("grey30"),
              plot.title = element_text(size = 22, face = "bold"),
              axis.text = element_text(size=18, face="bold"))
    } else {
      # If no group is specified
      ggPlateList[[ggIndex]] <- final_layout %>%
        filter(plate == p) %>%
        ggplot(aes(x=columns, y=row)) +
        geom_tile(aes(x=columns, y=row),  fill = "grey80", colour = "grey30", size = 1.5, show.legend = FALSE) +
        geom_text(aes(label=top), size=4, colour="grey30", fontface="bold") +
        scale_x_discrete(name="", limits=c(1:12), expand=c(0.01,0.01), position="top") +
        scale_y_discrete(name="", limits=toupper(letters[8:1]),expand=c(0.08,0.08)) +
        ggtitle(paste("Sample Overview for Plate",p))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_rect("grey30"),
              plot.title = element_text(size = 22, face = "bold"),
              axis.text = element_text(size=18, face="bold"))
    }
  }

  # Save layouts as a pdf
  ggsave("C:/Users/ljsinke/Documents/5 DIMENSION/Bioconductor Package/sample_sheets_final.pdf",
         paper = "a4", height = 11, width = 8,
         marrangeGrob(grobs = ggPlateList, nrow=2, ncol=1, top=NULL))

  return(ggPlateList)

}
