#' Sample Sheet Generation
#'
#' This function will generate visually intuitive plate layouts
#' for the wet lab, with the option to colour code different types
#' of samples (e.g. for studies investigating multiple tissues).
#'
#' @param omixr_layout Randomized sample list
#' @param group Colour-coding indicator
#'
#' @return PDF of sample layout in working directory
#'
#' @importFrom dplyr select filter
#' @importFrom ggplot2 ggplot geom_tile geom_text scale_fill_brewer scale_x_discrete ggsave scale_y_discrete ggtitle theme element_blank element_text element_rect ggplotGrob
#' @importFrom gridExtra marrangeGrob
#' @importFrom tidyselect everything all_of
#' @import magrittr
#' @export

omixr_sheet <- function(omixr_layout = omixr_layout, group = NULL) {

  # Initialize variables
  sample_id <- NULL
  plate <- NULL
  columns <- NULL
  bottom <- NULL
  top <- NULL

  # Set labels

  omixr_layout <- omixr_layout %>%
      select(
        top = sample_id,
        bottom = ifelse(is.null(group), NA, all_of(group)),
        everything()
      )

  # Create list of plate layouts
  ggPlateList <- lapply(1:max(omixr_layout$plate), function(x) {
      ggPlate <- omixr_layout %>%
        filter(plate == x) %>%
        ggplot(aes(x=columns, y=row)) +
        geom_tile(aes(x = columns, y = row, fill = bottom), colour = "grey20", size = 1.5, show.legend = FALSE) +
        coord_equal() +
        geom_text(aes(label=ifelse(is.na(bottom),"",as.character(bottom))), colour = "grey30", size = 3.5, nudge_y = 0.2) +
        geom_text(aes(label=ifelse(is.na(top),"",as.character(top))), colour = "grey30", fontface = "bold", size = 4, nudge_y = -0.1) +
        scale_fill_brewer(palette = "Set3") +
        scale_x_discrete(name = "", limits = c(min(omixr_layout$columns):max(omixr_layout$columns)), position = "top", expand = c(0,0)) +
        scale_y_discrete(name = "", limits = toupper(letters[max(as.numeric(omixr_layout$row)):min(as.numeric(omixr_layout$row))]), expand = c(0,0)) +
        ggtitle(paste("Sample Overview for Plate", x)) +
        theme(plot.title = element_text(size = 22, face = "bold", hjust = 0.5, colour = "grey30"),
              axis.text = element_text(size = 18, face = "bold"),
              axis.ticks = element_blank(),
              panel.grid = element_blank(),
              panel.background = element_rect("grey50"),
              panel.border = element_rect("grey20", fill = NA, size = 3))
    return(ggPlate)
  })

  glist <- lapply(ggPlateList, ggplotGrob)

  # Save layouts as a pdf
  ggsave("./omixr_sample_sheets.pdf",
         paper = "a4", height = 11, width = 8,
         marrangeGrob(grobs = glist, nrow=2, ncol=1, top=NULL))
}
