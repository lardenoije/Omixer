#' Sample List Regeneration
#'
#' Regenerate an omixr-produced randomized sample list quickly,
#' by providing the seed returned from omixr_rand().
#'
#' @param df Sample list
#' @param seed Seed
#' @param sample_id Name of sample ID variable
#' @param block Paired sample identifier
#' @param plate_layout Pre-specified or custom plate layout
#' @param rand_vars Randomization variables (default: all except sample ID and block)
#' @param tech_vars Technical covariates
#' @param mask Wells to be left empty
#'
#' @return Chosen layout as a data frame
#'
#' @importFrom dplyr select group_by slice ungroup bind_rows filter full_join arrange n
#' @importFrom readr write_delim write_csv write_csv2
#' @importFrom tidyselect everything all_of
#' @importFrom ggplot2 aes geom_point scale_fill_distiller scale_size_continuous coord_equal
#' @importFrom grid unit
#' @import magrittr
#' @export

omixr_specific <- function(df, seed = seed, sample_id = sample_id, block = block, plate_layout = "Illumina96", rand_vars, tech_vars, mask = 0) {

  # Initialize variables
  perm_var <- NULL
  plate <- NULL
  well <- NULL
  rand_var <- NULL
  tech_var <- NULL
  corr_val <- NULL
  corr_p <- NULL

  # Print information
  print("Initializing randomization package.")

  # Set up predefined plate layout
  if(plate_layout == "Illumina96") {
    plate_layout <- data.frame(
      well = rep(1:96, 7),
      plate = rep(1:7, each=96),
      row = factor(rep(1:8, 84), labels = toupper(letters[1:8])),
      columns = rep(rep(1:12, each=8), 7))
  }

  # Define sample ID, blocks, and permutation variables
  # Check that sample ID exists in the provided data frame
  if(block %in% colnames(df)) {
    df <- df %>%
      select(sample_id = all_of(sample_id),
             block = all_of(block),
             perm_var = all_of(block),
             everything())
  } else if(sample_id %in% colnames(df)) {
    df <- df %>%
      select(sample_id = all_of(sample_id),
             perm_var = all_of(sample_id),
             everything())
  } else {
    stop("Sample ID not found in provided sample list. Please include a Sample ID.")
  }

  # Use seed to generate perm_set of the specified layout
  set.seed(seed)

  perm_set <- sample(unique(df$perm_var))

  # Print information
  print("Randomized sample list is being regenerated.")

  # Create randomized sample layout using perm_set
  set.seed(seed)

  df_shuffle <- df %>%
    group_by(perm_var) %>%
    slice(sample(1:n())) %>%
    ungroup()

  df_rand <- data.frame()

  for(i in 1:length(perm_set)) {
    df_rand <- bind_rows(df_rand, df_shuffle %>% filter(perm_var == perm_set[i]))
    }

  # Set up mask
  plate_layout$mask <- mask

  # Filter masked wells from plate layout
  plate_layout_masked <- plate_layout %>% filter(mask == 0)

  # Check if number of unmasked wells is equal to number of samples
  if(nrow(plate_layout_masked) != nrow(df_rand)) {
    stop("Number of unmasked wells must equal number of samples")
  }

  # Print information
  print("Combining sample list with plate layout.")

  # Combine randomized sample lists with plate layout
  omixr_layout <- cbind(df_rand, plate_layout_masked)

  # If randomization variables are not specified, then use all columns except ID and block
  if(!exists("rand_vars")){
    rand_vars = colnames(df_rand)[!colnames(df_rand) %in% c("sample_id", "block_var")]
  }

  # Print information
  print("Calculating correlations.")

  # Save correlation estimates and p-values
  corr_df <- data.frame()

  for(i in rand_vars){
    for(j in tech_vars){
      corr <- omixr_corr(omixr_layout[, i], omixr_layout[, j])
      corr_df <- bind_rows(corr_df, data.frame(rand_var = i, tech_var = j, corr_val = corr$corr_val, corr_p = corr$corr_p))
      }
    }

  # Create data frame with seed, absolute sum of correlations, and whether p-value is under 0.05
  corr_sum <- data.frame(seed = seed, abs_sum = sum(abs(corr_df$corr_val)), p_test = any(corr_df$corr_p < 0.05))

  # Rejoin layout with masked wells
  omixr_layout <- full_join(omixr_layout, plate_layout, by = c("well", "plate", "row", "columns", "mask"))
  omixr_layout <- omixr_layout %>% arrange(plate, well)

  # Print information
  print(paste("This sample layout was regenerated using a seed of", seed))

  # Visualize correlations
  print(ggplot(corr_df, aes(x = rand_var, y = tech_var)) +
          geom_tile(aes(fill = corr_val),
                    size = 3, colour = "white", show.legend = FALSE) +
          geom_text(aes(label = round(corr_val,3)),
                    colour = ifelse(corr_df$corr_val < mean(corr_df$corr_val), "white", "grey30"), fontface = "bold", nudge_y = 0.2, size = 5) +
          geom_text(aes(label = paste("p =",round(corr_p, 3))),
                    nudge_y = -0.2, colour = ifelse(corr_df$corr_val < mean(corr_df$corr_val), "white", "grey30")) +
          scale_fill_distiller(palette = "YlGnBu") +
          scale_x_discrete(position = "top",
                           name = "Randomization variables \n",
                           label=function(x) abbreviate(x, minlength=6),
                           expand = c(0,0)) +
          scale_y_discrete(name = "Technical \n covariates",
                           label=function(x) abbreviate(x, minlength=6),
                           expand = c(0,0)) +
          ggtitle("Correlations present in the chosen layout") +
          coord_equal() +
          theme(plot.title = element_text(hjust = 0.5, size = 18),
                axis.title = element_text(face = "bold", size = 14),
                axis.ticks = element_blank(),
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(angle = 90, size = 12, vjust=1)))


  return(omixr_layout)
}



