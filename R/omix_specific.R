#' Sample List Generation
#'
#' This function takes a sample list and seed, and returns the
#' layout that would have been produced by omix_rand() with this seed.
#'
#' This is for regeneration of randomized sample lists if, for example,
#' adding new columns that would not affect optimal randomization.
#'
#' @param df Sample list
#' @param seed Seed returned by omix_rand()
#' @param sample_id Name of sample ID variable
#' @param block (Optional) Paired sample identifier
#' @param plate_layout Pre-specified or custom plate layout
#' @param rand_vars (Optional) Randomization variables (default: all except sample ID and blocking variable)
#' @param tech_vars Technical covariates
#' @param mask Wells to be left empty
#' @param save Method of saving output
#'
#' @return Data frame containing the chosen layout
#'
#' @importFrom dplyr select group_by slice ungroup bind_rows filter full_join arrange
#' @importFrom readr write_delim write_csv write_csv2
#' @importFrom tidyselect everything all_of
#' @importFrom ggplot2 aes geom_point scale_fill_distiller scale_size_continuous
#' @importFrom grid unit
#' @import magrittr
#' @export

omix_specific <- function(df,
                          seed = seed,
                          sample_id = sample_id,
                          block = block,
                          plate_layout = "Illumina96",
                          rand_vars,
                          tech_vars,
                          mask = 0,
                          save = "rdata"
) {

  # Initialize variables
  perm_var <- NULL
  p_test <- NULL
  corr_sum <- NULL
  plate <- NULL
  well <- NULL
  n <- NULL
  abs_sum <- NULL
  rand_var <- NULL
  tech_var <- NULL
  corr_val <- NULL

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
    print("Order of paired samples will be randomized.")
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
  sample_layout <- cbind(df_rand, plate_layout_masked)

  # If randomization variables are not specified, then use all columns except ID and block
  if(!exists("rand_vars")){
    rand_vars = colnames(df_rand_list[[1]])[!colnames(df_rand_list[[1]]) %in% c("sample_id", "block_var")]
  }

  # Print information
  print("Calculating correlations.")

  # Save correlation estimates and p-values
  corr_df <- data.frame()

  for(i in rand_vars){
    for(j in tech_vars){
      corr <- omix_corr(sample_layout[, i], sample_layout[, j])
      corr_df <- bind_rows(corr_df, data.frame(rand_var = i, tech_var = j, corr_val = corr$corr_val, corr_p = corr$corr_p))
      }
    }

  # Create data frame with seed, absolute sum of correlations, and whether p-value is under 0.05
  corr_sum <- data.frame(seed = seed, abs_sum = sum(abs(corr_df$corr_val)), p_test = any(corr_df$corr_p < 0.05))

  # Select the layout with minimum correlations and no significant correlations
  final_layout <- sample_layout

  # Rejoin layout with masked wells
  final_layout <- full_join(final_layout, plate_layout, by = c("well", "plate", "row", "columns", "mask"))
  final_layout <- final_layout %>% arrange(plate, well)

  # Print information
  print(paste("Selected layout created using a seed of:", seed))

  # Visualize correlations
  print(ggplot(corr_df, aes(x = rand_var, y = tech_var)) +
          geom_tile(colour = "white", size = 3, fill = "grey90") +
          geom_point(aes(size = corr_val, fill = corr_val), stroke = 2, color = "grey20", shape = 21, show.legend = FALSE) +
          geom_text(aes(label = round(corr_val,3)), colour = "grey20", fontface = "bold") +
          scale_fill_distiller(palette = "Spectral") +
          scale_size_continuous(range = c(15,25)) +
          scale_x_discrete(position = "top", name = "Randomization variables", label=function(x) abbreviate(x, minlength=8), expand = c(0,0)) +
          scale_y_discrete(name = "Technical covariates", label=function(x) abbreviate(x, minlength=8), expand = c(0,0)) +
          ggtitle("Correlations present in the chosen layout") +
          theme(plot.title = element_text(hjust = 0.5, size = 16),
                panel.background = element_blank(),
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(angle = 90, size = 12),
                plot.margin = unit(c(1,1,1,1), "cm")))

  return(final_layout)

  # Save layout and selected
  if(save == "csv2") {
    write_csv2(final_layout, "./final_layout.csv")
  } else if(save == "csv") {
    write_csv(final_layout, "./final_layout.csv")
  } else if(save == "tab") {
    write_delim(final_layout, "./final_layout.txt", delim="\t")
  } else {
    save(final_layout, file = "./final_layout.Rdata")
  }
}



