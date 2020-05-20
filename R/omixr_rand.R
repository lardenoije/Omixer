#' Multivariate Randomization
#'
#' As the main function of the omixr package, this function
#' outputs a randomized sample list that minimizes correlations
#' between biological factors and technical covariates.
#'
#' @param df Sample list
#' @param sample_id String specifying sample ID variable
#' @param block Paired sample identifier
#' @param iter_num Number of layouts to generate
#' @param plate_layout Pre-specified or custom plate layout
#' @param rand_vars Randomization variables
#' @param tech_vars Technical covariates
#' @param mask Wells to be left empty
#'
#' @return Selected randomized sample list as a data frame
#'
#' @importFrom dplyr select group_by slice ungroup bind_rows filter full_join arrange n
#' @importFrom readr write_delim write_csv write_csv2
#' @importFrom tidyselect everything all_of
#' @importFrom ggplot2 aes geom_point scale_fill_distiller scale_size_continuous coord_equal
#' @importFrom grid unit
#' @import magrittr
#' @export

omixr_rand <- function(df, sample_id = sample_id, block = block, iter_num = 10000, plate_layout = "Illumina96", rand_vars, tech_vars, mask = 0) {

  # Initialize variables
  perm_var <- NULL
  p_test <- NULL
  abs_sum <- NULL
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
    print("Order of paired samples will be randomized.")
  } else if(sample_id %in% colnames(df)) {
    df <- df %>%
      select(sample_id = all_of(sample_id),
             perm_var = all_of(sample_id),
             everything())
  } else {
    stop("Sample ID not found in provided sample list. Please include a Sample ID.")
  }

  # Create the specified number of permutations of block or ID
  perm_set <- lapply(1:iter_num, function(x) {
    set.seed(x)
    sample(unique(df$perm_var))
  })

  # Print information
  print("Randomization in progress. Please be patient.")

  # Create list of randomized sample layouts using perm_set
  df_rand_list <- lapply(1:iter_num, function(x) {
    set.seed(x)
    df_shuffle <- df %>%
      group_by(perm_var) %>%
      slice(sample(1:n())) %>%
      ungroup()
    df_rand <- data.frame()
    for(i in 1:length(perm_set[[x]])) {
      df_rand <- bind_rows(df_rand, df_shuffle %>% filter(perm_var == perm_set[[x]][i]))
    }
    df_rand$seed <- x
    return(df_rand)
  })

  # Print information
  print("Randomization complete.")

  # Set up mask
  plate_layout$mask <- mask

  # Filter masked wells from plate layout
  plate_layout_masked <- plate_layout %>% filter(mask == 0)

  # Check if number of unmasked wells is equal to number of samples
  if(nrow(plate_layout_masked) != nrow(df_rand_list[[1]])) {
    stop("Number of unmasked wells must equal number of samples")
  }

  # Print information
  print("Combining sample list with plate layout.")

  # Combine randomized sample lists with plate layout
  sample_layout_list <- lapply(1:length(df_rand_list), function(x) {
    sample_layout <- cbind(df_rand_list[[x]], plate_layout_masked)
    return(sample_layout)
  })

  # If randomization variables are not specified, then use all columns except ID and block
  if(!exists("rand_vars")){
    rand_vars = colnames(df_rand_list[[1]])[!colnames(df_rand_list[[1]]) %in% c("sample_id", "block_var")]
  }

  # Print information
  print("Calculating correlations.")

  # Save correlation estimates and p-values
  corr_df_list <- lapply(1:length(sample_layout_list), function(x){
    sample_layout <- sample_layout_list[[x]]
    corr_df <- data.frame()
    for(i in rand_vars){
      for(j in tech_vars){
        corr <- omixr_corr(sample_layout[, i], sample_layout[, j])
        corr_df <- bind_rows(corr_df, data.frame(rand_var = i, tech_var = j, corr_val = corr$corr_val, corr_p = corr$corr_p))
      }
    }
    return(corr_df)
  })

  # Create data frame with seed, absolute sum of correlations, and whether p-value is under 0.05
  corr_sum <- data.frame()

  for(i in 1:length(corr_df_list)){
    corr_df <- corr_df_list[[i]]
    corr_sum <- bind_rows(corr_sum, data.frame(seed = i, abs_sum = sum(abs(corr_df$corr_val)), p_test = any(corr_df$corr_p < 0.05)))
  }

  # Find the layout with the minimum sum of correlations and no significant correlation test
  chosen_seed <- (corr_sum %>% filter(p_test == FALSE) %>% filter(abs_sum == min(abs_sum)))$seed

  # Check number of optimized layouts
  if(length(chosen_seed)==0){
    print("All randomized layouts provided sufficient evidence for correlation between technical covariates and randomization variables.")
    print("Please increase iteration number or reduce number of randomization and/or technical covariates.")
  } else if(length(chosen_seed) > 1) {
    chosen_seed <- chosen_seed[1]
    print("Several layouts were equally optmized, so only the first will be selected.")
  } else {
  }

  # Save correlations for chosen layout
  corr_select <- corr_df_list[[chosen_seed]]

  # Select the layout with minimum correlations and no significant correlations
  omixr_layout <- sample_layout_list[[chosen_seed]]

  # Rejoin layout with masked wells
  omixr_layout <- full_join(omixr_layout, plate_layout, by = c("well", "plate", "row", "columns", "mask"))
  omixr_layout <- omixr_layout %>% arrange(plate, well)

  # Print information
  print(paste("Selected layout created using a seed of:", omixr_layout$seed[1]))
  print("To recreate this layout, please use the omixr_specific() function with this seed.")

  # Visualize correlations
  print(ggplot(corr_select, aes(x = rand_var, y = tech_var)) +
    geom_tile(aes(fill = corr_val),
              size = 3, colour = "white", show.legend = FALSE) +
    geom_text(aes(label = round(corr_val,3)),
              colour = ifelse(corr_select$corr_val < mean(corr_select$corr_val), "white", "grey30"), fontface = "bold", nudge_y = 0.2, size = 5) +
    geom_text(aes(label = paste("p =",round(corr_p, 3))),
              nudge_y = -0.2, colour = ifelse(corr_select$corr_val < mean(corr_select$corr_val), "white", "grey30")) +
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



