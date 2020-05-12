#' Sample Randomization
#'
#' This function takes a sample list, with defined sample ID,
#' and randomizes it a specified number of times.
#'
#' It has an option to specify paired samples (blocks), which
#' will be kept adjacent to one another but with their order
#' shuffled.
#'
#' Randomized sample layouts are then combined with either a
#' pre-specified or custom plate layout, with the option to
#' mask specified wells.
#'
#' The function then chooses the layout that minimizes
#' correlation between randomization variables and technical
#' covariates, so long as no test of correlation reached
#' significance (p<0.05) for that layout.
#'
#' @param df Sample list
#' @param sample_num Name of sample ID variable
#' @param block_var (Optional) Paired sample identifier
#' @param iter_num (Optional) Number of layouts to generate (default: 10,000)
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
#' @importFrom tidyselect everything
#' @import magrittr
#' @export

omix_rand <- function(df,
                      sample_num = sample_num,
                      block_var = block_var,
                      iter_num = 10000,
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

  # Relabel columns based on input
  data_frame <- df %>%
    select(sample_num = sample_num,
           block_var = block_var,
           everything())

  # Set up permutation variable (block or ID)
  if("block_var" %in% colnames(data_frame)) {
    print("Order of paired samples will be randomized.")
    data_frame$perm_var <- data_frame$block_var
  } else {
    data_frame$perm_var <- data_frame$sample_num
  }

  # Create the specified number of permutations of block or ID
  perm_set <- lapply(1:iter_num, function(x) {
    set.seed(x)
    sample(unique(data_frame$perm_var))
  })

  # Initialize variables
  df_rand_list <- list()
  count <- 1

  # Print information
  print("Randomization in progress. This may take a while.")

  # Shuffle blocks and then randomize the list a specified number of times
  for(i in 1:iter_num) {
    df_rand <- data.frame()
    set.seed(i)
    df_shuffle <- data_frame %>%
      group_by(perm_var) %>%
      slice(sample(1:n())) %>%
      ungroup()
    for(j in 1:length(perm_set[[i]])) {
      df_rand <- bind_rows(df_rand, df_shuffle %>%
                             filter(perm_var==perm_set[[i]][j]))
    }
    df_rand_list[[count]] <- df_rand
    count <- count + 1
  }

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

  # Initialize variables
  sample_layout <- list()
  count <- 1

  # Print information
  print("Combining sample list with plate layout.")

  # Combine randomized sample lists with plate layout
  for(i in 1:length(df_rand_list)){
    sample_layout[[count]] <- cbind(
      df_rand_list[[i]],
      plate_layout_masked)
    count <- count + 1
  }

  # If randomization variables are not specified, then use all columns except ID and block
  if(!exists("rand_vars")){
    rand_vars = colnames(df_rand_list[[1]])
  }

  rand_vars <- rand_vars[!rand_vars %in%
                           c("sample_num", "block_var")]

  # Print information
  print("The following columns will have their correlations with technical covariates minimized: ")
  print(rand_vars)

  # Initialize variables
  corr_matrix <- matrix(nrow = length(sample_layout),
                        ncol = length(rand_vars) *
                          length(tech_vars))
  corr_p <- corr_matrix

  # Print information
  print("Calculating correlations.")

  # For each layout, calculate correlations between randomization variables and technical covariates
  for(i in 1:length(sample_layout)) {
    df <- sample_layout[[i]]
    count <- 1
    for(j in rand_vars) {
      for(k in tech_vars) {
        corr_matrix[i,count] <- omix_corr(df[,j], df[,k])$corr_val
        corr_p[i,count] <- omix_corr(df[,j], df[,k])$corr_p
        count <- count + 1
      }
    }
  }

  # Calculate the absolute sum of correlations and set up row ID
  corr_matrix <- data.frame(corr_matrix)
  corr_matrix$corr_sum <- rowSums(abs(corr_matrix))
  corr_matrix$n <- 1:nrow(corr_matrix)

  # Filter out rows where a correlation test returned a significant (< 0.05) p-value
  for(i in 1:length(sample_layout)){
    corr_matrix$p_test[i] <- any(corr_p[i,] < 0.05)
  }

  # Find the layout with the minimum sum of correlations and no significant correlation test
  corr_select <- corr_matrix %>%
    filter(p_test == FALSE) %>%
    filter(corr_sum == min(corr_sum))

  # Check number of optimized layouts
  if(nrow(corr_select)==0){
    print("No randomized layout failed to provide insufficient evidence of correlation between technical covariates and randomization variables.")
    print("Please increase iteration number or reduce number of randomization and/or technical covariates.")
  } else if(nrow(corr_select > 1)) {
    corr_select <- corr_select[1, ]
    print("Several layouts were equally optmized, so only the first will be selected.")
  } else {
  }

  # Select the layout with minimum correlations and no significant correlations
  final_layout <- sample_layout[[corr_select$n]]

  # Rejoin layout with masked wells
  final_layout <- full_join(final_layout, plate_layout, by = c("well", "plate", "row", "columns", "mask"))
  final_layout <- final_layout %>% arrange(plate, well)

  # Print information
  print(paste("Maximum correlation:", round(max(abs(corr_select[1:(length(rand_vars)*length(tech_vars))])), 4)))

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



