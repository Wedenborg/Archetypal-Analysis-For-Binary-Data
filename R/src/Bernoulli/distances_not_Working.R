library(dplyr)
library(tibble)

calculate_manhattan_distances <- function(X, Z, col_names,age_gender_df) {
  # Define the function to calculate the Manhattan distance
  dist_to_vector_manhattan <- function(vec, mat) {
    dist <- colSums(abs(mat-rep(vec,dim(mat)[2])))

    # Return the distance vector
    return(dist)
  }

  # Get the list of all targets
  targets <- col_names

  # Initialize a list to store distances for all targets
  all_distances <- list()
  # Convert row names to a column in X and ensure X is a data frame
  X <- as.data.frame(X)
  X <- rownames_to_column(X, var = "subject_id")

  # Convert row names to a column in age_gender_df
  age_gender_df <- rownames_to_column(age_gender_df, var = "subject_id")

  # Merge X with age_gender_df to include gender information
  X <- X %>%
    left_join(age_gender_df, by = "subject_id")

  # Print the merged data frame to check the merge
  # Split the data into men, women, and all
  X_men <- X %>% filter(GENDER == "Male")
  X_women <- X %>% filter(GENDER == "Female")



  calculate_distances <- function(X_subset) {


    # Loop through each target
    for (target in targets) {
      # Extract the columns that have the target present
      X_target <- X_subset[X_subset[,target ] == 1,]
      #X_target <- X_subset %>% filter(target == 1)
      ## Assuming X_target is your data frame
      X_target_for_dist <- X_target %>%
        select(-1, -((ncol(X_target)-1):ncol(X_target)))


      # Initialize a matrix to store distances for this target
      target_distances <- matrix(ncol = ncol(Z), nrow = nrow(X_target_for_dist))

      # Loop through each column in Z
      for (i in 1:ncol(Z)) {
        # Ensure Z[, i] is a column vector
        vec <- as.matrix(Z[, i])
        if (is.null(dim(vec))) {
          vec <- matrix(vec, ncol = 1)
        }

        # Calculate distances for each column in X_target

        distances <- dist_to_vector_manhattan(vec, t(X_target_for_dist))

        # Store the distances in the matrix
        target_distances[, i] <- distances
      }

      # Store the matrix in the list
      target_distances <- cbind(target_distances,X_target %>% select(GENDER))
      all_distances[[target]] <- target_distances
    }

    # Create a data frame to store all distances
    large_df <- data.frame()

    # Loop through the list of distances and create a data frame
    for (target in targets) {
      # Create a data frame with the distance values and the target
      df <- data.frame(target = rep(target, nrow(all_distances[[target]])),
                       all_distances[[target]])
      #print(dim(df))
      #print(colnames(df)[1:ncol(Z)+1])
      # Set the column names
      colnames(df)[1:ncol(Z)+1] <- paste0("dist_z", 1:ncol(Z))

      # Append the data frame to the large data frame
      large_df <- rbind(large_df, df)
    }

    # Return the large data frame
    return(large_df)
  }

  # Calculate distances for men, women, and all
  distances_all <- calculate_distances(X)
  distances_men <- calculate_distances(X_men)
  distances_women <- calculate_distances(X_women)


  # Return the large data frames for men, women, and all
  return(list(all = distances_all, men = distances_men, women = distances_women))
}

# Example usage:
# Assuming X, Z, and col_names are already defined
# large_df <- calculate_manhattan_distances(X, Z, col_names,age_gender_df)
# print(large_df)
