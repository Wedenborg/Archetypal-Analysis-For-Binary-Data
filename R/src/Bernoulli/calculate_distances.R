library(dplyr)

calculate_manhattan_distances <- function(X, Z, col_names,age_gender_df) {
  # Define the function to calculate the Manhattan distance
  dist_to_vector_manhattan <- function(vec, mat) {
    # Ensure vec is a column vector
    vec <- as.matrix(vec)
    if (is.null(dim(vec))) {
      vec <- matrix(vec, ncol = 1)
    }

    vec <- array(vec)

    # Ensure mat is a matrix
    mat <- as.matrix(mat)

    # Check if the number of rows in vec matches the number of rows in mat
    if (nrow(vec) != nrow(mat)) {
      stop("The number of rows in vec does not match the number of rows in mat.")
    }

    # Calculate the Manhattan distance between each column of mat and vec
    dist <- colSums(abs(sweep(mat, 2, vec, "-")))

    # Return the distance vector
    return(dist)
  }

  # Get the list of all targets
  targets <- col_names

  # Initialize a list to store distances for all targets
  all_distances <- list()

  # Loop through each target
  for (target in targets) {
    # Extract the columns that have the target present
    X_target <- X[, X[target, ] == 1, drop = FALSE]


    # Initialize a matrix to store distances for this target
    target_distances <- matrix(ncol = ncol(Z), nrow = ncol(X_target))

    # Loop through each column in Z
    for (i in 1:ncol(Z)) {
      # Ensure Z[, i] is a column vector
      vec <- as.matrix(Z[, i])
      if (is.null(dim(vec))) {
        vec <- matrix(vec, ncol = 1)
      }

      # Ensure X_target has the correct dimensions
      if (nrow(X_target) != nrow(vec)) {
        stop(paste("Dimension mismatch for target:", target, "and column:", i))
      }

      # Calculate distances for each column in X_target
      distances <- dist_to_vector_manhattan(vec, X_target)

      # Store the distances in the matrix
      target_distances[, i] <- distances
    }

    # Store the matrix in the list
    all_distances[[target]] <- target_distances
  }

  # Create a data frame to store all distances
  large_df <- data.frame()

  # Loop through the list of distances and create a data frame
  for (target in targets) {
    # Create a data frame with the distance values and the target
    df <- data.frame(target = rep(target, nrow(all_distances[[target]])),
                     all_distances[[target]])

    # Set the column names
    colnames(df)[-1] <- paste0("dist_z", 1:ncol(Z))

    # Append the data frame to the large data frame
    large_df <- rbind(large_df, df)
  }

  # Return the large data frame
  return(large_df)
}

# Example usage:
# Assuming X, Z, and col_names are already defined
# large_df <- calculate_manhattan_distances(X, Z, col_names)
# print(large_df)
