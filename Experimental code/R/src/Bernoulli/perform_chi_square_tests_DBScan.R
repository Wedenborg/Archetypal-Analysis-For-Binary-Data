perform_chi_square_tests <- function(contingency_tables) {
  # Initialize a list to store the Chi-square results
  chi_square_results <- list()

  # Loop through each cluster
  for (cluster in names(contingency_tables)) {
    # Get the contingency table for the current cluster
    contingency_table <- contingency_tables[[cluster]]

    # Get the number of features
    num_features <- nrow(contingency_table)

    # Initialize a matrix to store the p-values
    p_values <- matrix(NA, nrow = num_features, ncol = num_features)

    # Loop through each pair of features
    for (i in 1:(num_features - 1)) {
      for (j in (i + 1):num_features) {
        # Extract the counts for the current pair of features
        counts_11 <- contingency_table[i, j, 1]
        counts_10 <- contingency_table[i, j, 2]
        counts_01 <- contingency_table[i, j, 3]
        counts_00 <- contingency_table[i, j, 4]

        # Create a 2x2 contingency table for the Chi-square test
        contingency_2x2 <- matrix(c(counts_11, counts_10, counts_01, counts_00), nrow = 2, byrow = TRUE)

        # Check for zero counts or insufficient data
        if (any(contingency_2x2 == 0) || sum(contingency_2x2) < 5) {
          p_values[i, j] <- NA  # Set p-value to NA if conditions are not met
          p_values[j, i] <- NA  # Symmetric matrix
        } else {
          # Perform the Chi-square test
          chi_square_test <- chisq.test(contingency_2x2)

          # Store the p-value in the matrix
          p_values[i, j] <- chi_square_test$p.value
          p_values[j, i] <- chi_square_test$p.value  # Symmetric matrix
        }
      }
    }

    # Set the row and column names for better readability
    rownames(p_values) <- rownames(contingency_table)
    colnames(p_values) <- colnames(contingency_table)

    # Save the p-values in the results list
    chi_square_results[[cluster]] <- p_values
  }

  # Return the Chi-square results
  return(chi_square_results)
}
