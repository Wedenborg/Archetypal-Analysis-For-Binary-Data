calculate_contingency_tables <- function(X, eps = 0.5, minPts = 5, n_neighbors = 2, n_components = 2) {
  # Install and load necessary packages
  if (!requireNamespace("umap", quietly = TRUE)) {
    install.packages("umap")
  }
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    install.packages("dbscan")
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    install.packages("reshape2")
  }
  library(umap)
  library(dbscan)
  library(reshape2)
  library(ggplot2)
  library(viridis)
  library(dplyr)
  source("utils.R")
  # Apply UMAP
  umap_result <- umap(X, n_neighbors = n_neighbors, n_components = n_components, metric = "manhattan")
  # Save the cluster labels



  # Apply DBSCAN to the UMAP-reduced data
  dbscan_result <- hdbscan(umap_result$layout,  minPts = minPts)
  cluster_labels <- dbscan_result$cluster
  umap_df <- as.data.frame(umap_result$layout)
  umap_df <- bind_cols(umap_df, age_gender_df)

  cluster_labels_df <- as.data.frame(cluster_labels)
  umap_df <- bind_cols(umap_df,cluster_labels_df)

  # Add the relevant feature names as a new column in proj_h_df
  umap_df$features <- I(relevant_features)  # Use I() to prevent coercion to factor
  hex_data <- create_hex_data_universal(umap_df, 30, "V1", "V2", "AGE", "GENDER", "cluster_labels", "features")

  cluster_charB <- cluster_analysis(umap_df, X, "AGE", "GENDER", "cluster_labels","DxB") ## Set to before
  cluster_charB <- cluster_charB %>%
    mutate(across(where(is.list), ~ purrr::map(.x, as.list)))

  cluster_charA <- cluster_analysis(umap_df, X, "AGE", "GENDER", "cluster_labels","DxA") ## Set to After
  cluster_charA <- cluster_charA %>%
    mutate(across(where(is.list), ~ purrr::map(.x, as.list)))

  cluster_charT <- cluster_analysis(umap_df, X, "AGE", "GENDER", "cluster_labels","TrA") ## Set to Treatment
  cluster_charT <- cluster_charT %>%
    mutate(across(where(is.list), ~ purrr::map(.x, as.list)))

  age_bins <- umap_df %>%
    group_by(cluster_labels, GENDER) %>%
    summarise(
      age_bins = list(create_dynamic_bins(AGE)),
      .groups = 'drop'
    )




  # Get unique cluster labels (excluding noise points, typically labeled as 0)
  unique_clusters <- unique(cluster_labels)
  unique_clusters <- unique_clusters[unique_clusters != 0]


  # Initialize lists to store the results
  contingency_tables <- list()

  # Loop through each cluster
  for (cluster in unique_clusters) {
    # Extract the data for the current cluster

    cluster_data <- X[cluster_labels == cluster, ]

    # Get the number of features
    num_features <- ncol(cluster_data)

    # Initialize a 3D array to store the contingency tables
    contingency_table <- array(0, dim = c(num_features, num_features, 4))

    # Loop through each pair of features
    for (i in 1:(num_features - 1)) {
      for (j in (i + 1):num_features) {
        # Count the occurrences of the four combinations (11, 10, 01, 00)
        counts_11 <- sum(cluster_data[, i] == 1 & cluster_data[, j] == 1)
        counts_10 <- sum(cluster_data[, i] == 1 & cluster_data[, j] == 0)
        counts_01 <- sum(cluster_data[, i] == 0 & cluster_data[, j] == 1)
        counts_00 <- sum(cluster_data[, i] == 0 & cluster_data[, j] == 0)

        # Store the counts in the 3D array
        contingency_table[i, j, 1] <- counts_11
        contingency_table[i, j, 2] <- counts_10
        contingency_table[i, j, 3] <- counts_01
        contingency_table[i, j, 4] <- counts_00

        # Symmetric matrix
        contingency_table[j, i, 1] <- counts_11
        contingency_table[j, i, 2] <- counts_01
        contingency_table[j, i, 3] <- counts_10
        contingency_table[j, i, 4] <- counts_00
      }
    }

    # Set the row and column names for better readability
    rownames(contingency_table) <- colnames(cluster_data)
    colnames(contingency_table) <- colnames(cluster_data)


    # Save the contingency table
    contingency_tables[[as.character(cluster)]] <- contingency_table
  }

  # Return the results as a list
  return(list(contingency_tables, hex_data,cluster_charB,cluster_charA,cluster_charT,age_bins))
}
