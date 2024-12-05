library(dplyr)
library(tidyr)
library(rlang)


# Function to create dynamic bins

# Apply dynamic binning to each cluster and gender to ensure that there is
#always more than 5 and less than 10 in each bin

#The breaks vector will always have one more element than the counts vector
##because it includes the start and end points of the bins.
create_dynamic_bins <- function(ages, min_size = 5, max_size = 10) {
  n <- length(ages)
  if (n < min_size) {
    return(list(breaks = c(min(ages), max(ages)), counts = c(n)))
  }

  sorted_ages <- sort(ages)
  bins <- list()
  current_bin <- c(sorted_ages[1])

  for (i in 2:n) {
    if (length(current_bin) < max_size) {
      current_bin <- c(current_bin, sorted_ages[i])
    }
    if (length(current_bin) >= min_size) {
      bins[[length(bins) + 1]] <- current_bin
      current_bin <- c()
    }
  }

  if (length(current_bin) > 0) {
    bins[[length(bins) + 1]] <- current_bin
  }

  # Merge bins if any bin has less than min_size
  while (any(sapply(bins, length) < min_size)) {
    small_bins <- which(sapply(bins, length) < min_size)
    for (bin_idx in small_bins) {
      if (bin_idx == 1) {
        bins[[2]] <- c(bins[[1]], bins[[2]])
        bins <- bins[-1]
      } else if (bin_idx == length(bins)) {
        bins[[length(bins) - 1]] <- c(bins[[length(bins) - 1]], bins[[length(bins)]])
        bins <- bins[-length(bins)]
      } else {
        bins[[bin_idx - 1]] <- c(bins[[bin_idx - 1]], bins[[bin_idx]])
        bins[[bin_idx + 1]] <- c(bins[[bin_idx + 1]], bins[[bin_idx]])
        bins <- bins[-bin_idx]
      }
    }
  }

  breaks <- c(min(ages))
  counts <- integer(length(bins))

  for (i in 1:length(bins)) {
    breaks <- c(breaks, max(bins[[i]]))
    counts[i] <- length(bins[[i]])
  }

  return(list(breaks = breaks, counts = counts))
}




# Function calculate age and gender distribution for each cluster
cluster_analysis <- function(cluster_data, X, age_col, gender_col, cluster_col, diagnosis_prefix) {
  # Calculate age and gender distribution for each cluster
  cluster_summary <- cluster_data %>%
    group_by(!!sym(cluster_col)) %>%
    summarise(
      mean_age = mean(!!sym(age_col), na.rm = TRUE),
      age_sd = sd(!!sym(age_col), na.rm = TRUE),
      gender_distribution = list(table(!!sym(gender_col)))
    )

  # Merge X with cluster assignments
  X_cluster_df <- cbind(X, cluster = cluster_data[[cluster_col]])
  X_cluster_df <- as.data.frame(X_cluster_df)

  # Calculate the number of patients with each diagnosis in each cluster
  diagnosis_counts <- X_cluster_df %>%
    group_by(cluster) %>%
    summarise(across(starts_with(diagnosis_prefix), ~sum(. == 1, na.rm = TRUE)))

  # Calculate the total number of patients with each diagnosis across all clusters
  total_diagnosis_counts <- colSums(X)

  ## Normalize the diagnosis counts to get probabilities
  #diagnosis_probability <- diagnosis_counts %>%
  #  mutate(across(starts_with(diagnosis_prefix), ~ . / total_diagnosis_counts[cur_column()]))

  # Combine cluster summary and diagnosis probabilities
  cluster_char <- bind_cols(cluster_summary, diagnosis_counts, .name_repair = "unique")

  return(cluster_char)
}

# Example usage:
# cluster_char <- cluster_analysis(cluster_data, X, "AGE", "GENDER",""cluster, "C0")




create_hex_data_universal <- function(cluster_data, num_bins, x_col, y_col, age_col, gender_col, cluster_col, features_col) {
  # Create bins for x and y coordinates
  bin_width_x <- (max(cluster_data[[x_col]]) - min(cluster_data[[x_col]])) / num_bins
  bin_width_y <- (max(cluster_data[[y_col]]) - min(cluster_data[[y_col]])) / num_bins

  cluster_data <- cluster_data %>%
    mutate(
      x_bin = floor((!!sym(x_col) - min(!!sym(x_col))) / bin_width_x),
      y_bin = floor((!!sym(y_col) - min(!!sym(y_col))) / bin_width_y),
      x_center = min(!!sym(x_col)) + (x_bin + 0.5) * bin_width_x,
      y_center = min(!!sym(y_col)) + (y_bin + 0.5) * bin_width_y
    )

  hex_data <- cluster_data %>%
    group_by(x_bin, y_bin) %>%
    summarise(
      count = n(),
      features = list(unlist(!!sym(features_col))),
      x = first(x_center),
      y = first(y_center),
      gender_counts = list(table(!!sym(gender_col))),  # Count each GENDER in the bin
      mean_age = round(mean(!!sym(age_col), na.rm = TRUE), 0),  # Mean age in the bin
      sd_age = sd(!!sym(age_col), na.rm = TRUE),  # Standard deviation of age in the bin
      most_common_cluster = names(sort(table(!!sym(cluster_col)), decreasing = TRUE))[1]  # Most common cluster in the bin
    ) %>%
    ungroup()

  # Dynamically create columns for each feature
  features <- unique(unlist(hex_data$features))
  for (feature in features) {
    hex_data <- hex_data %>%
      mutate(!!paste0("feature_", feature, "_count") := sapply(features, function(f) sum(f == feature)))
  }

  ## Manually force entries with less than 5 to be 5 and set sd to 0
  hex_data$sd_age[hex_data$count < 5] <- 0
  hex_data$count[hex_data$count < 5] <- 5

  return(hex_data)
}

# Example usage:
# hex_data <- create_hex_data_universal(cluster_data, 30, "V1", "V2", "AGE", "GENDER", "cluster", "features")


