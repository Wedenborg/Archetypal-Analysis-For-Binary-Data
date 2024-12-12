# Load necessary packages
library(dplyr)
library(tidyr)
library(archetypes)

# Define the function
create_hex_data <- function(X, S, Z,age_gender_df, relevant_features, num_bins = 30, selected_feature = NULL) {



  params <- t(Z)  #z
  coefs <- t(S) # t(ae) #alphas

  order = NULL
  labels_cex = 1
  labels = NULL
  show_labels = TRUE
  points_col = "#00000044"
  points_pch = 19
  points_cex = 1
  projection = simplex_projection
  show_points = TRUE
  show_circle = TRUE
  circle_col = "lightgray"
  show_edges = TRUE
  edges_col = "lightgray"
  show_direction = TRUE
  direction_length = 1
  directions_col = points_col

  object=c()
  object$alphas = coefs
  radius <- 10
  proj_z <- projection(params, r = radius - 1)

  proj_h <- coefs %*% proj_z

  # Extract the relevant feature names for each sample


  # Convert proj_h to a data frame
  proj_h_df <- as.data.frame(proj_h)
  proj_h_df <- bind_cols(proj_h_df, age_gender_df)

  # Add the relevant feature names as a new column in proj_h_df
  proj_h_df$features <- I(relevant_features)  # Use I() to prevent coercion to factor

  # Create bins for x and y coordinates
  bin_width_x <- (max(proj_h_df$x) - min(proj_h_df$x)) / num_bins
  bin_width_y <- (max(proj_h_df$y) - min(proj_h_df$y)) / num_bins

  proj_h_df <- proj_h_df %>%
    mutate(
      x_bin = floor((x - min(x)) / bin_width_x),
      y_bin = floor((y - min(y)) / bin_width_y),
      x_center = min(x) + (x_bin+0.5) * bin_width_x,
      y_center = min(y) + (y_bin+0.5) * bin_width_y
    )

  hex_data <- proj_h_df %>%
    group_by(x_bin, y_bin) %>%
    summarise(
      count = n(),
      features = list(unlist(features)),
      x = first(x_center),
      y = first(y_center),
      gender_max = names(which.max(table(GENDER))),  # Count each GENDER in the bin
      cohort_max = names(which.max(table(TARGET_COHORT))),
      mean_age = round(mean(AGE, na.rm = TRUE),0),  # Mean age in the bin
      sd_age = sd(AGE, na.rm = TRUE)  # Standard deviation of age in the bin
    ,.groups = 'drop')


  ## manually force entries with less than 5 to be 5 and set sd to 0.
  hex_data$sd_age[hex_data$count<5] = 0
  hex_data$count[hex_data$count<5] = 5


  # Dynamically create columns for each feature
  features <- unique(unlist(hex_data$features))
  for (feature in features) {
    hex_data <- hex_data %>%
      mutate(!!paste0("feature_", feature, "_count") := sapply(features, function(f) sum(f == feature)))
  }


  # Select a specific feature to visualize if not provided
  #if (is.null(selected_feature)) {
  #  selected_feature <- feature_names[1]
  #}

  # Create a binary indicator for the presence of the selected feature
  #hex_data <- hex_data %>%
  #  mutate(feature_present = sapply(features, function(f) any(f == selected_feature)))

  return(hex_data)
}

