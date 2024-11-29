
# Define the function to create the plot
plot_distances <- function(large_df, target1, target2) {
  # Filter the data frame for the specified targets
  filtered_df <- large_df %>%
    filter(category %in% c(target1, target2))
  
  # Create the plot
  p <- ggplot(filtered_df, aes(x = category, y = distance, fill = category)) +
    geom_boxplot() +
    scale_fill_manual(values = c(target1 = "blue", target2 = "red")) +
    theme_minimal() +
    labs(title = paste("Distance Comparison between", target1, "and", target2),
         x = "Target",
         y = "Distance") +
    theme(axis.title = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = 0.5))
  
  # Print the plot
  print(p)
}