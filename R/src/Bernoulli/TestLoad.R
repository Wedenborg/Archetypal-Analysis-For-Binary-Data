source("NMI.R")
library(jsonlite)

# Define the values of n_reps
n_reps <- c(1, 2, 3)

# Initialize an empty list to store the data from each JSON file
data <- list()

# Iterate over the values of n_reps
for (n in n_reps) {
  # Construct the filename using paste0
  json_file <- paste0("../output_4_", n, ".json")

  # Load the JSON file into R
  data[[n]] <- fromJSON(json_file)
}
