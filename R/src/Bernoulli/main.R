# Load the necessary packages and functions
library(parallel)
source("R/read_and_process_data.R")

# Define the file path and read and process the data
file_path <- "files/drug_side_effects.csv"
dataX <- read_and_process_data(file_path)
X <- t(dataX)
X <- X[rowSums(X != 0) >= 20, ]
X <- X[, colSums(X != 0) >= 20]

# Define the range of n_arc values and the number of repetitions
n_arc_values <- c(2, 3)
n_reps <- 2

# Define the number of cores to use for parallel processing
num_cores <- detectCores() - 1

# Create a cluster with the specified number of cores
cl <- makeCluster(num_cores)

# Define a function to run the parallel loop
parallel_loop <- function(n_arc, n_reps) {
  # Import the ClosedFormArchetypalAnalysis function from the source file
  source("R/AABer.R")
  source("R/read_and_process_data.R")
  # Define the file path and read and process the data
  file_path <- "files/drug_side_effects.csv"
  dataX <- read_and_process_data(file_path)
  X <- t(dataX)
  X <- X[rowSums(X != 0) >= 20, ]
  X <- X[, colSums(X != 0) >= 20]



  # Repeat the function call n_reps times
  outputs <- lapply(1:n_reps, function(i) ClosedFormArchetypalAnalysis(X, n_arc))
  print(n_arc)

  # Return the outputs for this n_arc value
  return(outputs)
}

# Run the parallel loop
results <- parLapply(cl, n_arc_values, parallel_loop, n_reps = n_reps)

# Stop the cluster
stopCluster(cl)

# Write the results to a JSON file
jsonlite::write_json(results, "results.json")
