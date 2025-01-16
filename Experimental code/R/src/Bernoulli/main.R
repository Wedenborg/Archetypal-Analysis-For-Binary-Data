# Load the necessary packages and functions
library(parallel)
library(foreach)
library(doParallel)

# source("R/read_and_process_data.R")
#
# # Define the file path and read and process the data
# file_path <- "files/drug_side_effects.csv"
# dataX <- read_and_process_data(file_path)
# X <- t(dataX)
# X <- X[rowSums(X != 0) >= 40, ]
# X <- X[, colSums(X != 0) >= 40]

X <- t(one_hot_encoded_CH)
X[c(1:prod(dim(X)))] <- as.integer(X[c(1:prod(dim(X)))])

# Define the range of n_arc values and the number of repetitions
n_arc_values <- c(2, 3,4)

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)


# Define a function to run the parallel loop
parallel_loop <- function(X,n_arc, n_reps) {
  # Import the ClosedFormArchetypalAnalysis function from the source file
  source("R/AA_AABer.R")
  # source("R/read_and_process_data.R")
  # # Define the file path and read and process the data
  # file_path <- "files/drug_side_effects.csv"
  # dataX <- read_and_process_data(file_path)
  # X <- t(dataX)
  # X <- X[rowSums(X != 0) >= 20, ]
  # X <- X[, colSums(X != 0) >= 20]


  # Repeat the function call n_reps times
  outputs <- lapply(1:n_reps, function(i) ClosedFormArchetypalAnalysis(X, n_arc))
  print(n_arc)

  # Return the outputs for this n_arc value
  return(outputs)
}
# Export the function to the cluster
clusterExport(cl, c("parallel_loop","X"))

# Run the parallel loop
results <- parLapply(cl, n_arc_values, function(arc) parallel_loop(X,arc,n_reps = 2))

# Stop the cluster
stopCluster(cl)

# Write the results to a JSON file
jsonlite::write_json(results, "results.json")
