read_and_process_data <- function(file_path) {
  # Read the CSV file
  data <- read.csv(file_path, header = TRUE, sep = ",", row.names = NULL)

  # Convert the data to a matrix and exclude the first column
  data_matrix <- as.matrix(data)


  return(data_matrix)
}



