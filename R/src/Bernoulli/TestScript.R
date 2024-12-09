# Load necessary libraries
library(topicmodels)
library(tm)
library(stm)
library(wordcloud)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(umap)
library(dbscan)
library(cluster)
library(ggplot2)
library(jsonlite)
library(umap)
library(sbm)
source("utils.R")
source("read_and_process_data.R")
source("AABer.R")
source("NMI.R")
source("create_hex_data.R")
source("distances_not_Working.R")
source("calculate_contingency_tables.R")

# Define the file path and read and process the data
file_path <- "../files/drug_side_effects.csv"

# Read and process the data
dataX <- read_and_process_data(file_path)



# Filter the matrix
X <- dataX[rowSums(dataX != 0) >= 50, ]
X <- X[, colSums(dataX != 0) >= 50]



X <- as.data.frame(X)

X <- X[1:100,1:150 ]

# Keep the column names
col_names <- colnames(X)
row_names <- X[,1] ### MÅSKE Denne her alt efter hvad row_names er??


col_names <- col_names[col_names != "drug_name"]
fake_ages <- sample(18:65, nrow(X), replace = TRUE)
fake_genders <- sample(c("Male", "Female"), nrow(X), replace = TRUE)
X$AGE <- fake_ages
X$GENDER <- fake_genders


age_gender_df <- X %>% select(AGE, GENDER)
rownames(age_gender_df) <- row_names
X <- X %>% dplyr::select(-c(AGE,GENDER))

relevant_features <- apply(X, 2, function(col) row_names[col == 1])

X <- t(X[, -1]) ## Vær lige obs på -1 her. Dette kommer an på hvordan dataen ser ud.


X <- apply(X, c(1, 2), as.numeric) ## Tror dette er så den er numerisk og ikke T/F


# Assign the column names back to the matrix
colnames(X) <- row_names

relevant_features <- apply(X, 2, function(col) col_names[col == 1])

Xt <- t(X)

# Define the range of n_arc values and the number of repetitions
#n_arc_values <- c(5,10,15,20,25,30,35,40,45,50)
n_arc_values <- c(2,3,4)
#n_reps <- c(1,2,3,4,5,6,7,8,9,10)
n_reps = c(1,2,3)
pairs <- combn(1:max(n_reps), 2)
NMI <- matrix(ncol = length(n_arc_values), nrow = dim(pairs)[2])
## AA
dir.create('AA_results', showWarnings = FALSE)

j = 1
S_list <- list()
for (n_arc in n_arc_values) {
  for (n in n_reps){
    res <- ClosedFormArchetypalAnalysis(X, n_arc)
    # Create a list to store the output for this iteration
    output_list <- list(
      n_arc = n_arc,
      n = n,
      res = res[-which(names(res)=="S")]
    )
    S_list[[n]]<- res$S
    # Write the output list to a JSON file
    ## OBS: S skal ikke gemmes i .json
    json_file <- paste0("AA_results/output_", n_arc, "_", n, ".json")
    write_json(output_list, json_file)
    hex_data <- create_hex_data(X, res$S,res$Z,age_gender_df, relevant_features, num_bins = 30, selected_feature = NULL)
    json_file <- paste0("AA_results/hex_data_", n_arc, "_", n, ".json")
    write_json(hex_data, json_file)

    ## OBS: large_da må kun returneres som plots og skal slettes til sidst (eller laves til bins - somehow)
    large_df <- calculate_manhattan_distances(Xt, res$Z, col_names,age_gender_df)
    for (name in names(large_df)) {
      file_name <-  paste0("AA_results/distances_",name,"_" ,n_arc, "_", n, ".csv")
      write.csv(large_df[[name]], file_name, row.names = FALSE,quote=FALSE)
    }
  }
  ## Calc NMI HERE (return NMI)

  idx1 <- pairs[1,]
  idx2 <- pairs[2,]
  for (i in seq_along(idx1)) {
    NMI[j,i] <- calcNMI(S_list[[idx1[i]]],S_list[[idx2[i]]])
  }

  j <- j+1

}
file_name <- paste0("AA_results/NMI.csv")
write.csv(NMI, file_name)


dir.create('LDA_results', showWarnings = FALSE)
## LDA
i=1
length_n_arc_list <- length(n_arc_values)
length_n_reps <- length(n_reps)

# Create an empty matrix with the specified dimensions
LL_LDA <- matrix(nrow = length_n_arc_list, ncol = length_n_reps)
for (n_arc in n_arc_values) {
  for (n in n_reps){
    print('LDA')

    lda_model <- LDA(Xt, n_arc, method = "Gibbs")

    gamma <- lda_model@gamma
    dbscan_gamma <- hdbscan(gamma,  minPts = 10)
    #umap_dbscan_gamma <- umap(t(X), n_neighbors = 10, n_components = 2, metric = "cosine")
    umap_dbscan_gamma <- umap(gamma, n_neighbors = 10, n_components = 2, metric = "cosine")
    umap_dbscan_gamma_proj <- as.data.frame(umap_dbscan_gamma$layout)
    umap_dbscan_gamma_proj$cluster <- dbscan_gamma$cluster
    umap_dbscan_gamma_proj <- bind_cols(umap_dbscan_gamma_proj, age_gender_df)
    umap_dbscan_gamma_proj$features <- I(relevant_features)


    hex_data <- create_hex_data_universal(umap_dbscan_gamma_proj, 30, "V1", "V2", "AGE", "GENDER", "cluster", "features")

    json_file <- paste0("LDA_results/hex_data_", n_arc, "_", n, ".json")
    write_json(hex_data, json_file)


    #ggplot(hex_data, aes(x = x_bin, y = y_bin, fill = most_common_cluster)) +
    #  geom_hex()

    # Plot the UMAP results using ggplot2
    #ggplot(umap_dbscan_gamma_proj, aes(x = V1, y = V2, color = as.factor(cluster))) +
    #  geom_point(size = 2) +
    #  scale_color_manual(values = c("grey","red", "blue", "green")) +
    #  labs(title = "2D UMAP Plot", x = "UMAP1", y = "UMAP2", color = "Cluster") +
    #  theme_minimal()


    # Calculate age and gender distribution for each cluster
    cluster_summary <- umap_dbscan_gamma_proj %>%
      group_by(cluster) %>%
      summarise(
        mean_age = mean(AGE, na.rm = TRUE),
        age_sd = sd(AGE, na.rm = TRUE),
        gender_distribution = list(table(GENDER))
      )

    # Merge X_df with cluster assignments
    X_cluster_df <- cbind(t(X), cluster = umap_dbscan_gamma_proj$cluster)
    X_cluster_df <- as.data.frame(X_cluster_df)

    ####### OBS HER: SKAL RETTES OG KOPIERES! (noget med ends_with "Before"/"After"/"Treatment")

    cluster_charB <- cluster_analysis(umap_dbscan_gamma_proj, Xt, "AGE", "GENDER","cluster" ,"DxB") ## Set to before
    cluster_charB <- cluster_charB %>%
      mutate(across(where(is.list), ~ purrr::map(.x, as.list)))

    file_name <-  paste0("LDA_results/cluster_char_Before", n_arc, "_", n, ".json")
    write_json(cluster_charB, file_name)

    cluster_charA <- cluster_analysis(umap_dbscan_gamma_proj, Xt, "AGE", "GENDER","cluster" ,"DxA") ## Set to after
    cluster_charA <- cluster_charA %>%
      mutate(across(where(is.list), ~ purrr::map(.x, as.list)))

    file_name <-  paste0("LDA_results/cluster_char_After", n_arc, "_", n, ".csv")
    write.csv(cluster_charA, file_name, row.names=FALSE, quote=FALSE)

    cluster_charT <- cluster_analysis(umap_dbscan_gamma_proj, Xt, "AGE", "GENDER","cluster", "TrA") ## Set to Treatment
    cluster_charT <- cluster_charT %>%
      mutate(across(where(is.list), ~ purrr::map(.x, as.list)))

    file_name <-  paste0("LDA_results/cluster_char_Treat", n_arc, "_", n, ".csv")
    write.csv(cluster_charT, file_name, row.names=FALSE, quote=FALSE)

    #######


    LL_LDA[i,n] <- lda_model@loglikelihood
    # Get the word distributions for each topic
    word_distributions <- as.data.frame(lda_model@beta)
    #word_distributions<-10^t(lda_model@beta)


    frex_score <- calcfrex(lda_model@beta,wordcounts = colSums(Xt))
    rownames(frex_score) <- rownames(X)
    file_name <-  paste0("LDA_results/frex_score", n_arc, "_", n, ".csv")
    write.csv(frex_score, file_name, row.names=FALSE, quote=FALSE)

    # Ensure the word distributions data frame has proper row names
    colnames(word_distributions) <- rownames(X)
    file_name <-  paste0("LDA_results/word_distributions", n_arc, "_", n, ".csv")
    write.csv(word_distributions, file_name, row.names=FALSE, quote=FALSE)

    # Apply dynamic binning to each cluster and gender to ensure that there is
    #always more than 5 and less than 10 in each bin

    #The breaks vector will always have one more element than the counts vector
    ##because it includes the start and end points of the bins.

    age_bins <- umap_dbscan_gamma_proj %>%
      group_by(cluster, GENDER) %>%
      summarise(
        age_bins = list(create_dynamic_bins(AGE)),
        .groups = 'drop'
      )
    json_file <- paste0("LDA_results/age_bins_", n_arc, "_", n, ".json")
    write_json(age_bins, json_file)

  }
  i = i+1

}


##########################TEST

## UMAP + HDBScan

print('UMAP')

dir.create('HDBScan_results', showWarnings = FALSE)
results <- calculate_contingency_tables(Xt, eps = 0.5, minPts = 5, n_neighbors = 2, n_components = 2)

json_file <- paste0("HDBScan_results/hdbscan_Contingency", n_arc, "_", n, ".json") ## OBS: Testing ranges?
write_json(results[[1]], json_file)

json_file <- paste0("HDBScan_results/hdbscan_hex_data", n_arc, "_", n, ".json") ## OBS: Testing ranges?
write_json(results[[2]], json_file)

json_file <- paste0("HDBScan_results/hdbscan_cluster_char_Before", n_arc, "_", n, ".json")
write_json(results[[3]], json_file)

json_file <- paste0("HDBScan_results/hdbscan_cluster_char_After", n_arc, "_", n, ".json")
write_json(results[[4]], json_file)


json_file <- paste0("HDBScan_results/hdbscan_cluster_char_Treatment", n_arc, "_", n, ".json")
write_json(results[[5]], json_file)

json_file <- paste0("HDBScan_results/hdbscan_cluster_char_Age_bins", n_arc, "_", n, ".json")
write_json(results[[6]], json_file)























###################

# Assuming your data matrix is named 'data' and you have a target variable 'target'
# data <- your_binary_data_matrix
# target <- your_categorical_target_variable

# Apply UMAP
umap_result <- umap(data, n_neighbors = n_arc, n_components = 2)

# Apply DBSCAN to the UMAP-reduced data
dbscan_result <- dbscan(umap_result$layout, eps = 0.5, minPts = 5)

# Save the cluster labels
cluster_labels <- dbscan_result$cluster

# Save the core points
core_points <- dbscan_result$core_sample_indices

# Save the noise points
noise_points <- which(cluster_labels == 0)  # Assuming noise points are labeled as 0

# Optionally, calculate and save the cluster centers
calculate_cluster_centers <- function(data, labels) {
  unique_labels <- unique(labels)
  unique_labels <- unique_labels[unique_labels != 0]  # Exclude noise points
  centers <- sapply(unique_labels, function(label) {
    colMeans(data[labels == label, , drop = FALSE])
  }, simplify = FALSE)
  do.call(rbind, centers)
}

cluster_centers <- calculate_cluster_centers(umap_result$layout, cluster_labels)

# Optionally, calculate and save the silhouette score
if (length(unique(cluster_labels)) > 1) {
  sil_result <- silhouette(cluster_labels, dist(umap_result$layout))
  silhouette_score <- mean(sil_result[, "sil_width"])
} else {
  silhouette_score <- NA  # Silhouette score is not meaningful if there is only one cluster or all points are noise
}

# Calculate the number of samples in each cluster
cluster_sizes <- table(cluster_labels)

# Create a contingency table for the Chi-square test
contingency_table <- table(cluster_labels, target)

# Perform the Chi-square test
chi_square_result <- chisq.test(contingency_table)

# Save the results to a list
dbscan_results <- list(
  cluster_labels = cluster_labels,
  core_points = core_points,
  noise_points = noise_points,
  cluster_centers = cluster_centers,
  silhouette_score = silhouette_score,
  cluster_sizes = cluster_sizes,
  chi_square_result = chi_square_result
)


#############################################
#SBM

