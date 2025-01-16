age_gender_df <- X %>% select(AGE, GENDER,TARGET_COHORT)
rownames(age_gender_df) <- row_names
X <- X %>% dplyr::select(-c(AGE,GENDER,TARGET_COHORT))
# create_clustering_results(cdm,source_file_path,output_path)



create_clustering_results <- function(cdm,source_file_path,output_path){
  
  max_rows_in <- 5000
  
  # collect the one-hot-encoded data
  dataX <- tbl(cdmCon(cdm), "one_hot_encoded_data") %>% 
    filter(row_number()<1001) %>%  # XX!! take this one out when done debugging and ready for sending to EHDEN
    select(-subject_id) %>% 
    collect() %>% 
    mutate(indexnumber=row_number()) %>% 
    as.data.frame()

  # NB: XX!! remove this ! (or introduce a cap here due to processing times!)
  dataX <- dataX[c(1:min(c(max_rows_in,dim(dataX)[1]))),] 

  row_names <- dataX[,"indexnumber"] # Keep the index for row names
  
  age_gender_df <- dataX %>% select(AGE, GENDER, TARGET_COHORT) # keep age and gender covariates
  rownames(age_gender_df) <- row_names
  
  # remove additional covariates and index from
  dataX <- dataX %>% dplyr::select(-c(AGE,GENDER, TARGET_COHORT, indexnumber))
  
  col_names <- colnames(dataX) # Keep the column names without covariates

  X <- t(dataX) 
  
  X <- apply(X, c(1, 2), as.numeric) 
  
  # Assign the column names back to the matrix
  colnames(X) <- row_names
  
  relevant_features <- apply(X, 2, function(col) col_names[col == 1])
  
  Xt <- t(X)

  # Archetypal Analysis (AA)
  output_path_AA <- file.path(output_path,"Archetypal")
  dir.create(output_path_AA,showWarnings = FALSE)
    
  # Define the range of n_arc values and the number of repetitions NB: XX! Remember to include more arcs and reps! 
  #n_arc_values <- c(5,10,15,20,25,30,35,40,45,50)
  #n_reps <- c(1,2,3,4,5,6,7,8,9,10)
  n_arc_values <- c(2,3,4)
  n_reps <- c(1:10)
  
  pairs <- combn(1:max(n_reps), 2) # get all pairs of the repetitions to calculate MNI
  NMI <- matrix(nrow = length(n_arc_values), ncol = dim(pairs)[2])
  
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
      json_file <- file.path(output_path_AA,paste0("output_", n_arc, "_", n, ".json"))
      write_json(output_list, json_file)
      hex_data <- create_hex_data(X, res$S,res$Z,age_gender_df, relevant_features, num_bins = 30, selected_feature = NULL)
      json_file <- file.path(output_path_AA,paste0("hex_data_", n_arc, "_", n, ".json"))
      jsonlite::write_json(hex_data, json_file)
      
      
      ## OBS: large_da må kun returneres som plots og skal slettes til sidst (eller laves til bins - somehow)
      large_df <- calculate_manhattan_distances(Xt, res$Z, col_names,age_gender_df)
      arc_distances <- aggregate_distances(large_df$all)
      file_name <-  file.path(output_path_AA,paste0("distances_" ,n_arc, "_", n, ".csv"))
      write.csv(arc_distances, file_name, row.names = FALSE,quote=FALSE)
    }
    ## Calc NMI HERE (return NMI)
    
    idx1 <- pairs[1,]
    idx2 <- pairs[2,]
    for (i in seq_along(idx1)) {
      NMI[j,i] <- calcNMI(S_list[[idx1[i]]],S_list[[idx2[i]]])
    }
    
    j <- j+1
    
  }
  file_name <- file.path(output_path_AA,paste0("NMI.csv"))
  write.csv(NMI, file_name)
  
  
  
  ## LDA
  output_path_LDA <- file.path(output_path,"LDA")
  dir.create(output_path_LDA,showWarnings = FALSE)

  i=1
  length_n_arc_list <- length(n_arc_values)
  length_n_reps <- length(n_reps)
  
  # Create an empty matrix with the specified dimensions
  LL_LDA <- matrix(nrow = length_n_arc_list, ncol = length_n_reps)
  for (n_arc in n_arc_values) {
    for (n in n_reps){
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
      
      json_file <- file.path(output_path_LDA,paste0("hex_data_", n_arc, "_", n, ".json"))
      write_json(hex_data, json_file)
      
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
      
      #######
      
      cluster_charB <- cluster_analysis(umap_dbscan_gamma_proj, Xt, "AGE", "GENDER","cluster" ,"DxB") ## Set to before
      cluster_charB <- cluster_charB %>%
        mutate(across(where(is.list), ~ purrr::map(.x, as.list)))
      
      file_name <-  file.path(output_path_LDA,paste0("cluster_char_Before", n_arc, "_", n, ".json"))
      write_json(cluster_charB, file_name)
      
      cluster_charA <- cluster_analysis(umap_dbscan_gamma_proj, Xt, "AGE", "GENDER","cluster" ,"DxA") ## Set to after
      cluster_charA <- cluster_charA %>%
        mutate(across(where(is.list), ~ purrr::map(.x, as.list)))
      
      file_name <-  file.path(output_path_LDA,paste0("cluster_char_After", n_arc, "_", n, ".csv"))
      write.csv(cluster_charA, file_name, row.names=FALSE, quote=FALSE)
      
      cluster_charT <- cluster_analysis(umap_dbscan_gamma_proj, Xt, "AGE", "GENDER","cluster", "TrA") ## Set to Treatment
      cluster_charT <- cluster_charT %>%
        mutate(across(where(is.list), ~ purrr::map(.x, as.list)))
      
      file_name <-  file.path(output_path_LDA,paste0("cluster_char_Treat", n_arc, "_", n, ".csv"))
      write.csv(cluster_charT, file_name, row.names=FALSE, quote=FALSE)
      
      #######
      
      
      LL_LDA[i,n] <- lda_model@loglikelihood
      # Get the word distributions for each topic
      word_distributions <- as.data.frame(lda_model@beta)
      #word_distributions<-10^t(lda_model@beta)
      
      
      frex_score <- calcfrex(lda_model@beta,wordcounts = colSums(Xt))
      rownames(frex_score) <- rownames(X)
      file_name <-  file.path(output_path_LDA,paste0("frex_score", n_arc, "_", n, ".csv"))
      write.csv(frex_score, file_name, row.names=FALSE, quote=FALSE)
      
      # Ensure the word distributions data frame has proper row names
      colnames(word_distributions) <- rownames(X)
      file_name <-  file.path(output_path_LDA,paste0("word_distributions", n_arc, "_", n, ".csv"))
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
      json_file <- file.path(output_path_LDA,paste0("age_bins_", n_arc, "_", n, ".json"))
      jsonlite::write_json(age_bins, json_file)
      
    }
    i = i+1
    
  }
  
  ## UMAP + HDBScan
  output_path_HDBScan_UMAP <- file.path(output_path,"HDBScan_UMAP")
  dir.create(output_path_HDBScan_UMAP,showWarnings = FALSE)

  results <- calculate_contingency_tables(Xt, eps = 0.5, minPts = 5, n_neighbors = 2, n_components = 2)
  
  json_file <- file.path(output_path_HDBScan_UMAP,paste0("hdbscan_Contingency", n_arc, "_", n, ".json")) ## OBS: Testing ranges?
  jsonlite::write_json(results[[1]], json_file)
  
  json_file <- file.path(output_path_HDBScan_UMAP,paste0("hdbscan_hex_data", n_arc, "_", n, ".json")) ## OBS: Testing ranges?
  jsonlite::write_json(results[[2]], json_file)
  
  json_file <- file.path(output_path_HDBScan_UMAP,paste0("hdbscan_cluster_char_Before", n_arc, "_", n, ".json"),)
  jsonlite::write_json(results[[3]], json_file)
  
  json_file <- file.path(output_path_HDBScan_UMAP,paste0("hdbscan_cluster_char_After", n_arc, "_", n, ".json"))
  jsonlite::write_json(results[[4]], json_file)
  
  
  json_file <- file.path(output_path_HDBScan_UMAP,paste0("hdbscan_cluster_char_Treatment", n_arc, "_", n, ".json"))
  jsonlite::write_json(results[[5]], json_file)
  
  json_file <- file.path(output_path_HDBScan_UMAP,paste0("hdbscan_cluster_char_Age_bins", n_arc, "_", n, ".json"))
  jsonlite::write_json(results[[6]], json_file)
  
  
  
  
  
}