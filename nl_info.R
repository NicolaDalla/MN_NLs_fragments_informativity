library(igraph)
library(tidyverse)
library(dplyr)
library(Spectra)
library(qgam)

##function for NLs informativity
#assuming that: msms is our MS/MS spectra in Spectra format; matrix: similaity natrix genrated form msms Spectra;
#perm: the number of permutaiton test for Z-score calcualtion; min_presence and max_presence: number of minimum and maximum node having a spectral features;
#Q and K: QGAM parameters (Q=quantile, K=k)

nl_info <- function(matrix, msms, perm=100, thr, min_presence=10, max_presence=length(msms), Q=0.8, K=10){
  
  #calculatin NLs
  neutral_loss <- function(x, precursorMz, ...) {
    x[, "mz"] <- precursorMz - x[, "mz"]
    x[order(x[, "mz"]), , drop = FALSE]
  }
  nl <- addProcessing(msms, neutral_loss,
                      spectraVariables = "precursorMz")
  msms$nl<- mz(nl)
  
  #listing the most frequent losses 
  nl_list <- msms$nl
  all <- sort(unlist(nl_list))
  
  #grouping the NLs
  db <- dbscan(as.matrix(all), eps = 0.001, minPts = 3)
  
  df <-data.frame(
    value = all,
    cluster = db$cluster
  )
  
  nldf <- df %>%
    filter(cluster != 0) %>%    
    group_by(cluster) %>%
    summarise(
      val = mean(value),
      freq = n()
    ) %>%
    arrange(desc(freq)) %>%
    filter(freq>min_presence & freq< max_presence)
  
  #creating the graph object
  
  create_edge_list <- function(similarity_matrix, thr, spectra) {
    non_zero_indices <- which(similarity_matrix > thr, arr.ind = TRUE)
    if (length(non_zero_indices) > 0 && is.null(dim(non_zero_indices))) {
      non_zero_indices <- matrix(non_zero_indices, ncol = 2, byrow = TRUE)
    }
    if (length(non_zero_indices) == 0) {
      return(data.frame(from = character(), to = character(), weight = numeric()))
    }
    non_zero_indices <- non_zero_indices[non_zero_indices[, 1] < non_zero_indices[, 2], , drop = FALSE]
    if (nrow(non_zero_indices) == 0) {
      return(data.frame(from = character(), to = character(), weight = numeric()))
    }
    edges <- data.frame(
      from = spectra$id[non_zero_indices[, 1]],
      to = spectra$id[non_zero_indices[, 2]],
      weight = similarity_matrix[non_zero_indices]
    )
    return(edges)
  }
  
  edges <- create_edge_list(matrix,thr, msms)
  involved_nodes <- unique(c(edges$from, edges$to))
  idx <- match(involved_nodes, msms$id)
  nodes <- data.frame(
    id = msms$id[idx]
  )
  graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
  
  nl_list <- msms$nl[match(V(graph)$name, msms$id)]
  
  #df for results
  results <- data.frame(
    nli = nldf$val, 
    number = NA_real_,
    z_score = NA_real_
  )
  
  #calculating unsuipervized community
  cl <- cluster_louvain(graph)
  V(graph)$unsupervised <- membership(cl)
  
  #calculating metrics fo each NL
  for (i in seq_len(nrow(results))) {
    nli <- results$nli[i]
    has_nl <- vapply(seq_along(nl_list), function(j) {
      vals <- as.numeric(as.character(nl_list[[j]]))
      any(abs(vals - nli) <= 0.005, na.rm = TRUE)
    }, logical(1))
    
    n_group <- sum(has_nl)
    results$number[i] <- n_group
    
    if (sum(n_group) >= min_presence & sum(n_group) <= max_presence) {
      
      #MODULARITY z-score
      V(graph)$comm <- ifelse(has_nl, 1, 2)
      mod_real <- modularity(graph, membership = V(graph)$comm)
      
      rand_mod <- replicate(perm, {
        rand_nodes <- sample(vcount(graph), sum(has_nl))
        comm <- rep(2, vcount(graph))
        comm[rand_nodes] <- 1
        modularity(graph, membership = comm)
      })
      
      if (is.na(mod_real) || sd(rand_mod, na.rm = TRUE) == 0) {
        z <- NA
      } else {
        
        z <- (mod_real - abs(mean(rand_mod, na.rm = TRUE))) / sd(rand_mod, na.rm = TRUE)
        
      }
      
      results$z_score[i] <- z
      
    } 
  }
  results<- results[!is.na(results$z_score),]
  
  #compute QGAM model and calculate the residuals
  results$log10 <- log10(results$number)
  qg <- qgam(
    z_score ~ s(log10, k = K),
    data = results,
    qu = Q
  )
  
  q_pred <- predict(qg)
  results$excess <- results$z_score - q_pred
  results <- results %>%
    select(-log10)
  return(results)
}

#plotting the top results(according to QGAM residuals)
#assuming rs is the data frame generated from 'nl_info' function, variable needs to be the same used in precedent function
#n: number of the top spectral features to plot

plot_nl <- function(rs, msms, matrix, thr, n){
  #neutrla losses 
  #calculatin NLs
  neutral_loss <- function(x, precursorMz, ...) {
    # Compute neutral loss
    x[, "mz"] <- precursorMz - x[, "mz"]
    
    # Sort by neutral loss m/z
    x[order(x[, "mz"]), , drop = FALSE]
  }
  nl <- addProcessing(msms, neutral_loss,
                      spectraVariables = "precursorMz")
  msms$nl<- mz(nl)
  
  nl_list <- msms$nl
  #base plot 
  edges <- create_edge_list(matrix,thr, msms)
  
  involved_nodes <- unique(c(edges$from, edges$to))
  idx <- match(involved_nodes, msms$id)
  
  nodes <- data.frame(
    id = msms$id[idx]
  )
  graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
  
  top_nl <- rs %>%
    arrange(desc(excess)) %>%
    slice_head(n=n)
  
  #plottign unsupervized
  cl <- cluster_louvain(graph)
  V(graph)$unsupervised <- membership(cl)
  lay <- layout_with_fr(graph)
  
  plot(graph, 
       vertex.color = V(graph)$unsupervised, 
       vertex.label = NA,                
       vertex.size = 5, 
       layout = lay,  
       main = "Graph colored by Louvain clusters")
  #plotting nls
  nl_list <- msms$nl[match(V(graph)$name, msms$id)]
  
  for (x in 1:n) {
    nli<- top_nl$nli[x]
    # find which vertices contain this neutral loss
    has_nl <- vapply(seq_along(nl_list), function(j) {
      vals <- as.numeric(as.character(nl_list[[j]]))
      any(abs(vals - nli) <= 0.005, na.rm = TRUE)
    }, logical(1))
    
    # color vertices: red if contains NL, gray otherwise
    V(graph)$color <- ifelse(has_nl, "red",  rgb(135, 206, 235, 150, maxColorValue = 255))
    
    # plot separately
    plot(
      graph,
      vertex.size = 6,
      layout = lay,    
      vertex.label = NA,
      main = paste("Neutral loss:", round(nli,4))
    )
  }
}















