library(igraph)
library(tidyverse)
library(dplyr)
library(Spectra)
library(qgam)

##function for MN spectral features informativity 
#we can investigate neutral losses or fragment ions, according to "method" parameters which can be set to "FIs" or "NLs"
#assuming that: msms is our MS/MS spectra in Spectra format; matrix: similarity matrix generated form msms Spectra;
#perm: the number of permutation test for Z-score calculation; min_presence and max_presence: number of minimum and maximum node having a spectral features;
#Q and K: QGAM parameters (Q=quantile, K=k)
#eps and minPts: DBscan parameters, where: eps: the maximum distance betweeen 2 points to be considered of the same group; 
#                                          minPts: minimum indifiduals required to create a group


MN_info <- function(method, matrix, msms, perm=100, thr, min_presence=10, max_presence=length(msms), Q=0.5, K=5, eps=0.001, minPts= 3){
  #for NLs investigation
  if(method=="NLs"){
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
    db <- dbscan(as.matrix(all), eps = eps, minPts = minPts)
    
    df <-data.frame(
      value = all,
      cluster = db$cluster
    )
    
    df <- df %>%
      filter(cluster != 0) %>%    
      group_by(cluster) %>%
      summarise(
        val = mean(value),
        freq = n()
      ) %>%
      arrange(desc(freq)) %>%
      filter(freq>min_presence & freq< max_presence)
  }
  
  #for FIs investigation
  if(method=="FIs"){
    #listing the most frequent fragments
    mz_list <- msms$mz
    all <- sort(unlist(mz_list))
    
    #grouping the m/z
    db <- dbscan(as.matrix(all), eps = 0.001, minPts = 3)
    
    df <-data.frame(
      value = all,
      cluster = db$cluster
    )
    
    df <- df %>%
      filter(cluster != 0) %>%    
      group_by(cluster) %>%
      summarise(
        val = mean(value),
        freq = n()
      ) %>%
      arrange(desc(freq)) %>%
      filter(freq>min_presence & freq< max_presence)
  }
  
  #creating the igraph object
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
  
  if(method=="NLs"){
    list <- msms$nl[match(V(graph)$name, msms$id)]}
  if(method=="FIs"){
    list <- msms$mz[match(V(graph)$name, msms$id)]}
  
  #df for results
  results <- data.frame(
    value = df$val, 
    number = NA_real_,
    z_score = NA_real_
  )
  
  #calculating unsuipervized community
  cl <- cluster_louvain(graph)
  V(graph)$unsupervised <- membership(cl)
  
  #calculating metrics fo each NL or mz
  for (i in seq_len(nrow(results))) {
    nli <- results$value[i]
    has_nl <- vapply(seq_along(list), function(j) {
      vals <- as.numeric(as.character(list[[j]]))
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
  results$log_num <- log10(results$number)
  qg <- qgam(
    z_score ~ s(log_num , k = K),
    data = results,
    qu = Q
  )
  
  q_pred <- predict(qg)
  results$excess <- results$z_score - q_pred
  results <- results %>%
    dplyr::select(-log_num)
  return(results)
}
