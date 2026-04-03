source("../fit_density_pca.R")

simulate_data_misspec <- function(n = 30, m = 20){
  x_grid <- seq(0.1,0.9,0.005)
  
  
  #######################################
  # true densities
  true_observed_densities <- lapply(1:n, function(i){
    shape1 <- rbeta(1, 0.5, 0.5) + 1
    shape2 <- runif(1, 0.5, 5)
    data.frame("x" = x_grid, "y" = dbeta(x_grid, shape1, shape2))
  })
  true_observed_clr_densities <- sapply(1:n, function(i){
    clr_trafo(true_observed_densities[[i]])[,2]
  })
  
  ######################################
  x_data <- lapply(1:n, function(i){
    probs <- true_observed_densities[[i]][,2]
    x_grid <- true_observed_densities[[i]][,1]
    sample(x_grid, m, replace = TRUE, prob = probs)
  })
  ######################################
  oracle_pca <- prcomp(data.frame(t(true_observed_clr_densities)))
  #compute normalizing constants
  constant <- apply(oracle_pca$rotation, 2,  function(g){
    L_2_norm(cbind(x_grid, g))
  })
  oracle_pca$rotation <- t(t(oracle_pca$rotation)/constant)
  oracle_pca$sdev <- oracle_pca$sdev*constant
  
  return(list("x_data" = x_data, "x_grid" = x_grid, "true_densities" = true_observed_densities, "pca" = oracle_pca))
}

