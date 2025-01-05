source("../fit_density_pca.R")

simulate_data <- function(n = 30, m = 20, var_1 = 0.5, var_2 = 0.2){
  x_grid <- seq(0,1,0.005)
  mu <- data.frame("x" = x_grid, "y" = -20*(x_grid - 0.5)^2)
  clr_mean <- center_function(mu)
  
  pc_1 <- data.frame("x" = x_grid, "y" = 0.2*sin(10*(x_grid-0.5)))
  pc_1 <- center_function(pc_1)
  pc_1[,2] <- pc_1[,2]/L_2_norm(pc_1)
  #####
  pc_2 <- data.frame("x" = x_grid, "y" = 0.1*cos(2*pi*(x_grid -0.5)))
  pc_2 <- center_function(pc_2)
  pc_2[,2] <- pc_2[,2]/L_2_norm(pc_2)
  
  #######################################
  # true densities
  true_observed_clr_densities <- sapply(1:n, function(i){
    clr_mean[,2] + rnorm(1, 0, sqrt(var_1))*pc_1[,2] + rnorm(1, 0, sqrt(var_2))*pc_2[,2]
  })
  true_observed_densities <- lapply(1:n, function(i){
    clr_density <- data.frame(x_grid, true_observed_clr_densities[,i])
    inverse_clr_trafo(clr_density)
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

