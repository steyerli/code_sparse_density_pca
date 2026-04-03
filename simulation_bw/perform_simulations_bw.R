# set working directory to source file directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../simulation/simulate_data.R")
source("../fit_density_pca.R")
source("../simulation/fit_pre_smooth_pca.R")

################################################################################
m <- 40
adjust <- 1

lapply(1:100, function(i){
  seed_i <- 1000 + i
  set.seed(seed_i)
  simulated_data <- simulate_data(m = m)
  pre_smooth_pca <- fit_pre_smooth_pca(simulated_data$x_data, adjust = adjust)
  density_pca <- fit_density_pca(simulated_data$x_data, max_iter = 50, adjust = adjust)
  saveRDS(list(simulated_data, density_pca, pre_smooth_pca), 
          file = paste0("simulation_results/simulation_", 10*adjust,"_", seed_i))
})

#######################
m <- 40
adjust <- 1.5

lapply(1:100, function(i){
  seed_i <- 1000 + i
  set.seed(seed_i)
  simulated_data <- simulate_data(m = m)
  pre_smooth_pca <- fit_pre_smooth_pca(simulated_data$x_data, adjust = adjust)
  density_pca <- fit_density_pca(simulated_data$x_data, max_iter = 50, adjust = adjust)
  saveRDS(list(simulated_data, density_pca, pre_smooth_pca), 
          file = paste0("simulation_results/simulation_", 10*adjust,"_", seed_i))
})

#######################
m <- 40
adjust <- 2

lapply(1:100, function(i){
  seed_i <- 1000 + i
  set.seed(seed_i)
  simulated_data <- simulate_data(m = m)
  pre_smooth_pca <- fit_pre_smooth_pca(simulated_data$x_data, adjust = adjust)
  density_pca <- fit_density_pca(simulated_data$x_data, max_iter = 50, adjust = adjust)
  saveRDS(list(simulated_data, density_pca, pre_smooth_pca), 
          file = paste0("simulation_results/simulation_", 10*adjust,"_", seed_i))
})

#######################
m <- 40
adjust <- 2.5

lapply(1:100, function(i){
  seed_i <- 1000 + i
  set.seed(seed_i)
  simulated_data <- simulate_data(m = m)
  pre_smooth_pca <- fit_pre_smooth_pca(simulated_data$x_data, adjust = adjust)
  density_pca <- fit_density_pca(simulated_data$x_data, max_iter = 50, adjust = adjust)
  saveRDS(list(simulated_data, density_pca, pre_smooth_pca), 
          file = paste0("simulation_results/simulation_", 10*adjust,"_", seed_i))
})

#######################

m <- 40
adjust <- 3

lapply(1:100, function(i){
  seed_i <- 1000 + i
  set.seed(seed_i)
  simulated_data <- simulate_data(m = m)
  pre_smooth_pca <- fit_pre_smooth_pca(simulated_data$x_data, adjust = adjust)
  density_pca <- fit_density_pca(simulated_data$x_data, max_iter = 50, adjust = adjust)
  saveRDS(list(simulated_data, density_pca, pre_smooth_pca), 
          file = paste0("simulation_results/simulation_", 10*adjust,"_", seed_i))
})
