source("../simulation/simulate_data.R")
source("../fit_density_pca.R")
source("../simulation/fit_pre_smooth_pca.R")
source("fit_density_fpca_qui.R")

################################################################################
# m <- 20
# lapply(1:100, function(i){
#   seed_i <- 1000 + i
#   set.seed(seed_i)
#   simulated_data <- simulate_data(n = 40, m = m)
#   pre_smooth_pca <- fit_pre_smooth_pca(simulated_data$x_data, bw = 0.12)
#   compo_pca <- fit_compositional_spline_pca(simulated_data$x_data)
#   density_pca <- fit_density_pca(simulated_data$x_data, max_iter = 50, dim_reduction = 0.0001, bw = 0.12)
#   saveRDS(list(simulated_data, density_pca, pre_smooth_pca, compo_pca), 
#           file = paste0("simulation_results/simulation_", m,"_", seed_i))
# })

################################################################################
m <- 40
lapply(1:100, function(i){
  seed_i <- 1000 + i
  set.seed(seed_i)
  simulated_data <- simulate_data(n = 40, m = m)
  # keep for first half of densities only 20 observations
  simulated_data$x_data[1:20] <- lapply(simulated_data$x_data[1:20], "[", 1:20)
  # fit density pca using all three methods
  pre_smooth_pca <- fit_pre_smooth_pca(simulated_data$x_data)
  density_fpca_qui <- fit_density_fpca_qui(simulated_data$x_data)
  density_pca <- fit_density_pca(simulated_data$x_data, max_iter = 50)
  
  saveRDS(list(simulated_data, density_pca, pre_smooth_pca, density_fpca_qui), 
          file = paste0("simulation_results/simulation_", m,"_", seed_i))
})

#######################
m <- 80
lapply(1:100, function(i){
  seed_i <- 1000 + i
  set.seed(seed_i)
  simulated_data <- simulate_data(m = m)
  pre_smooth_pca <- fit_pre_smooth_pca(simulated_data$x_data, bw = 0.08)
  compo_pca <- fit_compositional_spline_pca(simulated_data$x_data)
  density_pca <- fit_density_pca(simulated_data$x_data, max_iter = 50, bw = 0.08)
  saveRDS(list(simulated_data, density_pca, pre_smooth_pca, compo_pca), 
          file = paste0("simulation_results/simulation_", m,"_", seed_i))
})

#######################
m <- 160
lapply(1:100, function(i){
  seed_i <- 1000 + i
  set.seed(seed_i)
  simulated_data <- simulate_data(m = m)
  pre_smooth_pca <- fit_pre_smooth_pca(simulated_data$x_data, bw = 0.07)
  compo_pca <- fit_compositional_spline_pca(simulated_data$x_data)
  density_pca <- fit_density_pca(simulated_data$x_data, max_iter = 50, bw = 0.07)
  saveRDS(list(simulated_data, density_pca, pre_smooth_pca, compo_pca), 
          file = paste0("simulation_results/simulation_", m,"_", seed_i))
})




