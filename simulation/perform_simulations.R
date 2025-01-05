source("simulate_data.R")
source("../fit_density_pca.R")
source("fit_compositional_spline_pca.R")
source("fit_pre_smooth_pca.R")
################################################################################
m <- 20
lapply(1:100, function(i){
  seed_i <- 1000 + i
  set.seed(seed_i)
  simulated_data <- simulate_data(m = m)
  pre_smooth_pca <- fit_pre_smooth_pca(simulated_data$x_data, bw = 0.12)
  compo_pca <- fit_compositional_spline_pca(simulated_data$x_data)
  density_pca <- fit_density_pca(simulated_data$x_data, max_iter = 50, dim_reduction = 0.0001, bw = 0.12)
  saveRDS(list(simulated_data, density_pca, pre_smooth_pca, compo_pca), 
          file = paste0("simulation_results/simulation_", m,"_", seed_i))
})

################################################################################
m <- 40
lapply(1:100, function(i){
  seed_i <- 1000 + i
  set.seed(seed_i)
  simulated_data <- simulate_data(m = m)
  pre_smooth_pca <- fit_pre_smooth_pca(simulated_data$x_data)
  compo_pca <- fit_compositional_spline_pca(simulated_data$x_data)
  density_pca <- fit_density_pca(simulated_data$x_data, max_iter = 50)
  saveRDS(list(simulated_data, density_pca, pre_smooth_pca, compo_pca), 
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




