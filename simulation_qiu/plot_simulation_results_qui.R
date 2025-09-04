library(reshape2)
library(ggplot2)
library(cowplot)
source("../fit_density_pca.R")

################################################################################
# plot simulation results for one example with m = 40 observations per density
################################################################################
# read simulation results
sim_results <- readRDS("simulation_results/simulation_40_1001")

type_levels <- c("oracle", "latent density",
                 "two-step, \nkernel density", 
                 "two-step, \ndensity fpca Qui")
################################################################################
### make covariance plots
get_cov_from_pca <- function(pca){
  pca$rotation%*%diag(pca$sdev^2)%*%t(pca$rotation)
}

cov_data_list <- lapply(1:4, function(i){
  data.frame(expand.grid(sim_results[[i]]$x_grid, sim_results[[i]]$x_grid), 
                        cov = melt(get_cov_from_pca(sim_results[[i]]$pca))[,3],
                        type = type_levels[i])
})

cov_plot_data <- do.call("rbind", cov_data_list)
cov_plot_data$type <- factor(cov_plot_data$type, levels = type_levels) 

g1 <- ggplot(data = cov_plot_data, aes(x = Var1, y = Var2, fill = cov)) + geom_raster() +
  coord_equal() + xlab("x") + ylab("x") +
  scale_fill_gradientn(colors = c("blue", "white", "red", "darkred")) +
  facet_grid(cols = vars(type)) + theme(strip.background = element_blank(), strip.text = element_blank())
g1
################################################################################
### make mean plus pc plot

mean_plot_data_list <- lapply(1:4, function(i){
  mean <- data.frame("x" = sim_results[[i]]$x_grid, y = sim_results[[i]]$pca$center, 
                     g = "mean", type =  type_levels[i])
  pc_1  <- data.frame("x" = sim_results[[i]]$x_grid, y = sim_results[[i]]$pca$rotation[,1])
  pc_1[,2] <- sign(pc_1[50,2])*pc_1[,2]/L_2_norm(pc_1)
  pc_1 <- cbind(pc_1, g = "PC 1", type = type_levels[i]) 
  pc_2  <- data.frame("x" = sim_results[[i]]$x_grid, y = sim_results[[i]]$pca$rotation[,2])
  pc_2[,2] <- sign(pc_2[100,2])*pc_2[,2]/L_2_norm(pc_2)
  pc_2 <- cbind(pc_2, g = "PC 2", type = type_levels[i])                    

  rbind(mean, pc_1, pc_2)
})

mean_plot_data <- do.call("rbind", mean_plot_data_list)
mean_plot_data$type <- factor(mean_plot_data$type, levels = type_levels) 

g2 <- ggplot(data = mean_plot_data, aes(x = x, y = y, group = g, col = g)) + geom_path(linewidth = 1) +
  facet_grid(cols = vars(type)) + scale_color_manual(values = c("black", "red", "blue")) +
  theme(strip.background = element_blank(), strip.text = element_blank()) + ylab("g(x)")
g2

################################################################################
### make predicted densities plot

density_plot_data_list <- lapply(1:4, function(j){
  n <- length(sim_results[[1]]$true_densities)
  density_data <- lapply(1:n, function(i){
    if(j == 1){
      density_i <- cbind(sim_results[[j]]$true_densities[[i]])  
    } else if (j == 2){
      prediction <- predict_latent_densities(sim_results[[j]])
      density_i <- cbind(inverse_clr_trafo(cbind("x" = prediction$x_grid, 
                                      "y" = prediction$clr_densities[,i])))
    } else {
      density_i <- cbind(sim_results[[j]]$densities[[i]])  
    }
    cbind(density_i, "idx" = i, "type" = type_levels[j], 
          "m" = length(sim_results[[1]]$x_data[[i]]))
  })
  do.call("rbind", density_data)
})

density_plot_data<- do.call("rbind", density_plot_data_list)
density_plot_data$type <- factor(density_plot_data$type, levels = type_levels) 
density_plot_data$idx <- factor(density_plot_data$idx)
density_plot_data$m <- factor(density_plot_data$m)

g3 <- ggplot(data = density_plot_data, aes(x = x, y = y, group = idx, color = m)) + geom_path() +
  facet_grid(cols = vars(type)) + ylab("density") + scale_color_manual(values = c("chocolate", "cadetblue")) +
  theme_bw()

g3
################################################################################
plot_grid(g3, g2, g1, ncol = 1, align = "v", axis = "lr")


################################################################################
# plot simulation performance for different m
################################################################################
m_vec <- c(20, 40, 80, 160)
dist_data_list <- lapply(m_vec, function(m){
  file_names <- list.files("simulation_results/", pattern = paste0("simulation_", m),
                           full.names = TRUE)
  dist_data_list <- lapply(file_names, function(file){
    simulation_result <- readRDS(file)
    cov_1 <- get_cov_from_pca(simulation_result[[1]]$pca)
    mean_1 <- simulation_result[[1]]$pca$center
    x_grid_1 <- simulation_result[[1]]$x_grid
    
    dist_data_list <- lapply(2:4, function(j){
      cov_j <- get_cov_from_pca(simulation_result[[j]]$pca)
      mean_j <- simulation_result[[j]]$pca$center
      x_grid_j <- simulation_result[[j]]$x_grid
      idx_j <- round(approx(y = seq_along(x_grid_1), x = x_grid_1, xout = x_grid_j)$y)
      dist_cov <- data.frame("x" = sqrt(sum((cov_1[idx_j, idx_j] - cov_j)^2)/
                               ((length(idx_j)*(max(x_grid_j) - min(x_grid_j)))^2)),
                             method = type_levels[j], fun = "cov")
      dist_mean <- data.frame("x" = sqrt(sum((mean_1[idx_j] -mean_j)^2)/
                                (length(idx_j)*(max(x_grid_j) - min(x_grid_j)))),
                              method = type_levels[j], fun = "mean")
      rbind(dist_cov, dist_mean)
    })
    dist_data <- do.call("rbind", dist_data_list)
  })
  dist_data <- do.call("rbind", dist_data_list)
  cbind(dist_data, m = m)
})

dist_data<- do.call("rbind", dist_data_list)
dist_data$m <- as.factor(dist_data$m)
dist_data$method <- factor(dist_data$method, levels = type_levels) 

g_mean <- ggplot(data = dist_data[dist_data$fun == "mean",], aes(x = m, y = x, fill = method)) + geom_boxplot() +
  ylab("dist mean") + xlab(expression(m[i])) + ggtitle("Distance of mean estimate to oracle estimate") +
  theme(legend.position = "top") + scale_fill_brewer(palette="BuPu")

g_cov <- ggplot(data = dist_data[dist_data$fun == "cov",], aes(x = m, y = x, fill = method)) + geom_boxplot() +
  ylab("dist covariance") + xlab(expression(m[i])) + ggtitle("Distance of covariance estimate to oracle estimate" ) +
  theme(legend.position = "top") + scale_fill_brewer(palette="BuPu")

plot_grid(g_mean, g_cov)




