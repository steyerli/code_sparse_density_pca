fit_pre_smooth_pca <- function(x_data, x_grid = seq(min(unlist(x_data)), max(unlist(x_data)), length = 200),
                               bw = "nrd", adjust = 2){
  # kernel density estimates
  densities_estimated <- lapply(1:length(x_data), function(i){
    density <- try(density(x_data[[i]], from = min(x_grid), to = max(x_grid), 
                       kernel = "gaussian", bw, adjust = adjust, n = length(x_grid)),
                   silent = TRUE)
    if(inherits(density, "try-error")){
      bw <- (max(x_grid) -min(x_grid))/10
      density <- try(density(x_data[[i]], from = min(x_grid), to = max(x_grid), 
                             kernel = "gaussian", bw, adjust = adjust, n = length(x_grid)))
    }
    # move away from the boundary the of Bayes space
    density$y <- density$y + max(0, (0.001/(max(x_grid) - min(x_grid)) - min(density$y))) 
    data.frame("x" = density$x, "y" = density$y)
  })
  # compute initial pca
  clr_densities_estimated <- lapply(densities_estimated, clr_trafo)
  densities_estimated <- lapply(clr_densities_estimated, inverse_clr_trafo)
  clr_densities <- do.call("rbind", sapply(clr_densities_estimated, '[', 2))
  pca <- prcomp(na.omit(clr_densities))
  return(list("pca" = pca, "x_grid" = x_grid, "densities" = densities_estimated, "clr_densities" = clr_densities_estimated))
}