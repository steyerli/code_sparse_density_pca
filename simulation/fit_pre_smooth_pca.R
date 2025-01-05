fit_pre_smooth_pca <- function(x_data, x_grid = seq(min(unlist(x_data)), max(unlist(x_data)), length = 200),
                               bw = (max(x_grid) - min(x_grid))/10){
  # kernel density estimates
  densities_estimated <- lapply(1:length(x_data), function(i){
    density <- density(x_data[[i]], from = min(x_grid), to = max(x_grid), 
                       kernel = "gaussian", bw, 
                       n = length(x_grid))
    data.frame("x" = density$x, "y" = density$y)
  })
  # compute initial pca
  clr_densities_estimated <- lapply(densities_estimated, clr_trafo)
  clr_densities <- do.call("rbind", sapply(clr_densities_estimated, '[', 2))
  pca <- prcomp(na.omit(clr_densities))
  return(list("pca" = pca, "x_grid" = x_grid, "densities" = densities_estimated, "clr_densities" = clr_densities_estimated))
}