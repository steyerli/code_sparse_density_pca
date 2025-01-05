library(splines)

fit_compositional_spline_pca <- function(x_data, n_knots = 5, ord = 4){
  knots <- seq(min(unlist(x_data)), max(unlist(x_data)), length = n_knots)
  x_grid = seq(min(unlist(x_data)), max(unlist(x_data)), length = 200)
  
  clr_densities_estimated <- lapply(1:length(x_data), function(i){
    hist <- hist(x_data[[i]], plot = FALSE, breaks = seq(min(knots), max(knots), length = n_knots + ord +2))
    midpoints <- hist$breaks[-1] - 0.5*diff(hist$breaks)
    hist$counts[hist$counts == 0] <- 0.5
    hist_data <- cbind(midpoints, log(hist$counts) - mean(log(hist$counts)))
    
    # fit spline
    knots_long <- knots[c(rep(1,ord -1), 1:length(knots), rep(length(knots), ord - 1))]
    bb <- splineDesign(knots_long, x = hist_data[,1], outer.ok = FALSE, ord = ord, derivs = 1)
    bb <- bb[,2:(ncol(bb) - 1)]
    spline_coefs <- lm(hist_data[,2] ~ bb -1)$coef
    bb_new <- splineDesign(knots_long, x = x_grid, ord = ord, derivs = 1)
    bb_new <- bb_new[,2:(ncol(bb_new) - 1)]
    
    data.frame("x" = x_grid, "y" = bb_new%*%spline_coefs)
  })
  densities_estimated <- lapply(clr_densities_estimated, inverse_clr_trafo)
  
  clr_densities <- do.call("rbind", sapply(clr_densities_estimated, '[', 2))
  pca <- prcomp(na.omit(clr_densities))
  return(list("pca" = pca, "x_grid" = x_grid, "densities" = densities_estimated, "clr_densities" = clr_densities_estimated))
} 