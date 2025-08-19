library(densityFPCA)

fit_density_fpca_qui <-function(x_data, x_grid = seq(min(unlist(x_data)), max(unlist(x_data)), length = 200),
                                bw = (max(x_grid) - min(x_grid))/10){
  ###############
  # functional component analysis computed as in the R-package densityPCA
  ################
  # kernel density estimates
  densities_estimated <- lapply(1:length(x_data), function(i){
    density <- density(x_data[[i]], from = min(x_grid), to = max(x_grid), 
                       kernel = "gaussian", bw, 
                       n = length(x_grid))
    data.frame("x" = density$x, "y" = density$y)
  })
  # transform the density functions into Hilbert space via centered log transformation
  clr_densities_estimated <- lapply(densities_estimated, clr_trafo)
  clr_densities_estimated_values <- lapply(clr_densities_estimated, function(dens) dens$y)
  
  # functional principal component analysis
  fpca.res <- fdapace::FPCA(Ly = clr_densities_estimated_values,
                            Lt = replicate(x_grid, n = length(clr_densities_estimated), simplify = FALSE),
                            optns = list(
                              error = TRUE, lean = TRUE, FVEthreshold = 1,
                              methodSelectK = 'FVE', plot = FALSE, useBinnedData = 'OFF'
                            ))
  pca <- list(sdev = fpca.res$lambda, 
              rotation = fpca.res$phi,
              center = fpca.res$mu)
  
  return(list("pca" = pca, "x_grid" = x_grid, "densities" = densities_estimated, "clr_densities" = clr_densities_estimated))
  
  
  #fpca.den.fam <- fpca2DenFam(fpca.res, control = list(num.k = 10))
  
}