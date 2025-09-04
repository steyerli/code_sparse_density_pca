library(densityFPCA)

fit_density_fpca_qui <-function(x_data, sample_size_large = 30, 
                                x_grid = seq(min(unlist(x_data)), max(unlist(x_data)), length = 200),
                                bw = (max(x_grid) - min(x_grid))/10){
  ###############
  # functional component analysis computed as in the R-package densityPCA
  ################
  # use only densities with many observations for pca computation
  is_large_sample <- sapply(x_data, function(x) length(x) >= sample_size_large)
  
  # kernel density estimates
  densities_estimated_large <- lapply(which(is_large_sample), function(i){
    density <- density(x_data[[i]], from = min(x_grid), to = max(x_grid), 
                       kernel = "gaussian", bw, 
                       n = length(x_grid))
    data.frame("x" = density$x, "y" = density$y)
  })
  # transform the density functions into Hilbert space via centered log transformation
  clr_densities_estimated_large <- lapply(densities_estimated_large, clr_trafo)
  
  # functional principal component analysis
  clr_densities_estimated_values <- lapply(clr_densities_estimated, function(dens) dens$y)
  
  fpca.res <- fdapace::FPCA(Ly = clr_densities_estimated_values,
                            Lt = replicate(x_grid, n = length(clr_densities_estimated), simplify = FALSE),
                            optns = list(
                              error = TRUE, lean = TRUE, FVEthreshold = 1,
                              methodSelectK = 'FVE', plot = FALSE, useBinnedData = 'OFF'
                            ))
  pca <- list(sdev = fpca.res$lambda, 
              rotation = fpca.res$phi,
              center = fpca.res$mu)
  
  # reconstruct the densities with small sample size 
  fpca.den.fam <- fpca2DenFam(fpca.res, control = list(num.k = 10))
  
  ls.fpca.esti <- fpcaEsti(
    mat.obsv = x_data[!is_large_sample],
    fpca.res = fpca.res,
    esti.method = c("FPCA_MLE"),
    control = list(num.k = "AIC",
      method = 'LBFGS', return.scale = 'origin'
    )
  )
  
  clr_densities_estimated_small <- lapply(1:nrow(ls.fpca.esti$res), function(i){
    data.frame("x" = ls.fpca.esti$grid, "y" = ls.fpca.esti$res[i,])
  })
  
  clr_densities_estimated <- c(clr_densities_estimated_small, clr_densities_estimated_large)
  densities_estimated <- lapply(clr_densities_estimated, inverse_clr_trafo)
  
  return(list("pca" = pca, "x_grid" = x_grid, "densities" = densities_estimated, "clr_densities" = clr_densities_estimated))
}