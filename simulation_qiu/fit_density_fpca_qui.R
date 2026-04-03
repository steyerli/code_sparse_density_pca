library(densityFPCA)

fit_density_fpca_qui <-function(x_data, sample_size_large = 30, 
                                x_grid = seq(min(unlist(x_data)), max(unlist(x_data)), length = 200),
                                bw = NULL){
  ###############
  # functional component analysis computed as in the R-package densityPCA
  ################
  # use only densities with many observations for pca computation
  is_large_sample <- sapply(x_data, function(x) length(x) >= sample_size_large)
  
  # kernel density estimates
  densities_estimated_large <- lapply(which(is_large_sample), function(i){
    if (is.null(bw)){
      bw_i <- max(diff(c(min(x_grid), sort(x_data[[i]]), max(x_grid))))/2
    } else {
      bw_i <- bw
    }
    density <- density(x_data[[i]], from = min(x_grid), to = max(x_grid), 
                       kernel = "gaussian", bw_i, 
                       n = length(x_grid))
    data.frame("x" = density$x, "y" = density$y)
  })
  # transform the density functions into Hilbert space via centered log transformation
  clr_densities_estimated_large <- lapply(densities_estimated_large, clr_trafo)
  
  # functional principal component analysis
  clr_densities_estimated_values <- lapply(clr_densities_estimated_large, function(dens) dens$y)
  
  fpca.res <- fdapace::FPCA(Ly = clr_densities_estimated_values,
                            Lt = replicate(x_grid, n = length(clr_densities_estimated_large), simplify = FALSE),
                            optns = list(
                              error = FALSE, lean = TRUE, FVEthreshold = 1,
                              usergrid = TRUE, 
                              methodSelectK = 'FVE', plot = FALSE, useBinnedData = 'OFF'
                            ))
  #normalise phi
  rotation <- apply(fpca.res$phi, 2, function(x) x / sqrt(sum(x^2)))
  
  pca <- list(sdev = sqrt(fpca.res$lambda/(x_grid[2] - x_grid[1])), 
              rotation = rotation,
              center = fpca.res$mu)
  
  # reconstruct the densities with small sample size 
  fpca.den.fam <- fpca2DenFam(fpca.res, control = list(num.k = 10))
  
  ls.fpca.esti <- fpcaEsti(
    mat.obsv = x_data,
    fpca.res = fpca.res,
    esti.method = c("FPCA_BLUP"),
    control = list(num.k = "AIC",
      method = 'LBFGS', return.scale = 'origin'
    )
  )
  
  clr_densities_estimated <- lapply(1:nrow(ls.fpca.esti$res), function(i){
    data.frame("x" = ls.fpca.esti$grid, "y" = ls.fpca.esti$res[i,])
  })
  
  densities_estimated <- lapply(clr_densities_estimated, inverse_clr_trafo)
  
  return(list("pca" = pca, "x_grid" = x_grid, "densities" = densities_estimated, "clr_densities" = clr_densities_estimated))
}