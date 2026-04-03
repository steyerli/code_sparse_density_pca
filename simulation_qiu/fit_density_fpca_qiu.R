library(densityFPCA)
# https://github.com/jiamingqiu/densityFPCA/blob/master/R/esti_fpca.R

fit_density_fpca_qiu <-function(x_data, sample_size_large = 30, 
                                x_grid = seq(min(unlist(x_data)), max(unlist(x_data)), length = 200),
                                bw = "nrd", adjust = 2, num.k = "AIC"){
  ###############
  # functional component analysis computed as in the R-package densityPCA
  ################
  # use only densities with many observations for pca computation
  is_large_sample <- sapply(x_data, function(x) length(x) >= sample_size_large)
  
  # kernel density estimates
  densities_estimated_large <- lapply(which(is_large_sample), function(i){
    density <- density(x_data[[i]], from = min(x_grid), to = max(x_grid), 
                       kernel = "gaussian", bw, adjust = adjust, n = length(x_grid))
    # move away from the boundary the of Bayes space
    density$y <- density$y + max(0, (0.001/(max(x_grid) - min(x_grid)) - min(density$y))) 
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
  #normalise phi assuming regular grid
  constant <- apply(fpca.res$phi, 2,  function(g){
    L_2_norm(cbind(x_grid, g))/(max(x_grid) - min(x_grid))
  })
  rotation <- t(t(fpca.res$phi)/constant)
  
  pca <- list(sdev = sqrt(fpca.res$lambda)*constant, 
              rotation = rotation,
              center = fpca.res$mu)
  
  # reconstruct the densities with small sample size 
  fpca.den.fam <- fpca2DenFam(fpca.res, control = list(num.k = 10))
  
  ls.fpca.esti <- fpcaEsti(
    mat.obsv = x_data,
    fpca.res = fpca.res,
    esti.method = c("FPCA_BLUP"),
    control = list(num.k = num.k,
      method = 'LBFGS', return.scale = 'origin'
    )
  )
  
  scores <- fpcaEsti(
    mat.obsv = x_data,
    fpca.res = fpca.res,
    esti.method = c("FPCA_BLUP"),
    control = list(num.k = num.k,
                   method = 'LBFGS', return.scale = 'parameter'
    )
  )$res*constant
  
  clr_densities_estimated <- lapply(1:nrow(ls.fpca.esti$res), function(i){
    data.frame("x" = ls.fpca.esti$grid, "y" = (max(x_grid) - min(x_grid))*ls.fpca.esti$res[i,])
  })
  
  densities_estimated <- lapply(clr_densities_estimated, inverse_clr_trafo)
  
  return(list("pca" = pca, "x_grid" = x_grid, "densities" = densities_estimated, 
              "clr_densities" = clr_densities_estimated, "scores" = scores))
}