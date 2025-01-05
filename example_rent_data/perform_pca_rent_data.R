# working directory is assumed to be the location of this file
source("../fit_density_pca.R")
library(ggplot2)
library(cowplot)
library(foreign)
library(manipulate)
library(sf)
library(geojsonio)
library(xtable)
library(stringr)

#########
### read in data ###
load("msp_paper.RData")
rents <-  msp.paper

rents$nmqm <- as.numeric(rents$nmqm)
rents$sb_nummer <- factor(rents$kfbezirk)

################################################################################
munich_sf <- geojson_sf("munich.json")
munich_sb <- unique(as.data.frame(munich_sf)[c("sb_nummer", "name")])
munich_sb$sb_nummer <- as.numeric(munich_sb$sb_nummer)
munich_sb <- munich_sb[order(munich_sb$sb_nummer), ]
row.names(munich_sb) <- NULL
munich_sf <- munich_sf[c("sb_nummer")]
munich_sf$sb_nummer <- as.numeric(munich_sf$sb_nummer)
munich_sb$m <- table(rents$sb_nummer)
#########
#print table
table <- xtable(munich_sb)
digits(table) <- 0
colnames(table) <- c("district $i$", "name", "$m_i$")
print(table, include.rownames = FALSE, sanitize.text.function=function(x){x})

################################################################################
# pre-smooth densities
# kernel density estimates
x_grid <- seq(min(rents$nmqm), max(rents$nmqm), length = 200)
densities_estimated <- lapply(1:25, function(i){
  density <- density(rents$nmqm[rents$sb_nummer == i], 
                     from = min(x_grid), to = max(x_grid), 
                     kernel = "gaussian", bw = 2, 
                     n = length(x_grid))
  data.frame("x" = density$x, "y" = density$y, "district" = i)
})
densities_estimated_plot_data <- do.call("rbind", densities_estimated)
densities_estimated_plot_data$district <- factor(densities_estimated_plot_data$district)
ggplot() + geom_path(aes(x = x, y = y, group = district, col = district), data = densities_estimated_plot_data)

# pca on pre-smoothed data
clr_densities_estimated <- lapply(densities_estimated, clr_trafo)
clr_densities_data <- do.call("rbind", clr_densities_estimated)
clr_densities_data$district <- factor(clr_densities_data$district)
ggplot() + geom_path(aes(x = x, y = y, group = district, col = district), data = clr_densities_data)

pre_smoothed_pca <- prcomp(na.omit(t(sapply(clr_densities_estimated, "[[", "y"))))
constant <- apply(pre_smoothed_pca$rotation, 2,  function(g){
  L_2_norm(cbind(x_grid, g))/(max(x_grid) - min(x_grid))
})
pre_smoothed_pca$sdev <- pre_smoothed_pca$sdev*constant
pre_smoothed_pca$sdev[1:10]
# choose orientation of principal components
constant <- constant*c(-1,1,-1, rep(1, length(pre_smoothed_pca$sdev) - 3))
pre_smoothed_pca$rotation <- t(t(pre_smoothed_pca$rotation)/constant)
pre_smoothed_pca$x <- t(t(pre_smoothed_pca$x)*constant)

plot_pca(pre_smoothed_pca, x_grid, dim = 3)
################################################################################
# plot pre-smoothed pca
munich_sf$score1_pre_smooth <- pre_smoothed_pca$x[munich_sf$sb_nummer, 1]
munich_sf$score2_pre_smooth <- pre_smoothed_pca$x[munich_sf$sb_nummer, 2]
munich_sf$score3_pre_smooth <- pre_smoothed_pca$x[munich_sf$sb_nummer, 3]
sdvs <- pre_smoothed_pca$sdev

i <- 1
pc1_plot_data <- rbind(cbind(inverse_clr_trafo(cbind(x_grid, pre_smoothed_pca$center)), "score" = 0, "clr_trafo" = FALSE),
                       cbind(inverse_clr_trafo(cbind( x_grid, pre_smoothed_pca$center - sdvs[i]*pre_smoothed_pca$rotation[, i])), 
                             "score" = -round(sdvs[i], 2), "clr_trafo" = FALSE),
                       cbind(inverse_clr_trafo(cbind(x_grid, pre_smoothed_pca$center + sdvs[i]*pre_smoothed_pca$rotation[, i])), 
                             "score" = round(sdvs[i], 2), "clr_trafo" = FALSE),
                       cbind("x" = x_grid, "y" = pre_smoothed_pca$center, "score" = 0, "clr_trafo" = TRUE),
                       cbind("x" = x_grid, "y" = pre_smoothed_pca$center - sdvs[i]*pre_smoothed_pca$rotation[, i], 
                             "score" = -round(sdvs[i], 2), "clr_trafo" = TRUE),
                       cbind("x" = x_grid, "y" = pre_smoothed_pca$center + sdvs[i]*pre_smoothed_pca$rotation[, i], 
                             "score" = round(sdvs[i], 2), "clr_trafo" = TRUE))
pc1_plot_data$score <- factor(pc1_plot_data$score)
i<- 2
pc2_plot_data <- rbind(cbind(inverse_clr_trafo(cbind(x_grid, pre_smoothed_pca$center)), "score" = 0, "clr_trafo" = FALSE),
                       cbind(inverse_clr_trafo(cbind( x_grid, pre_smoothed_pca$center - sdvs[i]*pre_smoothed_pca$rotation[, i])), 
                             "score" = -round(sdvs[i], 2), "clr_trafo" = FALSE),
                       cbind(inverse_clr_trafo(cbind(x_grid, pre_smoothed_pca$center + sdvs[i]*pre_smoothed_pca$rotation[, i])), 
                             "score" = round(sdvs[i], 2), "clr_trafo" = FALSE),
                       cbind("x" = x_grid, "y" = pre_smoothed_pca$center, "score" = 0, "clr_trafo" = TRUE),
                       cbind("x" = x_grid, "y" = pre_smoothed_pca$center - sdvs[i]*pre_smoothed_pca$rotation[, i], 
                             "score" = -round(sdvs[i], 2), "clr_trafo" = TRUE),
                       cbind("x" = x_grid, "y" = pre_smoothed_pca$center + sdvs[i]*pre_smoothed_pca$rotation[, i], 
                             "score" = round(sdvs[i], 2), "clr_trafo" = TRUE))
pc2_plot_data$score <- factor(pc2_plot_data$score)
i <- 3
pc3_plot_data <- rbind(cbind(inverse_clr_trafo(cbind(x_grid, pre_smoothed_pca$center)), "score" = 0, "clr_trafo" = FALSE),
                       cbind(inverse_clr_trafo(cbind( x_grid, pre_smoothed_pca$center - sdvs[i]*pre_smoothed_pca$rotation[, i])), 
                             "score" = -round(sdvs[i], 2), "clr_trafo" = FALSE),
                       cbind(inverse_clr_trafo(cbind(x_grid, pre_smoothed_pca$center + sdvs[i]*pre_smoothed_pca$rotation[, i])), 
                             "score" = round(sdvs[i], 2), "clr_trafo" = FALSE),
                       cbind("x" = x_grid, "y" = pre_smoothed_pca$center, "score" = 0, "clr_trafo" = TRUE),
                       cbind("x" = x_grid, "y" = pre_smoothed_pca$center - sdvs[i]*pre_smoothed_pca$rotation[, i], 
                             "score" = -round(sdvs[i], 2), "clr_trafo" = TRUE),
                       cbind("x" = x_grid, "y" = pre_smoothed_pca$center + sdvs[i]*pre_smoothed_pca$rotation[, i], 
                             "score" = round(sdvs[i], 2), "clr_trafo" = TRUE))
pc3_plot_data$score <- factor(pc3_plot_data$score)

################################################################################
g_11 <- ggplot() + geom_sf(data = munich_sf, aes(fill = score1_pre_smooth), alpha =0.7) + theme_void() +
  guides(fill=guide_colorbar(title= "score")) + ylab("") + xlab("") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  geom_sf_label(data = munich_sf[-c(3,17),], aes(label = sb_nummer), label.size  = NA, alpha = 0, size = 3) +
  theme(legend.key.height=unit(2,"pt"), legend.key.width = unit(30, "pt"), legend.position = "bottom") 
  

g_12 <- ggplot(data = pc1_plot_data[pc1_plot_data$clr_trafo == 1,], aes(x = x, y = y, group = score, color = score)) + 
  geom_path(linewidth = 0.5) + theme(legend.position = "top", legend.key.height=unit(4,"pt")) +
  guides(color= guide_legend(title= "score")) + 
  scale_color_manual(values = c("blue", "black", "red")) + ylab("clr density") + xlab("rent per sqm") + 
  ggtitle(paste0("PC 1, ", round((sdvs/sum(sdvs))[1]*100), "% of variability"))

g_13 <- ggplot(data = pc1_plot_data[pc1_plot_data$clr_trafo == 0,], aes(x = x, y = y, group = score, color = score)) + 
  geom_path(linewidth = 0.5) + theme(legend.position = "none") +
  scale_color_manual(values = c("blue", "black", "red")) + ylab("density") + xlab("rent per sqm") 
################################################################################
g_21 <- ggplot() + geom_sf(data = munich_sf, aes(fill = score2_pre_smooth), alpha =0.7) + theme_void() +
  guides(fill=guide_colorbar(title= "score")) + ylab("") + xlab("") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  geom_sf_label(data = munich_sf[-c(3,17),], aes(label = sb_nummer), label.size  = NA, alpha = 0, size = 3) +
  theme(legend.key.height=unit(2,"pt"), legend.key.width = unit(30, "pt"), legend.position = "bottom") 


g_22 <- ggplot(data = pc2_plot_data[pc2_plot_data$clr_trafo == 1,], aes(x = x, y = y, group = score, color = score)) + 
  geom_path(linewidth = 0.5) + theme(legend.position = "top", legend.key.height=unit(4,"pt")) +
  guides(color= guide_legend(title= "score")) + 
  scale_color_manual(values = c("blue", "black", "red")) + ylab("clr density") + xlab("rent per sqm") + 
  ggtitle(paste0("PC 2, ", round((sdvs/sum(sdvs))[2]*100), "% of variability"))

g_23 <- ggplot(data = pc2_plot_data[pc2_plot_data$clr_trafo == 0,], aes(x = x, y = y, group = score, color = score)) + 
  geom_path(linewidth = 0.5) + theme(legend.position = "none") +
  scale_color_manual(values = c("blue", "black", "red")) + ylab("density") + xlab("rent per sqm") 
################################################################################
g_31 <- ggplot() + geom_sf(data = munich_sf, aes(fill = score3_pre_smooth), alpha =0.7) + theme_void() +
  guides(fill=guide_colorbar(title= "score")) + ylab("") + xlab("") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  geom_sf_label(data = munich_sf[-c(3,17),], aes(label = sb_nummer), label.size  = NA, alpha = 0, size = 3) +
  theme(legend.key.height=unit(2,"pt"), legend.key.width = unit(30, "pt"), legend.position = "bottom") 


g_32 <- ggplot(data = pc3_plot_data[pc3_plot_data$clr_trafo == 1,], aes(x = x, y = y, group = score, color = score)) + 
  geom_path(linewidth = 0.5) + theme(legend.position = "top", legend.key.height=unit(4,"pt")) +
  guides(color= guide_legend(title= "score")) + 
  scale_color_manual(values = c("blue", "black", "red")) + ylab("clr density") + xlab("rent per sqm") + 
  ggtitle(paste0("PC 3, ", round((sdvs/sum(sdvs))[3]*100), "% of variability"))

g_33 <- ggplot(data = pc3_plot_data[pc3_plot_data$clr_trafo == 0,], aes(x = x, y = y, group = score, color = score)) + 
  geom_path(linewidth = 0.5) + theme(legend.position = "none") +
  scale_color_manual(values = c("blue", "black", "red")) + ylab("density") + xlab("rent per sqm") 
################################################################################
plot_grid(g_12, g_22, g_32, g_13, g_23, g_33, g_11, g_21, g_31, 
          nrow = 3, align = "v", axis = "lr", rel_heights = c(2.5,2,3))

################################################################################
################################################################################
# fit density pca, DO NOT RUN!
set.seed(123)
x_data <- lapply(1:25, function(i){rents$nmqm[rents$sb_nummer == i]})
density_pca <- fit_density_pca(x_data = x_data, dim_reduction = 0.001, r = 20, lambda = 1,
                               max_iter = 50, bw = 2)
saveRDS(density_pca, file = "density_pca_mietspiegel.RDS")
################################################################################
################################################################################
density_pca <- readRDS("density_pca_mietspiegel.RDS")
plot_pca(density_pca$pca, x_grid = density_pca$x_grid, dim = 3)
scores <- data.frame(t(predict_latent_densities(density_pca)$predicted_scores), "sb_nummer" = 1:25)

# choose orientation of principal components
constant <- c(-1,1, -1, rep(1, length(density_pca$pca$sdev) - 3))
density_pca$pca$rotation <- t(t(density_pca$pca$rotation)/constant)
scores[, 1:length(density_pca$pca$sdev)] <- t(t(scores[, 1:length(density_pca$pca$sdev)])*constant)
################################################################################

munich_sf$score1 <- scores[,1][munich_sf$sb_nummer]
munich_sf$score2 <- scores[,2][munich_sf$sb_nummer]
munich_sf$score3 <- scores[,3][munich_sf$sb_nummer]
sdvs <- density_pca$pca$sdev

i <- 1
pc1_plot_data <- rbind(cbind(inverse_clr_trafo(cbind(x_grid, density_pca$pca$center)), "score" = 0, "clr_trafo" = FALSE),
                       cbind(inverse_clr_trafo(cbind( x_grid, density_pca$pca$center - sdvs[i]*density_pca$pca$rotation[, i])), 
                             "score" = -round(sdvs[i], 2), "clr_trafo" = FALSE),
                       cbind(inverse_clr_trafo(cbind(x_grid, density_pca$pca$center + sdvs[i]*density_pca$pca$rotation[, i])), 
                             "score" = round(sdvs[i], 2), "clr_trafo" = FALSE),
                       cbind("x" = x_grid, "y" = density_pca$pca$center, "score" = 0, "clr_trafo" = TRUE),
                       cbind("x" = x_grid, "y" = density_pca$pca$center - sdvs[i]*density_pca$pca$rotation[, i], 
                             "score" = -round(sdvs[i], 2), "clr_trafo" = TRUE),
                       cbind("x" = x_grid, "y" = density_pca$pca$center + sdvs[i]*density_pca$pca$rotation[, i], 
                             "score" = round(sdvs[i], 2), "clr_trafo" = TRUE))
pc1_plot_data$score <- factor(pc1_plot_data$score)
i<- 2
pc2_plot_data <- rbind(cbind(inverse_clr_trafo(cbind(x_grid, density_pca$pca$center)), "score" = 0, "clr_trafo" = FALSE),
                       cbind(inverse_clr_trafo(cbind( x_grid, density_pca$pca$center - sdvs[i]*density_pca$pca$rotation[, i])), 
                             "score" = -round(sdvs[i], 2), "clr_trafo" = FALSE),
                       cbind(inverse_clr_trafo(cbind(x_grid, density_pca$pca$center + sdvs[i]*density_pca$pca$rotation[, i])), 
                             "score" = round(sdvs[i], 2), "clr_trafo" = FALSE),
                       cbind("x" = x_grid, "y" = density_pca$pca$center, "score" = 0, "clr_trafo" = TRUE),
                       cbind("x" = x_grid, "y" = density_pca$pca$center - sdvs[i]*density_pca$pca$rotation[, i], 
                             "score" = -round(sdvs[i], 2), "clr_trafo" = TRUE),
                       cbind("x" = x_grid, "y" = density_pca$pca$center + sdvs[i]*density_pca$pca$rotation[, i], 
                             "score" = round(sdvs[i], 2), "clr_trafo" = TRUE))
pc2_plot_data$score <- factor(pc2_plot_data$score)
i <- 3
pc3_plot_data <- rbind(cbind(inverse_clr_trafo(cbind(x_grid, density_pca$pca$center)), "score" = 0, "clr_trafo" = FALSE),
                       cbind(inverse_clr_trafo(cbind( x_grid, density_pca$pca$center - sdvs[i]*density_pca$pca$rotation[, i])), 
                             "score" = -round(sdvs[i], 2), "clr_trafo" = FALSE),
                       cbind(inverse_clr_trafo(cbind(x_grid, density_pca$pca$center + sdvs[i]*density_pca$pca$rotation[, i])), 
                             "score" = round(sdvs[i], 2), "clr_trafo" = FALSE),
                       cbind("x" = x_grid, "y" = density_pca$pca$center, "score" = 0, "clr_trafo" = TRUE),
                       cbind("x" = x_grid, "y" = density_pca$pca$center - sdvs[i]*density_pca$pca$rotation[, i], 
                             "score" = -round(sdvs[i], 2), "clr_trafo" = TRUE),
                       cbind("x" = x_grid, "y" = density_pca$pca$center + sdvs[i]*density_pca$pca$rotation[, i], 
                             "score" = round(sdvs[i], 2), "clr_trafo" = TRUE))
pc3_plot_data$score <- factor(pc3_plot_data$score)
################################################################################
g_11 <- ggplot() + geom_sf(data = munich_sf, aes(fill = score1), alpha =0.7) + theme_void() +
  guides(fill=guide_colorbar(title= "score")) + ylab("") + xlab("") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  geom_sf_label(data = munich_sf[-c(3,17),], aes(label = sb_nummer), label.size  = NA, alpha = 0, size = 3) +
  theme(legend.key.height=unit(2,"pt"), legend.key.width = unit(30, "pt"), legend.position = "bottom") 


g_12 <- ggplot(data = pc1_plot_data[pc1_plot_data$clr_trafo == 1,], aes(x = x, y = y, group = score, color = score)) + 
  geom_path(linewidth = 0.5) + theme(legend.position = "top", legend.key.height=unit(4,"pt")) +
  guides(color= guide_legend(title= "score")) + 
  scale_color_manual(values = c("blue", "black", "red")) + ylab("clr density") + xlab("rent per sqm") + 
  ggtitle(paste0("PC 1, ", round((sdvs/sum(sdvs))[1]*100), "% of variability"))

g_13 <- ggplot(data = pc1_plot_data[pc1_plot_data$clr_trafo == 0,], aes(x = x, y = y, group = score, color = score)) + 
  geom_path(linewidth = 0.5) + theme(legend.position = "none") +
  scale_color_manual(values = c("blue", "black", "red")) + ylab("density") + xlab("rent per sqm") 
################################################################################
g_21 <- ggplot() + geom_sf(data = munich_sf, aes(fill = score2), alpha =0.7) + theme_void() +
  guides(fill=guide_colorbar(title= "score")) + ylab("") + xlab("") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  geom_sf_label(data = munich_sf[-c(3,17),], aes(label = sb_nummer), label.size  = NA, alpha = 0, size = 3) +
  theme(legend.key.height=unit(2,"pt"), legend.key.width = unit(30, "pt"), legend.position = "bottom") 


g_22 <- ggplot(data = pc2_plot_data[pc2_plot_data$clr_trafo == 1,], aes(x = x, y = y, group = score, color = score)) + 
  geom_path(linewidth = 0.5) + theme(legend.position = "top", legend.key.height=unit(4,"pt")) +
  guides(color= guide_legend(title= "score")) + 
  scale_color_manual(values = c("blue", "black", "red")) + ylab("clr density") + xlab("rent per sqm") + 
  ggtitle(paste0("PC 2, ", round((sdvs/sum(sdvs))[2]*100), "% of variability"))

g_23 <- ggplot(data = pc2_plot_data[pc2_plot_data$clr_trafo == 0,], aes(x = x, y = y, group = score, color = score)) + 
  geom_path(linewidth = 0.5) + theme(legend.position = "none") +
  scale_color_manual(values = c("blue", "black", "red")) + ylab("density") + xlab("rent per sqm") 
################################################################################
g_31 <- ggplot() + geom_sf(data = munich_sf, aes(fill = score3), alpha =0.7) + theme_void() +
  guides(fill=guide_colorbar(title= "score")) + ylab("") + xlab("") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  geom_sf_label(data = munich_sf[-c(3,17),], aes(label = sb_nummer), label.size  = NA, alpha = 0, size = 3) +
  theme(legend.key.height=unit(2,"pt"), legend.key.width = unit(30, "pt"), legend.position = "bottom") 


g_32 <- ggplot(data = pc3_plot_data[pc3_plot_data$clr_trafo == 1,], aes(x = x, y = y, group = score, color = score)) + 
  geom_path(linewidth = 0.5) + theme(legend.position = "top", legend.key.height=unit(4,"pt")) +
  guides(color= guide_legend(title= "score")) + 
  scale_color_manual(values = c("blue", "black", "red")) + ylab("clr density") + xlab("rent per sqm") + 
  ggtitle(paste0("PC 3, ", round((sdvs/sum(sdvs))[3]*100), "% of variability"))

g_33 <- ggplot(data = pc3_plot_data[pc3_plot_data$clr_trafo == 0,], aes(x = x, y = y, group = score, color = score)) + 
  geom_path(linewidth = 0.5) + theme(legend.position = "none") +
  scale_color_manual(values = c("blue", "black", "red")) + ylab("density") + xlab("rent per sqm") 
################################################################################
plot_grid(g_12, g_22, g_32, g_13, g_23, g_33, g_11, g_21, g_31, 
          nrow = 3, align = "v", axis = "lr", rel_heights = c(2.5,2,3))

################################################################################
################################################################################
# predict densities
latent_clr_densities <- predict_latent_densities(density_pca)
latent_densities <- lapply(1:25, function(i){
  inverse_clr_trafo(cbind(latent_clr_densities$x_grid, latent_clr_densities$clr_densities[,i]))
})

#plot histogram
old_par <- par(mfrow = c(5,5), mar = c(2, 2, 2, 1))
invisible(lapply(1:25, function(i){
  hist <- hist(rents$nmqm[rents$sb_nummer == i],
       breaks = seq(0,32, by = 2),
       plot = FALSE)
  plot(rbind(latent_densities[[i]], cbind("x" = hist$mids, "y" = hist$density)), cex = 0,
       main = paste0(i, "-", str_trunc(munich_sb$name[i], width = 20)))
  plot(hist, add = TRUE, freq = FALSE, )
  lines(densities_estimated[[i]][,1:2], col = "blue", lwd = 2)
  lines(latent_densities[[i]], col = "red", lwd = 2)
}))
par(old_par)
