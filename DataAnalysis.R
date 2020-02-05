#Run data analysis on STEM data using Hieararchical model, Spatial linear regression, 
#   and simple linear regression


rm(list = ls())
library(fields)
source("MCMC_final.R")
load('AtomInfo.RData')
a <- nrow(A_sites)
b <- nrow(B_sites)

#Make mean B-sites and A-site initial estimates
brow <- sqrt(b)
startx <- min(B_sites[,1])
stopx <- max(B_sites[,1])
starty <- min(B_sites[,2])
stopy <- max(B_sites[,2])
xlength <- (stopx - startx)/(brow-1)
ylength <- (stopy - starty)/(brow-1)
xgrid <- c(startx, startx + (1:(brow-1))*xlength)
ygrid <- c(starty, starty + (1:(brow-1))*ylength)
Bmean <- expand.grid(xgrid, ygrid, KEEP.OUT.ATTRS = FALSE)

Amean <- matrix(nrow = a, ncol = 2)
for(i in 1:a)  Amean[i,] <- colMeans(Bmean[A_neighbors_ind[[i]], ])

D <- rdist(Amean)
iters <- 100000
model_num = "1_100kiters"
post_adv <- MCMC_boxes(Y, A_sites, B_sites, A_neighbors, A_neighbors_ind, Bmean, D = D,
                        sd_Asite = 1, sd_Bsite = 1,
                        iters = iters, MH_rho = .1, MH_r = .1, MH_rpix = .02, MH_rhopix = .02,
                        MH_Asites = .05, MH_Bsites = .05, MH_psiA = .01, MH_psiB = .01, model_num = model_num)

save(post_adv, file = paste0("post_adv_",model_num, "iters_", iters,  format(Sys.time(), "_%Y%m%d_%H%M_"),".RData"))

#For simple and spatial linear regression:
source("boxes_ancillary_functions.R")
D2 <- exp(-rdist(Amean)) #These functions require different distance matrix ('boxes_ancillary_functions.R')
post_int <- MCMC_together(A_sites, B_sites, deltax, deltay,
                          Psix, Psiy, D = D2, iters = iters, report_acceptance = TRUE)
post_simp <- MCMC_together_nocorr(A_sites, B_sites, deltax, deltay,
                                  Psix, Psiy, iters = iters)
save(post_int, file = paste0("post_int_", "iters_", iters,  format(Sys.time(), "_%Y%m%d_%H%M_"),".RData"))

save(post_simp, file = paste0("post_simp_", "iters_", iters,  format(Sys.time(), "_%Y%m%d_%H%M_"),".RData"))
