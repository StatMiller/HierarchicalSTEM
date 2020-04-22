#Simulate images and run MCMC for hierarchical spatial, and simple linear regression

rm(list = ls())
library(fields)
library(MASS)
library(coda)
library(parallel)
library(foreach)
library(doParallel)

source('boxes_ancillary_functions.R')
source('MCMC.R')
source('atom_grabber_single.R')

#setup
iters <- 20000
model_num <- "sim_sigma220"
nsimul <- 112

ncores <- 16
c1 <- makeCluster(ncores)
registerDoParallel(c1)

simul_boxes <- foreach(i = 1:nsimul, .packages = c("fields", "MASS", "coda")) %dopar% {
  trial <- i
  load("seeds.RData")
  curseed <- seed[i]
  set.seed(curseed)
  sim_num <- i # for saving traceplots
  
  #Number of A- and B-sites
  a <- 18^2 
  b <- 19^2 
  b_dist <- 40 
  
  #locations of b_sites 
  b_locs <- 1:sqrt(b)*b_dist 
  
  
  ##                                      ##
  #  Create true values for the parameters #
  ##
  
  #  Create Mean B-sites                  ##  
  mean_B <- expand.grid(x = b_locs, y = b_locs, KEEP.OUT.ATTRS = FALSE)
  sd_b <- 0.25 
  
  #jitter for actual B_sites - change sd as necessary
  B_sites <- matrix(nrow = b, ncol = 2)
  B_sites[,1] <- mean_B[,1] + rnorm(b, sd=sd_b)
  B_sites[,2] <- mean_B[,2] + rnorm(b, sd=sd_b)
  
  m <- ceiling(max(b_locs) + b_dist)
  n <- m^2
  x = matrix(1:n, m)
  #Organize B_sites as neighbors
  A_neighbors <- list()
  A_neighbors_ind <- list()
  brow <- sqrt(b)
  arow <- sqrt(a)
  for(i in 1:a){
    colnum <- ceiling(i/arow)
    mod <- i%%arow
    if(mod == 0) rownum <- arow
    else rownum <- mod
    ind <- (colnum-1)*brow + rownum
    ind2 <- ind + 1
    ind3 <- ind + brow
    ind4 <- ind3 + 1
    A_neighbors_ind[[i]] <- c(ind, ind2, ind3, ind4)
    A_neighbors[[i]] <- B_sites[A_neighbors_ind[[i]],]
    
  }
  ##Setup for drawing A-sites##
  betaA_mean <- 3060 
  betaB_mean <- 1425 
  beta_sd <- 150
  beta <- rep(NA, a+b+1)
  beta[1] <- 87 
  #A-site betas
  beta[2:(a+1)] <- rnorm(a, mean = betaA_mean, sd = beta_sd)
  #B-site betas
  beta[(a+2):(a+b+1)] <- rnorm(b, mean = betaB_mean, sd = beta_sd)
  
  alpha0 <- -.08  
  alpha1 <- -0.15 
  
  psi_A <- 4.3  
  psi_B <- 3.75 
  sd_a <- 0.4 
  mean_A <-  mean_Asites(A_neighbors, A_neighbors_ind, a, beta, alpha0, alpha1)
  
  d.atoms <- rdist(mean_A$unweighted)
  r <- 0.73 
  rho <- 100 
  Sigma_atoms <- sd_a^2*((1-r)*diag(a) + r*exp(-d.atoms/rho))
  
  A_sites <- matrix(nrow = a, ncol = 2)
  A_sites[,1] <- mvrnorm(1, mu = mean_A$mean_A[,1], Sigma = Sigma_atoms)
  A_sites[,2] <- mvrnorm(1, mu = mean_A$mean_A[,2], Sigma = Sigma_atoms)
  
  #Create X
  A_halfwidth <- 7 
  B_halfwidth <- 6 
  
  A_width = 2*A_halfwidth + 1
  B_width = 2*B_halfwidth + 1
  
  
  A_bins_x = array(NA, c(A_width, A_width, a))
  B_bins_x = array(NA, c(B_width, B_width, b))
  
  for(j in 1:a){
    tmp <- floor(A_sites[j,])
    A_bins_x[,,j] = x[(tmp[1]-A_halfwidth):(tmp[1]+A_halfwidth), 
                      (tmp[2]-A_halfwidth):(tmp[2]+A_halfwidth)]
  }
  
  for(j in 1:b){
    tmp <- floor(B_sites[j,]) 
    B_bins_x[,,j] = x[(tmp[1]-B_halfwidth):(tmp[1]+B_halfwidth), 
                      (tmp[2]-B_halfwidth):(tmp[2]+B_halfwidth)]
  }
  
  
  
  #Get boxes of x-coordinates for A-sites and B-sites
  An = A_width^2
  Bn = B_width^2
  N = An*a + Bn*b
  Ax = array(dim = c(An, 2, a))
  Ax2 = array(dim = c(An, 2, a))
  Bx = array(dim = c(Bn, 2, b))
  Bx2 = array(dim = c(Bn, 2, b))
  Ay = matrix(nrow = An, ncol = a)
  By = matrix(nrow = Bn, ncol = b)
  
  Xcoords.a = array(dim = c(An, 2, a))
  Xcoords.b = array(dim = c(Bn, 2, b))
  
  for(j in 1:a){
    Xcoords.a[,,j] = coord_finder(as.vector(A_bins_x[,,j]), m)
    Ax[,,j] = cbind(1, exp(-rowSums((t(t(Xcoords.a[,,j])-A_sites[j,]))^2)/(2*psi_A^2)))
  }
  
  for(j in 1:b){
    Xcoords.b[,,j] = coord_finder(as.vector(B_bins_x[,,j]), m)
    Bx[,,j] = cbind(1, exp(-rowSums((t(t(Xcoords.b[,,j])-B_sites[j,]))^2)/(2*psi_B^2)))
  }
  
  #Make covariance matrices 
  r_pix <- .57 
  rho_pix <- 5.5 
  sigma2 <- 140^2 
  
  d.a   <- expand.grid(1:A_width-1,1:A_width-1)
  d.a   <- as.matrix(dist(d.a))
  
  V.a <-  (1-r_pix)*diag(An) + r_pix*exp(-d.a/rho_pix)
  varA <- sigma2*V.a
  
  d.b   <- expand.grid(1:B_width-1,1:B_width-1)
  d.b   <- as.matrix(dist(d.b))
  
  V.b <-  (1-r_pix)*diag(Bn) + r_pix*exp(-d.b/rho_pix)
  varB <- sigma2*V.b
  
  #Get responses
  Y_A <- array(dim = c(A_width, A_width, a))
  Y_B <- array(dim = c(B_width, B_width, b))
  
  Y <- matrix(rnorm(m^2, mean = beta[1], sd = 5), nrow = m, ncol = m)
  for(j in 1:a){
    mu <- Ax[,,j]%*%c(beta[1], beta[j+1])
    Y_A[,,j] <- mvrnorm(1, mu = mu, Sigma = varA)
    Y[Xcoords.a[,,j]] <- Y_A[,,j]
  }
  for(j in 1:b){
    mu <- Bx[,,j]%*%c(beta[1], beta[a+j+1])
    Y_B[,,j] <- mvrnorm(1, mu = mu, Sigma = varB)
    Y[Xcoords.b[,,j]] <- Y_B[,,j]
  }
  
  guesses <- atom_grabber(Y_A, Y_B, Xcoords.a, Xcoords.b)
  A_sites_guess <- guesses$A_sites
  B_sites_guess <- guesses$B_sites
  intenseB <- guesses$intenseB
  
  A_neighbors_guess = list()
  for(j in 1:a){
    A_neighbors_guess[[j]] = B_sites_guess[A_neighbors_ind[[j]],]
  }
  
  
  #Bmean Guess
  brow <- sqrt(b)
  startx <- min(B_sites_guess[,1])
  stopx <- max(B_sites_guess[,1])
  starty <- min(B_sites_guess[,2])
  stopy <- max(B_sites_guess[,2])
  xlength <- (stopx - startx)/(brow-1)
  ylength <- (stopy - starty)/(brow-1)
  xgrid <- c(startx, startx + (1:(brow-1))*xlength)
  ygrid <- c(starty, starty + (1:(brow-1))*ylength)
  Bmean_guess <- expand.grid(xgrid, ygrid, KEEP.OUT.ATTRS = FALSE)
  Amean <- matrix(nrow = a, ncol = 2)
  for(i in 1:a){
    Amean[i,] <- colMeans(Bmean_guess[A_neighbors_ind[[i]], ])
  }
  d.guess <- rdist(Amean)
  
  disp <- displacement_bayes_split(intenseB, A_sites_guess, 
                                   A_neighbors_guess, A_neighbors_ind)
  deltax <- disp$x_displacement
  deltay <- disp$y_displacement
  Psix <- disp$psi_x
  Psiy <- disp$psi_y
  
  #Tuning parameters
  MH_psiA <-  .01
  MH_psiB <-  .01
  MH_Asites <-  .05
  MH_Bsites <-  .05
  MH_rho <- .05
  MH_r <- .05
  MH_rhopix <- .02
  MH_rpix <- .02
  
  truth <- c(sqrt(sigma2), psi_A, psi_B, alpha0, alpha1, r_pix, rho_pix, sd_a,
             sd_b, r, rho, betaA_mean, betaB_mean, beta_sd, beta_sd,
             beta, A_sites[,1], A_sites[,2],
             B_sites[,1], B_sites[,2]) #for putting horizontal lines in trace plots
  
  post_adv <-  MCMC_boxes(Y, A_sites_guess, B_sites_guess,  
                          A_neighbors_guess, A_neighbors_ind, Bmean_guess, D = d.guess,
                          sd_Asite = 1, sd_Bsite = 1,
                          iters = iters,
                          MH_rho = MH_rho, MH_r = MH_r,
                          MH_rpix = MH_rpix, MH_rhopix = MH_rhopix,
                          MH_Asites = MH_Asites, MH_Bsites = MH_Bsites,
                          MH_psiA = MH_psiA, MH_psiB = MH_psiB,
                          model_num = model_num, sim_num = sim_num, truth = truth,
                          sim = TRUE, A_halfwidth = 6, B_halfwidth = 5)
  
  post_int <-  MCMC_together(A_sites_guess, B_sites_guess, deltax, deltay,
                             Psix, Psiy, iters = iters)
  
  post_simp <-  MCMC_together_nocorr(A_sites_guess, B_sites_guess, deltax, deltay,
                                     Psix, Psiy, iters=iters)
  
  adv_CI <- HPDinterval(as.mcmc(post_adv$MCMC[(floor(iters/2)+1):iters,]))
  adv_df <- data.frame(mean = colMeans(post_adv$MCMC[(floor(iters/2)+1):iters,]),
                       median = apply(post_adv$MCMC[(floor(iters/2)+1):iters,], 2, median),
                       sd = apply(post_adv$MCMC[(floor(iters/2)+1):iters,], 2, sd),
                       lower = adv_CI[,1], upper = adv_CI[,2])
  
  acceptance <- post_adv$acceptance
  
  int_CI <- HPDinterval(as.mcmc(post_int[(floor(iters/2)+1):iters,]))
  int_df <- data.frame(mean = colMeans(post_int[(floor(iters/2)+1):iters,]),
                       median = apply(post_int[(floor(iters/2)+1):iters,], 2, median),
                       sd = apply(post_int[(floor(iters/2)+1):iters,], 2, sd),
                       lower = int_CI[,1], upper = int_CI[,2])
  
  simp_CI <- HPDinterval(as.mcmc(post_simp[(floor(iters/2)+1):iters,]))
  simp_df <- data.frame(mean = colMeans(post_simp[(floor(iters/2)+1):iters,]),
                        median = apply(post_simp[(floor(iters/2)+1):iters,], 2, median),
                        sd = apply(post_simp[(floor(iters/2)+1):iters,], 2, sd),
                        lower = simp_CI[,1], upper = simp_CI[,2])
  
  assign(paste0("simp_samp", i), post_simp)
  assign(paste0("int_samp", i), post_int)
  assign(paste0("adv_samp", i), post_adv$MCMC)
  
  assign(paste0("model_sum", i), list(adv_df = adv_df, int_df = int_df,
                                      simp_df = simp_df, acceptance = acceptance,
                                      truth = truth, seed = curseed))
  store <- list(get(paste0("simp_samp", i)), 
                get(paste0("int_samp", i)),
                get(paste0("adv_samp", i)),
                get(paste0("model_sum", i)))
  
  store1 <- store[[1]]
  store2 <- store[[2]]
  store3 <- store[[3]]
  store4 <- store[[4]]
  
  save(store4, file = paste0(getwd(), "/sim_results/Model", model_num, "Trial", trial, "summary.RData"))
}
