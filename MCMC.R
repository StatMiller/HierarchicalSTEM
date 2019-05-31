#MCMC for spatial Bayesian hierarchical model for STEM images as described in Miller et al. (2019).
#Inputs defined as:
#Y: STEM image
#A_sites, B_sites: estimated locations for atom column sites
#Mean B: expected location of B-sites based on crystal structure
#D distance matrix for A-sites
#c,d,e,f,g,h hyperparameters for Inverse Gamma distributions
#alpha0, alpha1, intercept and slope initial values for process layer
#iters, number of iterations
#MH_: tuning parameters for respective Metropolis steps
#Mean_range, sd_range: hyperparameters for spatial range parameters rho and rho_pix
#rho, r, tau, rho_pix, r_pix, psi_A, psi_B, sd_Asite, sd_Bsite: initial values for respective parameters
#mean_betaA, mean_BetaB, mean_betaint: initial values that are filled in by OLS estimates if NULL
#sd_alpha0, sd_alpha1, sd_betaint: hyperparameters 
#model_num, sim_num, values that go into naming files at end
#truth: true parameter values in case of simulation to put red horizontal lines in trace plot
#sim: if true, then put red horizontal lines at truth in trace plot

MCMC_boxes <- function(Y, A_sites, B_sites,
                       A_neighbors, A_neighbors_ind, mean_B, D,
                       c = .01, d = .01, e = .01, f = .01, g = .01, h = .01,
                       alpha0 = 0, alpha1 = 0,
                       A_halfwidth = 5, B_halfwidth = 4, 
                       iters = 10000,
                       MH_rhopix=0.005, MH_rpix=0.05, MH_Asites= 0.05, MH_Bsites = 0.05,
                       MH_rho=0.05, MH_r = 0.05,
                       MH_psiA = 0.05, MH_psiB = 0.05,
                       mean_range = 0, sd_range = 10,
                       rho=10,r=0.5,tau=1/.04,
                       rho_pix = 10,r_pix=0.5, 
                       mean_betaA = NULL, mean_betaB = NULL, mean_betaint = 0,
                       sd_betaA=NULL, sd_betaB = NULL, psi_A=5, psi_B = 3,
                       sd_Asite = 1, sd_Bsite = 1, sd_alpha0 = 1000, sd_alpha1 = 1000, sd_betaint = 1000,
                       model_num = 0, sim_num = 0, truth = rep(0,2068), sim =  FALSE){
  
  source("boxes_ancillary_functions.R")
  
  library(fields)
  library(emulator)
  size = nrow(Y)
  x <- matrix(1:size^2, size)
  #set initial values
  sigma2 <- 1/tau
  sdbetaint2 <- sd_betaint^2
  sigma2A <- sd_Asite^2
  sigma2B <- sd_Bsite^2
  
  #Create pixel bins for A and B sites
  a <- nrow(A_sites)
  b <- nrow(B_sites)
  AwB <- AwithB(A_neighbors_ind, b)
  
  A_width <- 2*A_halfwidth + 1
  B_width <- 2*B_halfwidth + 1
  
  A_bins_y = array(NA, c(A_width, A_width, a))
  A_bins_x <- array(NA, c(A_width, A_width, a))
  B_bins_y = array(NA, c(B_width, B_width, b))
  B_bins_x <- array(NA, c(B_width, B_width, b))
   
  for(i in 1:a){
    tmp <- round(A_sites[i,])
    A_bins_y[,,i] = Y[(tmp[1]-A_halfwidth):(tmp[1]+A_halfwidth),
                      (tmp[2]-A_halfwidth):(tmp[2]+A_halfwidth)]
    A_bins_x[,,i] <- x[(tmp[1]-A_halfwidth):(tmp[1]+A_halfwidth),
                       (tmp[2]-A_halfwidth):(tmp[2]+A_halfwidth)]
  }
  
  for(i in 1:b){
    tmp <- round(B_sites[i,])
    B_bins_y[,,i] = Y[(tmp[1]-B_halfwidth):(tmp[1]+B_halfwidth),
                      (tmp[2]-B_halfwidth):(tmp[2]+B_halfwidth)]
    B_bins_x[,,i] <- x[(tmp[1]-B_halfwidth):(tmp[1]+B_halfwidth),
                       (tmp[2]-B_halfwidth):(tmp[2]+B_halfwidth)]
  }
  
  An <- A_width^2
  Bn <- B_width^2
  N <- An*a + Bn*b
  Ax <- array(dim = c(An, 2, a))
  Bx <- array(dim = c(Bn, 2, b))
  Ay <- matrix(nrow = An, ncol = a)
  By <- matrix(nrow = Bn, ncol = b)
  
  Xcoords.a <- array(dim = c(An, 2, a))
  Xcoords.b <- array(dim = c(Bn, 2, b))
  
  #Making the covariance matrix
  d.a   <- expand.grid(1:A_width-1,1:A_width-1)
  d.a   <- as.matrix(dist(d.a))
  
  V.a <- (1-r_pix)*diag(An) + r_pix*exp(-d.a/rho_pix)
  Q.a <- solve(V.a)
  ld.a <- determinant(Q.a)$modulus[1]
  R.a <- matrix(nrow = An, ncol = a)
  
  d.b   <- expand.grid(1:B_width-1,1:B_width-1)
  d.b  <- as.matrix(dist(d.b))
  
  V.b <- (1-r_pix)*diag(Bn) + r_pix*exp(-d.b/rho_pix)
  Q.b <- solve(V.b)
  ld.b <- determinant(Q.b)$modulus[1]
  R.b <- matrix(nrow = Bn, ncol = b)
  
  V.atoms <- (1-r)*diag(a) + r*exp(-D/rho)
  Q.atoms <- solve(V.atoms)
  ld.atoms <- determinant(Q.atoms)$modulus[1]
  
  for(i in 1:a){
    Xcoords.a[,,i] <- coord_finder(as.vector(A_bins_x[,,i]), size)
    Ax[,,i] <- cbind(1, exp(-((Xcoords.a[,1,i]-A_sites[i,1])^2+(Xcoords.a[,2,i]-A_sites[i,2])^2)/(2*psi_A^2)))
    Ay[,i] <- as.vector(A_bins_y[,,i])  
  }
  
  for(i in 1:b){
    Xcoords.b[,,i] <- coord_finder(as.vector(B_bins_x[,,i]), size)
    Bx[,,i] <- cbind(1, exp(-((Xcoords.b[,1,i]-B_sites[i,1])^2+(Xcoords.b[,2,i]-B_sites[i,2])^2)/(2*psi_B^2)))
    By[,i] <- as.vector(B_bins_y[,,i])
  }
  
  #Set initial beta values
  beta <- rep(0,a+b+1)
  
  #Get prior means using OLS
  for(i in 1:a){
    ols <- lm(Ay[,i]~Ax[,2,i])
    beta[1] <- beta[1] + ols$coefficients[1]
    beta[i+1] <- ols$coefficients[2]
  }
  for(i in 1:b){
    ols <- lm(By[,i]~Bx[,2,i])
    beta[1] <- beta[1] + ols$coefficients[1]
    beta[i+a+1] <- ols$coefficients[2]
  }
  beta[1] <- beta[1]/(a+b)
  if(is.null(mean_betaA)) mean_betaA <- mean(beta[2:(a+1)])
  if(is.null(mean_betaB)) mean_betaB <- mean(beta[(a+2):(a+b+1)])
  if(is.null(sd_betaA)) sd_betaA <- sd(beta[2:(a+1)])
  sdbetaA2 <- sd_betaA^2
  if(is.null(sd_betaB)) sd_betaB <- sd(beta[(a+2):(a+b)])
  sdbetaB2 <- sd_betaB^2
  mean_betaA0 <- mean_betaA
  mean_betaB0 <- mean_betaB
  mean_A <- mean_Asites(A_neighbors, A_neighbors_ind, a, beta, alpha0, alpha1)
  sdbetaA_priora <- sdbetaA2/25^2 + 2
  sdbetaA_priorb <- sdbetaA2*(sdbetaA_priora-1)
  sdbetaB_priora <- sdbetaB2/25^2 + 2
  sdbetaB_priorb <- sdbetaB2*(sdbetaB_priora-1)

  
  #Define what we want to report 
  #Parameters - x and y coordinate for each atom  2*(a+b)
  #           - beta for each atom plus intercept (a+b)+1
  #           - pixel intensity variance........        1
  #           - bandwidths .....................        2
  #           - mean and sd's of betas .........        4 
  #           - regression parameters ..........        2
  #           - A-site and B-site variances ....        2
  #           - spatial correlation parameters .        4 
  #   --------------------------------------------------
  #                                             3*(a+b)+16
  keepers <- matrix(0, iters, 3*(a+b)+16)
  colnames(keepers) <- c("sigma","psi_A", "psi_B", "alpha_0", "alpha_1","r_pix","rho_pix",
                         "sigma_A", "sigma_B", "r", "rho", "mean_betaA", "mean_betaB", "sd_betaA", "sd_betaB",
                         paste0("beta",1:(a+b+1)), paste0("A-site xcoord", 1:a),
                         paste0("A-site y-coord", 1:a), paste0("B-site xcoord", 1:b),
                         paste0("B-site y-coord", 1:b))
  
  #keep track of time
  tick <- Sys.time()
  #acceptance rates
  accept_psi_A <- 0
  accept_psi_B <- 0
  accept_Asites <- 0
  accept_Bsites <- 0
  accept_betaB <- 0
  accept_corr <- 0
  accept_corratoms <- 0
  B_den <- 0
  
  #GO!!!
  
  for(j in 1:iters){
    
    #Gibbs sampler for betas at A-sites
    for(i in 1:a){
      
      VVV <- t(Ax[,2,i])%*%Q.a%*%Ax[,2,i]/sigma2 + 1/sdbetaA2
      MMM <- t(Ax[,2,i]%*%Q.a%*%(Ay[,i]-beta[1]))/sigma2 + mean_betaA/sdbetaA2
      MMM <- MMM/VVV
      VVV <- 1/sqrt(VVV)
      
      beta[i+1] <- rnorm(1, MMM, VVV)
      
    }
    
    #Gibbs sampler for betas at B-sites - Do this as Gibbs for candidate, but metropolis overall
    Lambda <- 1/sigma2A*Q.atoms
    for(i in 1:b){
      
      #re did full conditionals 
      VVV <- t(Bx[,2,i])%*%Q.b%*%Bx[,2,i]/sigma2 + 1/sdbetaB2
      MMM <- t(Bx[,2,i]%*%Q.b%*%(By[,i]-beta[1]))/sigma2 + mean_betaB/sdbetaB2
      MMM <- MMM/VVV
      VVV <- 1/sqrt(VVV)
      
      #candidate draw
      beta_can <- beta
      beta_can[a+i+1] <- rnorm(1, MMM, VVV)

      canmean <- By[,i] - (beta[1]*Bx[,1,i]+beta_can[a+i+1]*Bx[,2,i])
      curmean <- By[,i] - (beta[1]*Bx[,1,i]+beta[a+i+1]*Bx[,2,i])

      U1can <- quad.form(Q.b, canmean)
      U1cur <- quad.form(Q.b, curmean)
      U1 <- -.5*tau*(U1can - U1cur)
      
      ind <- AwB[[i]]
      num <- length(ind)
      canA_mean <- mean_Asites(A_neighbors[AwB[[i]]], A_neighbors_ind[AwB[[i]]], a, beta_can,
                               alpha0, alpha1)
      Sigma_tmp <- solve(Lambda[ind,ind])
      matprodx <- Sigma_tmp%*%(Lambda[ind,-ind]%*%as.matrix(A_sites[-ind,1] - mean_A$mean_A[-ind,1]))
      matprody <- Sigma_tmp%*%(Lambda[ind,-ind]%*%as.matrix(A_sites[-ind,2] - mean_A$mean_A[-ind,2]))
      
      canmux <- canA_mean$mean_A[,1] - matprodx
      mux <- mean_A$mean_A[ind,1] - matprodx
      
      canmuy <- canA_mean$mean_A[,2] - matprody
      muy <- mean_A$mean_A[ind,2] - matprody
      
      U2 <- -.5*(quad.form(Lambda[ind,ind], A_sites[ind,1] - canmux) -
                   quad.form(Lambda[ind,ind], A_sites[ind,1] - mux) +
                   quad.form(Lambda[ind,ind], A_sites[ind,2] - canmuy) -
                   quad.form(Lambda[ind,ind], A_sites[ind,2] - muy))

      U3can <- (beta_can[a+i+1] - mean_betaB)^2
      U3cur <- (beta[a+i+1] - mean_betaB)^2
      U3 <- -.5/sdbetaB2*(U3can - U3cur)
      U <- U1 + U2 + U3
      if(is.na(U)) print("NA for betaB first go")
      if(log(runif(1)) < U){
        beta[a+i+1] <- beta_can[a+i+1]
        mean_A$mean_A[ind,] <- canA_mean$mean_A
        accept_betaB <- accept_betaB+1
      }
    }
    
    #Gibbs sampler for beta[1] (intercept)
    VB1 <- (a*sum(Q.a)+b*sum(Q.b))/sigma2+1/sdbetaint2
    
    MB1 <- 0
    for(i in 1:a){
      MB1 <- MB1 + sum(Q.a%*%(Ay[,i] - beta[i+1]*Ax[,2,i]))
    }
    for(i in 1:b){
      MB1 <- MB1 + sum(Q.b%*%(By[,i] - beta[a+i+1]*Bx[,2,i]))
    }
    
    MB1 <- 1/sigma2*MB1 + 1/sdbetaint2*mean_betaint
    MB1 <- MB1/VB1
    VB1 <- 1/sqrt(VB1)
    beta[1] <- rnorm(1, MB1, VB1)  
    
    #Gibbs sampler for sigma^2 - update this to combine info from A and B sites
    for(i in 1:a){
      R.a[,i] <- Ay[,i]-Ax[,,i]%*%c(beta[1], beta[i+1])
    }
    
    for(i in 1:b){
      R.b[,i] <- By[,i]-Bx[,,i]%*%c(beta[1], beta[i+a+1])
    }
    SSA <- sum(quad.diag(Q.a, R.a))
    SSB <- sum(quad.diag(Q.b, R.b))
    tau <- rgamma(1, c + N/2, rate =  d + (SSA+SSB)/2)
    sigma2 <- 1/tau
    
    #Metropolis sampler for correlation parameters
    canrho_pix <- exp(rnorm(1,log(rho_pix),MH_rhopix))
    canr_pix   <- pnorm(rnorm(1,qnorm(r_pix),MH_r))
    
    canQ.a <- solve((1-canr_pix)*diag(An) + canr_pix*exp(-d.a/canrho_pix))
    canQ.b <- solve((1-canr_pix)*diag(Bn) + canr_pix*exp(-d.b/canrho_pix))
    
    canld.a  <- determinant(canQ.a)$modulus[1]
    canld.b  <- determinant(canQ.b)$modulus[1]
    
    canSSA <- sum(quad.diag(canQ.a, R.a))
    canSSB <- sum(quad.diag(canQ.b, R.b))
    
    canU <- dnorm(log(canrho_pix),mean_range, sd_range,log=TRUE)+
      dnorm(qnorm(canr_pix),0,1,log=TRUE)+0.5*(a*canld.a + b*canld.b)-0.5*tau*(canSSA+canSSB)
    curU <- dnorm(log(rho_pix),mean_range, sd_range,log=TRUE)+
      dnorm(qnorm(r_pix),0,1,log=TRUE)+0.5*(a*ld.a + b*ld.b)-0.5*tau*(SSA+SSB)
    
    U <- canU - curU
    prob <- log(runif(1))
    
    if(prob<U){
      rho_pix  <- canrho_pix
      r_pix    <- canr_pix
      Q.a  <- canQ.a
      Q.b  <- canQ.b
      ld.a <- canld.a
      ld.b <- canld.b
      accept_corr <-  accept_corr + 1
    }
    
    #Metropolis sampler for bandwidth 
    #psi_A
    log_canpsi_A <- rnorm(1,log(psi_A), MH_psiA)
    canpsi_A <- exp(log_canpsi_A)
    canAx <- array(dim = c(An, 2, a))
    canR.a <- matrix(nrow = An, ncol = a)
    
    for(i in 1:a){
      canAx[,,i] <- cbind(1, exp(-((Xcoords.a[,1,i]-A_sites[i,1])^2+(Xcoords.a[,2,i]-A_sites[i,2])^2)/(2*canpsi_A^2)))
      
      canR.a[,i] <- Ay[,i] - canAx[,,i]%*%c(beta[1], beta[i+1])
    }
    
    U <- -0.5*tau*sum(quad.diag(Q.a, canR.a)) + 0.5*tau*sum(quad.diag(Q.a,R.a))
    + dnorm(log_canpsi_A, 0, 10, log=TRUE) - dnorm(log(psi_A), 0, 10, log=TRUE)
    
    if(is.na(U)) print("NA for psiA")
    if(log(runif(1))< U){
      
      Ax <- canAx
      R.a <- canR.a
      psi_A <- canpsi_A
      accept_psi_A <- accept_psi_A + 1
    }
    
    #psi_B
    log_canpsi_B <- rnorm(1,log(psi_B), MH_psiB)
    canpsi_B <- exp(log_canpsi_B)
    canBx <- array(dim = c(Bn, 2, b))
    canR.b <- matrix(nrow = Bn, ncol = b)
    
    for(i in 1:b){
      canBx[,,i] <- cbind(1, exp(-((Xcoords.b[,1,i]-B_sites[i,1])^2+(Xcoords.b[,2,i]-B_sites[i,2])^2)/(2*canpsi_B^2)))
      
      canR.b[,i] <- By[,i] - canBx[,,i]%*%c(beta[1], beta[a+i+1])
    }
    
    U <- -0.5*tau*sum(quad.diag(Q.b, canR.b)) + 0.5*tau*sum(quad.diag(Q.b,R.b))
    + dnorm(log_canpsi_B, 0, 10, log=TRUE) - dnorm(log(psi_B), 0, 10, log=TRUE)
    if(is.na(U)) print("NA for betaB PsiB")
    if(log(runif(1))< U){
      Bx <- canBx
      R.b <- canR.b
      psi_B <- canpsi_B
      accept_psi_B <- accept_psi_B + 1
    }
    
    #sd_betaA 
    tau_betaA <- rgamma(1, sdbetaA_priora +a/2, 
                        sdbetaA_priorb + sum((beta[(2:(a+1))]- mean_betaA)^2)/2)
    sdbetaA2 <- 1/tau_betaA
    sd_betaA <- sqrt(sdbetaA2)
    
    #sd_betaB 
    tau_betaB <- rgamma(1, sdbetaB_priora+b/2, 
                        sdbetaB_priorb+ sum((beta[(a+2):(a+b+1)]- mean_betaB)^2)/2)
    sdbetaB2 <- 1/tau_betaB
    sd_betaB <- sqrt(sdbetaB2)

    #mean_betaA
    VA <- 1/1000^2 + a/sdbetaA2
    MA <- mean_betaA0/1000^2 + sum(beta[2:(a+1)])/sdbetaA2
    mean_betaA <- rnorm(1, MA/VA, sqrt(1/VA))
    
    #mean_betaB
    VB <- 1/1000^2 + b/sdbetaB2
    MB <- mean_betaB0/1000^2 + sum(beta[(a+2):(a+b+1)])/sdbetaB2
    mean_betaB <- rnorm(1, MB/VB, sqrt(1/VB))
    
    # A-site locations
    mean_A <- mean_Asites(A_neighbors, A_neighbors_ind, a, beta,
                          alpha0, alpha1)
    
    for(i in 1:a){
      for(k in 1:2){
        can_atom <- A_sites[i,k] + rnorm(1, sd=MH_Asites)
        if(k == 1) can_atom <- c(can_atom, A_sites[i,2])
        else can_atom <- c(A_sites[i,1], can_atom)
        can_X <- cbind(1, exp(-((Xcoords.a[,1,i]-can_atom[1])^2+(Xcoords.a[,2,i]-can_atom[2])^2)/(2*psi_A^2)))
        can_R <- Ay[,i] - can_X%*%c(beta[1], beta[i+1])
        candiff <- can_atom[k] - mean_A$mean_A[i,k]
        curdiff <- A_sites[i,k] - mean_A$mean_A[i,k]
        canU <- -.5*tau*quad.form(Q.a, can_R) - .5/sigma2A*candiff^2
        curU <- -.5*tau*quad.form(Q.a, R.a[,i]) - .5/sigma2A*curdiff^2

        U <- canU - curU
        if(is.na(U)) print("NA for A-sites")
        if(log(runif(1)) < U){
          A_sites[i,] <- can_atom
          Ax[,,i] <- can_X
          R.a[,i] <- can_R
          accept_Asites <- accept_Asites + 1
        }
      }
    }
    
    D <- rdist(A_sites)
    V.atoms <- (1-r)*diag(a) + r*exp(-D/rho)
    Q.atoms <- solve(V.atoms)
    ld.atoms <- determinant(Q.atoms)$modulus[1]
    
    ## B-site locations ##
    for(i in 1:b){
      ind <- AwB[[i]]
      for(k in 1:2){
        can_atom <- B_sites[i,k] + rnorm(1, sd=MH_Bsites)
        if(k == 1) can_atom <- c(can_atom, B_sites[i,2])
        else can_atom <- c(B_sites[i,1], can_atom)
        can_X <- cbind(1, exp(-((Xcoords.b[,1,i]-can_atom[1])^2+(Xcoords.b[,2,i]-can_atom[2])^2)/(2*psi_B^2)))
        
        can_R <- By[,i] - can_X%*%c(beta[1], beta[a+i+1])
        
        candiff <- can_atom[k] - mean_B[i,k]
        canU <- -.5*tau*quad.form(Q.b, can_R) - .5/sigma2B*candiff^2
        
        curdiff <- B_sites[i,k] - mean_B[i,k]
        curU <- -.5*tau*quad.form(Q.b, R.b[,i]) - .5/sigma2B*curdiff^2
        
        U1 <- canU - curU
        
        canB_sites <- B_sites
        canB_sites[i,] <- can_atom
 
        canA_neighbors <- A_neighbors
        for(l in ind){
          canA_neighbors[[l]] <- canB_sites[A_neighbors_ind[[l]],]
        }
        canA_mean <- mean_Asites(canA_neighbors[ind], A_neighbors_ind[ind], a, beta,
                                 alpha0, alpha1)
        Sigma_tmp <- solve(Lambda[ind,ind])
        
        if(k ==1){
          matprodx <- Sigma_tmp%*%(Lambda[ind,-ind]%*%as.matrix(A_sites[-ind,1] - mean_A$mean_A[-ind,1]))
          canmux <- canA_mean$mean_A[,1] - matprodx
          mux <- mean_A$mean_A[ind,1] - matprodx
          U2 <- -.5*(quad.form(Lambda[ind,ind], A_sites[ind,1] - canmux) -
                       quad.form(Lambda[ind,ind], A_sites[ind, 1] - mux))
        }else{
          matprody <- Sigma_tmp%*%(Lambda[ind,-ind]%*%as.matrix(A_sites[-ind,2] - mean_A$mean_A[-ind,2]))
          canmuy <- canA_mean$mean_A[,2] - matprody
          muy <- mean_A$mean_A[ind,2] - matprody
          U2 <- -.5*(quad.form(Lambda[ind,ind], A_sites[ind,2] - canmuy) - 
                       quad.form(Lambda[ind,ind], A_sites[ind,2] - muy))
        }
        
        U <- U1+U2
        if(is.na(U)) print("NA for B_sites")
        if(log(runif(1)) < U){
          B_sites[i,] <- can_atom
          Bx[,,i] <- can_X
          R.b[,i] <- can_R
          mean_A$mean_A[ind,] <- canA_mean$mean_A
          accept_Bsites <- accept_Bsites + 1
          A_neighbors <- canA_neighbors
        }
      }
    }
    
    #Update A_neighbors
    mean_A <- mean_Asites(A_neighbors, A_neighbors_ind, a, beta,
                          alpha0, alpha1) 
    
    ##sigma2A##
    Rx <- A_sites[,1] - mean_A$mean_A[,1]
    Ry <- A_sites[,2] - mean_A$mean_A[,2]
    SSx <- quad.form(Q.atoms, Rx)
    SSy <- quad.form(Q.atoms, Ry)
    tau2A <- rgamma(1,e+a, f+(SSx+SSy)/2)
    sigma2A <- 1/tau2A 
    sd_Asite <- sqrt(sigma2A) 
    
    Sigmainv <- Q.atoms/sigma2A
    deltax <- A_sites[,1] - mean_A$unweighted[,1]
    deltay <- A_sites[,2] - mean_A$unweighted[,2]
    Psix <- mean_A$weighted[,1] - mean_A$unweighted[,1]
    Psiy <- mean_A$weighted[,2] - mean_A$unweighted[,2]
    
    A0 <- 2*sum(Sigmainv) + 1/sd_alpha0^2
    B0 <- sum(Sigmainv%*%(deltax + deltay - alpha1*(Psix + Psiy)))
    alpha0 <- rnorm(1, B0/A0, sqrt(1/A0))
    
    #alpha1
    A1 <- quad.form(Sigmainv, Psix) + quad.form(Sigmainv, Psiy) +
      1/sd_alpha1^2
    B1 <- t(Psix)%*%Sigmainv%*%(deltax - alpha0) + 
      t(Psiy)%*%Sigmainv%*%(deltay - alpha0)
    
    alpha1 <- rnorm(1, B1/A1, sqrt(1/A1))
    
    # Correlation Parameters for A-sites 
    canrho <- exp(rnorm(1,log(rho),MH_rho))
    canr   <- pnorm(rnorm(1,qnorm(r),MH_r))
    
    canV.atoms <- (1-canr)*diag(a) + canr*exp(-D/canrho)
    canQ.atoms <- solve(canV.atoms)
    
    canld.atoms  <- determinant(canQ.atoms)$modulus[1]
    
    canSSx <- quad.form(canQ.atoms, Rx)
    canSSy <- quad.form(canQ.atoms, Ry)
    
    canU <- dnorm(log(canrho),mean_range, sd_range,log=TRUE)+
      dnorm(qnorm(canr),0,1,log=TRUE)+canld.atoms-0.5*tau2A*(canSSx+canSSy)
    
    curU <- dnorm(log(rho),mean_range, sd_range,log=TRUE)+
      dnorm(qnorm(r),0,1,log=TRUE)+ld.atoms-0.5*tau2A*(SSx+SSy)
    
    U <- canU - curU
    prob <- log(runif(1))
    if(prob<U){
      rho  <- canrho
      r    <- canr
      V.atoms <- canV.atoms
      Q.atoms  <- canQ.atoms
      ld.atoms <- canld.atoms
      accept_corratoms <- accept_corratoms + 1
    }
    
    # B-site variance 
    shapeB <- b + g
    differenceB <- B_sites - mean_B
    sum4 <- sum(differenceB[,1]^2+differenceB[,2]^2)
    rateB <- sum4/2+h
    sigma2B <- 1/rgamma(1, shape = shapeB, rate = rateB)
    sd_Bsite <- sqrt(sigma2B)
    
    keepers[j,] <- c(1/sqrt(tau),
                     psi_A, psi_B, alpha0, alpha1, r_pix,rho_pix, sd_Asite, sd_Bsite,
                     r, rho, mean_betaA, mean_betaB, sd_betaA, sd_betaB,
                     beta, A_sites[,1], A_sites[,2], B_sites[,1], B_sites[,2])
    # Report out time every 100 iterations
    if(j%%100 == 0){
      tock <- Sys.time()
      print(paste(j,"out of", iters, "iterations complete in", tock - tick))
    }
  }
  
  #Save trace plots and results
  if(sim){ 
    filename <- paste0(getwd(), "/traceplots/", "Model", model_num, "Sim", sim_num, ".pdf")
  } else filename <- paste0(getwd(), "/traceplots/", "Model", model_num, ".pdf")
  pdf(filename)
  par(mfrow=c(4,4))
  for(j in 1:ncol(keepers)){
    if(sim){
      ymin <- min(truth[j], min(keepers[,j]))
      ymax <- max(truth[j], max(keepers[,j]))
      plot(1:iters, keepers[,j],type="l", xlab="MCMC iteration",ylab="Sample", 
           main=colnames(keepers)[j], ylim = c(ymin, ymax))
      abline(h = truth[j], col = "red")
    }
    else plot(1:iters, keepers[,j],type="l", xlab="MCMC iteration",ylab="Sample", 
              main=colnames(keepers)[j])
  }
  dev.off()
  
  acceptance <- c(accept_corr/iters*100, accept_corratoms/iters*100, accept_psi_A/iters*100, accept_psi_B/iters*100,
                  accept_Asites/(2*a*iters)*100, accept_Bsites/(2*b*iters)*100, accept_betaB/(b*iters+B_den)*100)
  names(acceptance) <- c("correlation", "correlation_atoms", "psi_A", "psi_B",
                         "A_sites", "B_sites", "beta_B")
  
  
  return(list(MCMC = keepers, acceptance = acceptance))
}

