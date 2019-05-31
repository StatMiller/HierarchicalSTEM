#Find the correct coordinates for each pixel
coord_finder = function(X, m){
  xcoord = X%%m
  xcoord[which(xcoord == 0)] = m
  ycoord = ceiling(X/m)
  
  if(is.vector(X)){return(cbind(xcoord,ycoord))}
  else{return(array(c(xcoord,ycoord), dim = c(dim(xcoord),2)))}
}


#Define expit and logit
expit = function(x){return(1/(1+exp(-x)))}
logit = function(x){return(log(x/(1-x)))}

#find unweighted and weighted A-sites to get the mean A-sites
mean_Asites <-  function(A_neighbors, A_neighbors_ind, a, beta=NULL, alpha0, alpha1){
  n <- length(A_neighbors)
  mean_A <- matrix(nrow = n, ncol = 2)
  unweighted <- matrix(nrow = n, ncol = 2)
  weighted <- matrix(nrow = n, ncol = 2)
  for(i in 1:n){
    B <- A_neighbors[[i]]
    unweighted[i,] <- colMeans(B)
    num <- colSums(beta[1+a+A_neighbors_ind[[i]]]*B)
    den <- sum(beta[1+a+A_neighbors_ind[[i]]])
    weighted[i,] <- num/den
    mean_A[i,] <- rep(alpha0,2) + (1-alpha1)*unweighted[i,] + alpha1*weighted[i,]
  }
  
  return(list(mean_A = mean_A, unweighted = unweighted, weighted = weighted))
}

#Get the A-sites associated with a particular B-site for metropolis step of beta_B's
AwithB <- function(A_neighbors_ind, b){
  nei <- list()
  n <- length(A_neighbors_ind)
  for(i in 1:b){
    tmp <- NULL
    for(j in 1:n){
      if(any(A_neighbors_ind[[j]] == i)) tmp <- c(tmp, j)
    }
    nei[[i]] <- tmp
  }
  return(nei)
}

#Spatial linear regression MCMC
MCMC_together <- function(A_sites, B_sites, deltax, deltay,
                 Psix, Psiy, D = NULL,
                 c = .01, d = .01, 
                 alpha0 = 0, alpha1 = 0,
                 iters = 10000, burn = 1000,
                 MH_rho=0.05, MH_r=0.05, 
                 mean_range = 0, sd_range = 10,
                 rho=10,r=0.5,tau=1/.04,
                 sd_alpha0 = 1000, sd_alpha1 = 1000,
                 report_acceptance =FALSE){
  library(emulator)
  library(fields)
  
  a <- nrow(A_sites)
  if(is.null(D)) D <- exp(-rdist(A_sites)) 
  V <- (1-r)*diag(a) + r*D^(1/rho)
  nu2 <- 1/tau
  Q <- solve(V)
  ld <- determinant(Q)$modulus[1]
  
  #Make keepers
  keepers <-  matrix(0, iters, 5)
  colnames(keepers) <-  c("alpha0", "alpha1", "nu", "r", "rho")
  
  accept_corr <- 0
  tick <- Sys.time()
  
  #Start MCMC
  for(i in 1:iters){
    
    # Gibbs samplers for regression coefficients
    Sigmainv <- Q/nu2
    A0 <- 2*sum(Sigmainv) + 1/sd_alpha0^2
    B0 <- sum(Sigmainv%*%(deltax + deltay - alpha1*(Psix + Psiy)))
    
    alpha0 <- rnorm(1, B0/A0, sqrt(1/A0))
    
    A1 <- quad.form(Sigmainv, Psix) + quad.form(Sigmainv, Psiy) +
      1/sd_alpha1^2
    B1 <- t(Psix)%*%Sigmainv%*%(deltax - alpha0) + 
      t(Psiy)%*%Sigmainv%*%(deltay - alpha0)
    
    alpha1 <- rnorm(1, B1/A1, sqrt(1/A1))
    
    # Gibbs sampler for variance
    mux <- alpha0 + alpha1*Psix
    SSx <- quad.form(Q, deltax - mux)
    muy <- alpha0 + alpha1*Psiy
    SSy <- quad.form(Q, deltay - muy)
    SS <- SSx+SSy
    
    tau <- rgamma(1, c + a, d + SS/2)
    nu2 <- 1/tau
    nu <- sqrt(nu2)
    
    #Metropolis sampler for rho, r
    canrho <- exp(rnorm(1,log(rho),MH_rho))
    canr   <- pnorm(rnorm(1,qnorm(r),MH_r))
    canV <- (1-canr)*diag(a) + canr*D^(1/canrho)
    canQ <- solve(canV)
    canSS <- quad.form(canQ, deltax - mux) + quad.form(canQ, deltay - muy)
    canld <- determinant(canQ)$modulus[1]
    
    canU <- dnorm(log(canrho),mean_range, sd_range,log=TRUE)+
      dnorm(qnorm(canr),0,1,log=TRUE)+0.5*(2*canld)-0.5*tau*(canSS)
    
    orgU <- dnorm(log(rho),mean_range, sd_range,log=TRUE)+
      dnorm(qnorm(r),0,1,log=TRUE)+0.5*(2*ld)-0.5*tau*(SS)
    U <- canU - orgU
    prob = log(runif(1))
    
    if(prob<U){
      rho  <- canrho
      r    <- canr
      V <- canV
      Q <- canQ
      ld <- canld
      accept_corr = accept_corr + 1
    }
    
    if(i%%1000 == 0){
      tock = Sys.time()
      print(paste(i,"out of", iters, "iterations complete in", tock - tick))
    }
    
    keepers[i,] <- c(alpha0, alpha1, nu, r, rho)
    
  }
  
  if(report_acceptance) return(list(posterior = keepers, accept_per = accept_corr/iters*100))
  else return(keepers)
}

#Simple linear regression MCMC
MCMC_together_nocorr <- function(A_sites, B_sites, deltax, deltay,
                          Psix, Psiy,
                          c = .01, d = .01, 
                          alpha0 = 0, alpha1 = 0,
                          iters = 10000, burn = 1000,
                          tau=1/.04, sd_alpha0 = 1000, sd_alpha1 = 1000){
  library(emulator)
  library(fields)
  
  a <- nrow(A_sites)
  nu2 <- 1/tau

  
  #Make keepers
  keepers <-  matrix(0, iters, 3)
  colnames(keepers) <-  c("alpha0", "alpha1", "tau_A")
  
  accept_corr <- 0
  tick <- Sys.time()
  
  #Start MCMC
  for(i in 1:iters){
    
    # Gibbs samplers for regression coefficients
    A0 <- 2*a/nu2 + 1/sd_alpha0^2
    B0 <- 1/nu2*sum(c(deltax,deltay) - alpha1*c(Psix, Psiy))
    
    alpha0 <- rnorm(1, B0/A0, sqrt(1/A0))
    
    A1 <- sum(c(Psix,Psiy)^2)/nu2 +  1/sd_alpha1^2
    B1 <- 1/nu2*(t(Psix)%*%(deltax - alpha0) + 
      t(Psiy)%*%(deltay - alpha0))
    
    alpha1 <- rnorm(1, B1/A1, sqrt(1/A1))
    
    # Gibbs sampler for variance
    mux <- alpha0 + alpha1*Psix
    SSx <- sum((deltax-mux)^2)
    muy <- alpha0 + alpha1*Psiy
    SSy <- sum((deltay - muy)^2)
    SS <- SSx+SSy
    
    tau <- rgamma(1, c + a, d + SS/2)
    nu2 <- 1/tau
    nu <- sqrt(nu2)
    tauA <- nu
    
    if(i%%1000 == 0){
      tock = Sys.time()
      print(paste(i,"out of", iters, "iterations complete in", tock - tick))
    }
    
    keepers[i,] <- c(alpha0, alpha1, tauA)
    
  }
  return(keepers)
}

#Calculate displacement and weighted and unweighted intensties 
displacement_bayes_split = function(intenseB, A_sites, 
                                    A_neighbors, A_neighbors_ind){
  n <- nrow(A_sites)
  neighbor_num <- nrow(A_neighbors[[1]])
  
  Ux <- rep(NA, n)
  Uy <- rep(NA, n)
  Wx <- rep(NA, n)
  Wy <- rep(NA, n)
  intensity <- matrix(nrow = n, ncol = neighbor_num)
  
  for(i in 1:n){
    unweighted <- colMeans(A_neighbors[[i]])
    Ux[i] <- unweighted[1]
    Uy[i] <- unweighted[2]
    
    intensity[i,] <- intenseB[A_neighbors_ind[[i]]]
    Wx[i] <- (intensity[i,]%*%A_neighbors[[i]][,1])/sum(intensity[i,])
    Wy[i] <- (intensity[i,]%*%A_neighbors[[i]][,2])/sum(intensity[i,])
  }
  
  x_displacement <- A_sites[,1]-Ux
  y_displacement <- A_sites[,2]-Uy
  psi_x <- Wx-Ux
  psi_y <- Wy-Uy
  
  return(data.frame(x_displacement = x_displacement, y_displacement = y_displacement,
                    psi_x = psi_x, psi_y = psi_y))
}

#Plot STEM image with A- and B-sites
STEMplot <- function(STEM, Aguess, Bguess){
  m <- nrow(STEM)
  par(pty = "s")
  image(1:m, 1:m, STEM, axes = FALSE, xlab = "", ylab = "",
        col = grey(seq(0,1,length = 256)))
  points(Aguess, col = "red")
  points(Bguess, col = "blue")
}
