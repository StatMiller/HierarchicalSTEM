#Posterior Ellispse Drawer
rm(list = ls())
library(car)
source("boxes_ancillary_functions.R")
load("AtomInfo.RData")

#Load in MCMC results for real data
load("post_adv_2_100kitersiters_1e+05_20190413_0935_.RData")

m <- nrow(Y)
a <- nrow(A_sites)
b <- nrow(B_sites)
par(pty = "s")
image(452:514, 45:107, Y[452:514, 45:107], axes = FALSE, xlab = "", ylab = "",
      col = grey(seq(0,1,length = 256)))
afterburn <- 10001:100000
adv <- post_adv2$MCMC[afterburn,]

rm(post_adv2)

adv_mean <- colMeans(adv)
indA <- c(17:18, 35:36)
indB <- c(17:19, 36:38, 55:57)
Amean_coord <- round(matrix(adv_mean[702:(702+2*a-1)], ncol = 2))
Bmean_coord <- round(matrix(adv_mean[(702+2*a):length(adv_mean)], ncol = 2))


Apostx <- adv[,(702:(702+2*a-1))[indA]]
Aposty <- adv[,(702:(702+2*a-1))[indA+a]]
Bpostx <- adv[,((702+2*a):length(adv_mean))[indB]]
Bposty <- adv[,((702+2*a):length(adv_mean))[indB+b]]

par(pty = "s")
image(452:514, 45:107, Y[452:514, 45:107], axes = FALSE, xlab = "", ylab = "",
      col = grey(seq(0,1,length = 256)))
for(i in 1:4){
 Acoords <- cbind(Apostx[,i], Aposty[,i])
 lines(mixtools::ellipse(mu = round(colMeans(Acoords)), sigma = cov(Acoords)), type = 'l', col = "red")
 points(Amean_coord, col = "red", pch = '.')
}

for(i in 1:9){
  Bcoords <- cbind(Bpostx[,i], Bposty[,i])
  lines(mixtools::ellipse(mu = round(colMeans(Bcoords)), sigma = cov(Bcoords)), type = 'l', col = "blue")
  points(Bmean_coord, col = "blue", pch = '.')
}
