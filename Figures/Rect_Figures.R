#Draw atom column location estimates and windows around 
   #Columns
rm(list = ls())
source("boxes_ancillary_functions.R")
load("AtomInfo.RData")
A_halfwidth <- 6
B_halfwidth <- 5
rectA <- round(A_sites)
rectB <- round(B_sites)

m <- nrow(Y)
par(pty = "s")
lwd <- 1
image(1:m, 1:m, Y, axes = FALSE, xlab = "", ylab = "",
      col = grey(seq(0,1,length = 256)))
rectcoordsA = cbind(rectA[,1] - A_halfwidth, rectA[,2] - A_halfwidth, rectA[,1] + A_halfwidth, rectA[,2] + A_halfwidth)
rectcoordsB = cbind(rectB[,1] - B_halfwidth, rectB[,2] - B_halfwidth, rectB[,1] + B_halfwidth, rectB[,2] + B_halfwidth)

rect(rectcoordsA[,1], rectcoordsA[,2], rectcoordsA[,3], rectcoordsA[,4], border = "red", lwd = lwd)
rect(rectcoordsB[,1], rectcoordsB[,2], rectcoordsB[,3], rectcoordsB[,4], border = "blue", lwd = lwd)

indA <- c(17:18, 35:36)
indB <- c(17:19, 36:38, 55:57)
lwd <- 2

image(452:514, 45:107, Y[452:514, 45:107], axes = FALSE, xlab = "", ylab = "",
      col = grey(seq(0,1,length = 256)))
rect(rectcoordsA[indA,1], rectcoordsA[indA,2], rectcoordsA[indA,3], rectcoordsA[indA,4], border = "red", lwd = lwd)
rect(rectcoordsB[indB,1], rectcoordsB[indB,2], rectcoordsB[indB,3], rectcoordsB[indB,4], border = "blue", lwd = lwd)
points(round(A_sites[indA,]), col = "red", pch = 19)
points(round(B_sites[indB,]), col = "blue", pch = 19)


