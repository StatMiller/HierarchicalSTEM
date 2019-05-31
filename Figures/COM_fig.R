#Process layer figure
B_sites <- cbind(c(250,250,750,750),c(250,750,250,750))
Unweighted <- colMeans(B_sites)
Intensity <- c(4,1,1,1)
color <- grey(seq(0,1,length = 256))
cols <- c(color[1], color[128], color[128], color[128])
Weighted <- colSums(B_sites*Intensity)/sum(Intensity)
alpha_1 <- -1
A_site = Unweighted + alpha_1*(Weighted - Unweighted)

par(pty = "s", xaxt = "n", yaxt = "n")
#par(xaxt = "n", yaxt = "n")
plot(B_sites, xlim = c(200,800), ylim = c(200,800)
     ,pch = 19, cex = 10, xlab = "", ylab = "", col = cols)
points(Unweighted[1], Unweighted[2], pch = 3, cex = 5)
points(Weighted[1], Weighted[2], pch = 8, cex = 5)
points(A_site[1], A_site[2], pch = 1, cex = 12)
legend(x = 200, y = 200, bty = "n", xpd = T,
       legend = c("B-sites", "A-site", "Unweighted Mean B", "Weighted Mean B"), ncol = 2,
       pch = c(19, 1, 3, 8), cex = 1.5,text.width=100, x.intersp = .5)
labs <- c(expression(s[B1]), expression(s[B2]), expression(s[B3]), expression(s[B4]))
text(B_sites[,1], B_sites[,2], labs, col = "white", cex =3)
text(A_site[1],A_site[2], expression(s[A1]), cex = 3)
text(Unweighted[1], Unweighted[2]-55, expression(U[A1]), cex = 3)
text(Weighted[1], Weighted[2]-55, expression(W[A1]), cex = 3)