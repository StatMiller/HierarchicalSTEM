#Density plot figures
rm(list = ls())
library(ggplot2)
library(reshape2)
library(wesanderson)

#Load in the MCMC data
load("post_simp_iters_1e+05_20190416_2347_.RData")
load("post_int_iters_1e+05_20190416_2347_.RData")
load("post_adv_2_100kitersiters_1e+05_20190413_0935_.RData")

afterburn <- 10001:100000
simp_thin <- post_simp[afterburn,]
int_thin <- post_int$posterior[afterburn,]
adv_thin <- post_adv2$MCMC[afterburn,]
rm(post_adv2)
rm(post_int)
rm(post_simp)

x <- data.frame(SimpleLR = simp_thin[,1], SpatialLR = int_thin[,1], BayesHier = adv_thin[,4])
data <- melt(x)

p1 <- ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.5, color = "NA") +
  labs(x = expression(alpha[0]), y = "Density", fill = "Model:")  + theme_bw() + 
  theme(text = element_text(size = 20)) + 
  scale_fill_manual(values=wes_palette(n=3, name="BottleRocket2")) + xlim(-.4,.4)+ 
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), legend.position = "none")

x <- data.frame(SimpleLR = simp_thin[,2], SpatialLR = int_thin[,2], BayesHier = adv_thin[,5])
data <- melt(x)

p2 <- ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.5, color = "NA") +
  labs(x = expression(alpha[1]), y = "Density", fill = "Model:")  + theme_bw() + 
  theme(text = element_text(size = 20)) + 
  scale_fill_manual(values=wes_palette(n=3, name="BottleRocket2"))+ 
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), legend.position = "none")

x <- data.frame(SimpleLR = simp_thin[,3], SpatialLR = int_thin[,3], BayesHier = adv_thin[,8])
data <- melt(x)

p3 <- ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.5, color = "NA") +
  labs(x = expression(sigma["A"]), y = "Density", fill = "Model:")  + theme_bw() + 
  theme(text = element_text(size = 20)) + 
  scale_fill_manual(values=wes_palette(n=3, name="BottleRocket2")) + xlim(0.25,.7)+ 
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), legend.position = "none")


x <- data.frame(SpatialLR = int_thin[,4], BayesHier = adv_thin[,10])
data <- melt(x)

p4 <- ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.5, color = "NA") +
  labs(x = expression(r), y = "Density", fill = "Model:")  + theme_bw() + 
  theme(text = element_text(size = 20)) + 
  scale_fill_manual(values=wes_palette(n=3, name="BottleRocket2")[2:3])+ 
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), legend.position = "none")

x <- data.frame(SpatialLR = int_thin[,5], BayesHier = adv_thin[,11])
data <- melt(x)

p5 <- ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.5, color = "NA") +
  labs(x = expression(rho), y = "Density", fill = "Model:")  + theme_bw() + 
  theme(text = element_text(size = 20)) + 
  scale_fill_manual(values=wes_palette(n=3, name="BottleRocket2")[2:3]) + xlim(0,500) + 
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), legend.position = "none")

library(ggpubr)
ggarrange(p1,p2,p3,p4,p5,ncol=3,nrow=2,common.legend=TRUE, legend="bottom")
        