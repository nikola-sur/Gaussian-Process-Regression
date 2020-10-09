# STAT 850 Final Project R Code
# Nikola Surjanovic, Fall 2019

rm(list=ls())
library(geoR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(sp)
library(gstat)
library(scales) # for "comma"
library(magrittr)
# library(devtools)
# install_version("sp", version = "0.8-7", repos = "http://cran.us.r-project.org")
# install_version("geoR", version = "1.6-11", repos = "http://cran.us.r-project.org")


# Motivation for kriging (1-dimensional) ---------------------
set.seed(214)
sim11 <- grf(1000, grid="reg", ny=1, cov.pars=c(1, 0.25), nug=0, cov.model="matern", kappa=0.5, mean=3)
# sim11 <- grf(1000, grid="reg", ny=1, cov.pars=c(1, 0.25), nug=0, cov.model="matern", kappa=2.5, mean=3)
image(sim11, type="l", ylim=range(sim11$data), cex.lab = 1.4, cex.axis=1.4)

sim11 %>% 
  as.data.frame() %>%
  ggplot(aes(x=x, y=data)) +
  geom_line() + ylab('y') +
  theme(axis.text=element_text(size=24),
       axis.title=element_text(size=24,face="bold"))



inds <- seq(1, 1000, 100)
temp_data <- data.frame(x=sim11[[1]][inds, 1], y=sim11$data[inds])
plot(x=temp_data$x, y=temp_data$y, col='black', type='p', lwd='6',
     main="Observed Gold Concentrations", xlab="Location", ylab="Concentration")
plot(x=temp_data$x, y=temp_data$y, col='black', type='p', lwd='6',
     main="Observed Gold Concentrations", xlab="Location", ylab="Concentration", ylim=c(1, 5.5))
lines(x=sim11[[1]][ ,1], y=sim11$data)
temp_lm <- lm(y ~ x, data=temp_data)
abline(temp_lm, col='red', lwd=3)

# Fitting a 1-D krige
gold_data <- cbind(temp_data, 1)
gold_data <- gold_data[ ,c(1,3,2)]
gold_data <- as.geodata(gold_data)
gold_ml <- likfit(gold_data, cov.model="matern", kappa=1.5, ini=c(0.9, 0.2), nugget=0.2, fix.nugget=FALSE)
gold_ml2 <- likfit(gold_data, cov.model="matern", kappa=1.5, ini=c(0.9, 0.2), nugget=0.5, fix.nugget=TRUE)

gold_coords <- seq(0, 1, by=0.001)
gold_gr <- pred_grid(coords = gold_coords, y.coords=1, by=0.001)
gold_KC <- krige.control(obj.model = gold_ml)
gold_KC2 <- krige.control(obj.model = gold_ml2)
gold_OC <- output.control(n.pred = 10, simul = TRUE)
gold_pred <- krige.conv(gold_data, locations = gold_gr, krige=gold_KC)
set.seed(2838321)
gold_pred_OC <- krige.conv(gold_data, locations = gold_gr, krige=gold_KC, out=gold_OC)
gold_pred2 <- krige.conv(gold_data, locations = gold_gr, krige=gold_KC2)

plot(x=temp_data$x, y=temp_data$y, col='black', type='p', lwd='6',
     main="Observed Gold Concentrations", xlab="Location", ylab="Concentration")
lines(x=gold_coords[-1], y=as.numeric(gold_pred$predict[-1]), lwd=3, col='red')

plot(x=temp_data$x, y=temp_data$y, col='black', type='p', lwd='6',
     main="Observed Gold Concentrations", xlab="Location", ylab="Concentration")
lines(x=gold_coords[-1], y=as.numeric(gold_pred2$predict[-1]), lwd=3, col='red')

krige_curve_df <- data.frame(x=gold_coords[-1],
                             curve1=as.numeric(gold_pred$predict[-1]),
                             curve2 = as.numeric(gold_pred2$predict[-1]))
krige_points_df <- data.frame(x=temp_data$x, y=temp_data$y)


ggplot(krige_curve_df, aes(x=x)) +
  theme(legend.position = 'none',
        axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold")) + 
  geom_point(data=krige_points_df, aes(x=x,y=y), size=3.5) + 
  ylab('y') + ylim(c(1, 5.5)) + xlim(c(0,1))

ggplot(krige_curve_df, aes(x=x)) +
  geom_line(aes(y=curve1, colour='red'), size=1.5) +
  theme(legend.position = 'none',
        axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold")) + 
  geom_point(data=krige_points_df, aes(x=x,y=y), size=3.5) + 
  ylab('y') + ylim(c(1, 5.5))

ggplot(krige_curve_df, aes(x=x)) +
  geom_line(aes(y=curve2, colour='red'), size=1.5) +
  theme(legend.position = 'none',
        axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold")) + 
  geom_point(data=krige_points_df, aes(x=x,y=y), size=3.5) + 
  ylab('y') + ylim(c(1, 5.5))

ggplot(krige_curve_df, aes(x=x)) +
  geom_line(data=data.frame(x=sim11[[1]][ ,1], y=sim11$data),
            aes(y=y), size=0.5, colour='black') +
  theme(legend.position = 'none',
        axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold")) + 
  geom_point(data=krige_points_df, aes(x=x,y=y), size=3.5) + 
  ylab('y') + ylim(c(1, 5.5))


plot(x=gold_coords, y=gold_pred_OC$simulations[, 1], type='l')
gold_sims <- as.data.frame(gold_pred_OC$simulations)
colnames(gold_sims) <- 1:10
gold_sims <- tidyr::gather(gold_sims)
gold_sims <- cbind(rep(gold_coords, 10), gold_sims)
gold_sims <- gold_sims[, c(1,3,2)]
colnames(gold_sims) <- c("x", "y", "sim_num")
gold_sims$sim_num <- as.numeric(gold_sims$sim_num)

# For the plots below, I used 1:1, 1:2, 1:5, and 1:10
ggplot(subset(gold_sims, sim_num %in% 1:10), aes(x=x, y=y, group=sim_num)) +
  geom_line(aes(alpha=0.001)) +
  theme(legend.position = "none") +
  xlab("Location") +
  ylab("Concentration") + ylim(0,6) +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"))



# Plots of several Matern correlation/covariance functions -------------------
x <- seq(0, 1, l = 101)
plot(x, cov.spatial(x, cov.model = "mat", kappa = 0.5, cov.pars=c(1,0.25)), type="l", xlab="u", 
     ylab=expression(rho(u)), ylim=c(0,1), main="Matern Correlation Functions")
lines(x, cov.spatial(x, cov.model = "mat", kappa = 1.5, cov.pars=c(1,0.16)), lty=2)
lines(x, cov.spatial(x, cov.model = "mat", kappa = 2.5, cov.pars=c(1,0.13)), lty=3)
legend(x=0.5, y=0.9, legend=c("kappa=0.5, phi=0.25","kappa=1.5, phi=0.16","kappa=2.5, phi=0.13"),
       lty=c(1,2,3))

matern_df <- data.frame(u = rep(x, 3),
           rho = c(cov.spatial(x, cov.model = "mat", kappa = 0.5, cov.pars=c(1,0.25)), 
                   cov.spatial(x, cov.model = "mat", kappa = 1.5, cov.pars=c(1,0.16)),
                   cov.spatial(x, cov.model = "mat", kappa = 2.5, cov.pars=c(1,0.13))),
           type = rep(c("0.5, 0.25","1.5, 0.16","2.5, 0.13"), each=length(x)))
ggplot(matern_df, aes(x=u, y=rho)) +
  geom_line(aes(colour=type)) +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"),
        legend.text=element_text(size=24),
        legend.title=element_text(size=24)) +
  labs(colour = expression(paste(kappa,", ",phi))) +
  ylab(expression(rho(u)))




# 2-D examples ----------------------------------------
set.seed(159)
image(grf(100^2, grid = "reg", cov.pars = c(1, 0.25)), col=gray(seq(1,0,l=51)), xlab="", ylab="")
set.seed(159)
image(grf(100^2, grid = "reg", cov.pars = c(1, 0.13)), col=gray(seq(1,0,l=51)), xlab="", ylab="")
image(grf(100^2, grid = "reg", cov.pars = c(1, 0.13), cov.model="matern", kappa=2.5), 
      col=gray(seq(1,0,l=51)), xlab="", ylab="")

data(SIC)
ml <- likfit(sic.all, ini = c(100, 40), nug = 10, lambda = 0.5, kappa=1)

gr <- pred_grid(sic.borders, by = 7.5)
KC <- krige.control(obj.model = ml)
OC <- output.control(n.pred = 1000, simul = TRUE, thres = 250)
set.seed(2419)
pred <- krige.conv(sic.all, loc = gr, borders = sic.borders, krige=KC, out=OC)

image(pred, col = gray(seq(1, 0.1, l = 21)), x.leg=c(0,350), y.leg=c(-60,-30))
image(pred, loc = gr, val = sqrt(pred$krige.var), col=gray(seq(1,0.1,l=21)), x.leg=c(0,350), y.leg=c(-60,-30))



