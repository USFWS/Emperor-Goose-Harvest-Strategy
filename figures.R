# code for NAAG figures
# Must run theta.logistic.r code to format data
library(tidyverse)
out <- out.eps
out
df <- data.frame(Year = 1985:2022, N = out$mean$N.est, upper = out$q97.5$N.est, 
                 lower = out$q2.5$N.est)
ggplot(data=ykd) + 
  geom_pointrange(aes(x=Year, y=itotal, ymin=itotal-2*itotal.se, 
                      ymax=itotal+2*itotal.se, col=Observer))  + 
  # geom_pointrange(data = har1, aes(x=Year+0.2, y=Harvest, ymin=lower, ymax=upper, 
  #                                  color=Season, shape = Type)) + 
  # geom_point(data = har2, aes(x=Year, y=PermitHar, color = Season, shape = Type)) + 
  geom_pointrange(data = df, aes(x=Year, y = N, ymin = lower, ymax = upper)) + 
  labs(y = "Indicated Total Birds/Harvest", title = "Model Eps") + 
  geom_pointrange(data = index, aes(x=Year+0.2, y = itotal, ymin = itotal - 2*itotal.se, 
                                    ymax = itotal + 2*itotal.se), color = "red")
#plot observer effects
df <- data.frame(Observer = levels(factor(ykd$Observer)), Effect = out$mean$alpha1, 
                 upper = out$q2.5$alpha1, lower = out$q97.5$alpha1) 
ggplot(data = df) + geom_pointrange(aes(x = Observer, y = Effect, ymin=lower, ymax = upper, 
                                    color = Observer)) + 
  geom_hline(aes(yintercept=0))

# #add closure threshold line
# ggplot(data=filter(ykd, Year < 2017)) + 
#   geom_pointrange(aes(x=Year, y=itotal, ymin=itotal-2*itotal.se, ymax=itotal+2*itotal.se)) + 
#   geom_hline(yintercept = 23000, size=2, color="red") + 
#   geom_pointrange(data = filter(har.all, Year < 2017), 
#                   aes(x=Year, y=Harvest, ymin=lower, ymax=upper, color=`Spring Season`)) + 
#   labs(x="Year", y="Number") + 
#   theme(legend.position="top", 
#         legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=15)) 
# ggsave("figures/2016-har-index-with-closure.png")
# 
# #to 2022
# ggplot(data=ykd) + 
#   geom_pointrange(aes(x=Year, y=itotal, ymin=itotal-2*itotal.se, ymax=itotal+2*itotal.se)) + 
#   #geom_hline(yintercept = 23000) + 
#   geom_pointrange(data = har.all, 
#                   aes(x=Year, y=Harvest, ymin=lower, ymax=upper, color=`Spring Season`)) + 
#   labs(x="Year", y="Number") + 
#   geom_segment(aes(x=2016.5, xend=2022, y=23000, yend=23000), color = "red", size=2)+
#   theme(legend.position="top", 
#         legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=15))
# ggsave("figures/2022-har-index.png")
# ## add the harvest mean posterior
# out <- readRDS("out1.harOriginal.RDS")
# df2 <- data.frame(x=c(2016.5, 2016.5, 2022, 2022), 
#                   y=c(out$q2.5$mu.green, out$q97.5$mu.green, out$q97.5$mu.green, out$q2.5$mu.green),
#                   y2=c(out$q2.5$mu.green-2*out$q97.5$sigma.har, out$q97.5$mu.green+2*out$q97.5$sigma.har, 
#                        out$q97.5$mu.green+2*out$q97.5$sigma.har, out$q2.5$mu.green-2*out$q97.5$sigma.har))
# ggplot(data=ykd) + 
#   geom_pointrange(aes(x=Year, y=itotal, ymin=itotal-2*itotal.se, ymax=itotal+2*itotal.se)) + 
#   #geom_hline(yintercept = 23000) + 
#   geom_pointrange(data = har.all, 
#                   aes(x=Year, y=Harvest, ymin=lower, ymax=upper, color=`Spring Season`)) + 
#   labs(x="Year", y="Number") + 
#   geom_segment(aes(x=2016.5, xend=2022, y=out$mean$mu.green, yend=out$mean$mu.green), 
#                color=3, size=2)+
#   geom_segment(aes(x=2016.5, xend=2022, y=23000, yend=23000), color = "red", size=2)+
#   geom_polygon(data=df2, aes(x=x, y=y), fill=3, alpha=0.5)+
#   geom_polygon(data=df2, aes(x=x, y=y2), fill=3, alpha=0.25)+
#   theme(legend.position="top", 
#         legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=15))
# ggsave("figures/2022-har-index-harpost1.png")
# out <- readRDS("out2.harOriginal.RDS")
# df2 <- data.frame(x=c(2016.5, 2016.5, 2022, 2022), 
#                   y=c(out$q2.5$mu.green, out$q97.5$mu.green, out$q97.5$mu.green, out$q2.5$mu.green),
#                   y2=c(out$q2.5$mu.green-2*out$q97.5$sigma.har, out$q97.5$mu.green+2*out$q97.5$sigma.har, 
#                        out$q97.5$mu.green+2*out$q97.5$sigma.har, out$q2.5$mu.green-2*out$q97.5$sigma.har))
# 
# ggplot(data=ykd) + 
#   geom_pointrange(aes(x=Year, y=itotal, ymin=itotal-2*itotal.se, ymax=itotal+2*itotal.se)) + 
#   #geom_hline(yintercept = 23000) + 
#   geom_pointrange(data = har.all, 
#                   aes(x=Year, y=Harvest, ymin=lower, ymax=upper, color=`Spring Season`)) + 
#   labs(x="Year", y="Number") + 
#   geom_segment(aes(x=2016.5, xend=2022, y=out$mean$mu.green, yend=out$mean$mu.green), 
#                color=3, size=2)+
#   geom_polygon(data=df2, aes(x=x, y=y), fill=3, alpha=0.5)+
#   geom_polygon(data=df2, aes(x=x, y=y2), fill=3, alpha=0.25)+
#   theme(legend.position="top", 
#         legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=15))
# ggsave("figures/2022-har-index-harpost2.png")
# ################################################################################
# #harvest priors
# x = seq(1, 20000, by=1)
# df = data.frame(Season = c(rep("Open", 10000), rep("Closed", 10000)), 
#                 Harvest = c(rlnorm(10000, log(11250), 0.2), rlnorm(10000, log(4081), 0.1)))
# ggplot(df, aes(x = Harvest, fill=Season, color= Season)) + 
#   geom_density(alpha=0.5) + 
#   theme(legend.position="top", 
#         legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=0))
# ggsave("figures/2016-har-priors.png")
# ################################################################################
# ## 2022 posterior results
# ## Demographic parameters
# out <- readRDS("out1.harOriginal.RDS")
# ##r_max
# x = rlnorm(10000, log(0.09), 0.3)
# x <- x[x < 0.2]
# df2 <- data.frame(r_max=x, Type = "Prior")
# df3 <- data.frame(r_max=out$sims.list$r.max, Type = "Posterior") 
# ggplot(df2, aes(x = r_max, fill=Type, color= Type)) + 
#   geom_density(alpha = 0.5) +
#   geom_density(data=df3, alpha=0.75) + 
#   theme(legend.position="top", 
#         legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=0))
# ggsave("figures/r_max_post.png")
# ## K
# df2 <- data.frame(K=runif(10000, 50000, 400000), Type = "Prior")
# df3 <- data.frame(K=out$sims.list$CC, Type = "Posterior") 
# ggplot(df2, aes(x = K, fill=Type, color= Type)) + 
#   geom_density(alpha = 0.5) +
#   geom_density(data=df3, alpha=0.75) + 
#   theme(legend.position="top", 
#         legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=0))
# ggsave("figures/K_post.png")
# ## Theta
# x=rlnorm(100000, log(1.76), 0.25)
# x <- x[x > 1 & x < 3]
# df2 <- data.frame(theta=x, Type = "Prior")
# df3 <- data.frame(theta=out$sims.list$theta, Type = "Posterior") 
# ggplot(df2, aes(x = theta, fill=Type, color= Type)) + 
#   geom_density(alpha = 0.5) +
#   geom_density(data=df3, alpha=0.75) + 
#   theme(legend.position="top", 
#         legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=0))
# ggsave("figures/theta_post.png")
# ## q
# x=runif(10000, 0.01, 1)
# df2 <- data.frame(q=x, Type = "Prior")
# df3 <- data.frame(q=out$sims.list$q, Type = "Posterior") 
# ggplot(df2, aes(x = q, fill=Type, color= Type)) + 
#   geom_density(alpha = 0.5) +
#   geom_density(data=df3, alpha=0.75) + 
#   theme(legend.position="top", 
#         legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=0))
# ggsave("figures/q_post.png")
# ## N 2022
# df3 <- data.frame(N=as.vector(out$sims.list$N.tot[,38]), Type = "Posterior") 
# ggplot(df3, aes(x = N, fill=Type, color= Type)) + 
#   geom_density(alpha = 0.5) +
#   theme(legend.position="top", 
#         legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=0))
# ggsave("figures/N2022_post.png")
# # process SD
# x=runif(10000, 0, 0.3)
# df2 <- data.frame(Sigma=x, Type = "Prior")
# df3 <- data.frame(Sigma=out$sims.list$sigma.proc, Type = "Posterior") 
# ggplot(df2, aes(x = Sigma, fill=Type, color= Type)) + 
#   geom_density(alpha = 0.5) +
#   geom_density(data=df3, alpha=0.75) + 
#   theme(legend.position="top", 
#         legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=0))
# ggsave("figures/Sigma_post.png")
# ## Harvest
# out <- readRDS("out1.harOriginal.RDS")
# df2 <- data.frame(Harvest=out$sims.list$mu.green[1:10000], Type = "Posterior", 
#                Season = "Open")
# #add 2022 posterior from original prior
# postHar <- ggplot(df, aes(x = Harvest, fill=Season, color= Season)) + 
#   geom_density(alpha=0.1) + 
#   geom_density(data=df2, alpha=0.5) + 
#   theme(legend.position="top", 
#         legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=0))
# print(postHar)
# ggsave("figures/postHar.1.png")
# #add posterior for out2
# out <- readRDS("out2.harOriginal.RDS")
# df3 <- data.frame(Harvest=out$sims.list$mu.green[1:10000], Type = "Posterior", 
#                   Season = "Open")
# 
# ggplot(df, aes(x = Harvest, fill=Season, color= Season)) + 
#   geom_density(alpha=0.1) + 
#   geom_density(data=df2, alpha=0.5) + 
#   geom_density(data=df3, alpha=1) + 
#   theme(legend.position="top", 
#         legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=0))
# ggsave("figures/postHar.2.png")
# ################################################################################
# ##Make population simulation figure
# out <- readRDS("out1.harOriginal.RDS")
# source("simulate.pop.R")
# Nsamples <- 100
# Tmax <- 100
# pick <- sample(1:length(out$sims.list$r.max), Nsamples)
# results <- matrix(NA, Tmax, Nsamples)
# for(i in 1:Nsamples){
#   results[,i] <- project.pop2(n1 = out$sims.list$N.tot[i,], 
#                               r = out$sims.list$r.max[i], 
#                               theta = out$sims.list$theta[i],
#                               K = out$sims.list$CC[i], 
#                               Hgreen = out$sims.list$mu.green[i],
#                               Hred = out$sims.list$m.har[i],
#                               sdH = out$sims.list$sigma.har[i],
#                               sdpop = out$sims.list$sigma.proc[i],
#                               q = out$sims.list$q[i],
#                               total = FALSE)$pop
# }
# df <- data.frame(Time = 1:Tmax, results) 
# gplot <- ggplot()
# for(i in 1:100){ #show just 100 samples
#   df2 <- data.frame(Time=df$Time, Pop=df[,i+1])
#   gplot <- gplot + geom_line(data=df2, aes(x=Time, y=Pop))
# }
# gplot <- gplot +
#   labs(x="Year", y="Total Birds") + 
#   geom_hline(yintercept = 23000, color="red", size=2) +
#   theme(legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=15))
# print(gplot)
# ggsave("figures/population.sims.png")
# ################################################################################
# ##Make MSY figures
# df <- readRDS("optim1.harOriginal.RDS")
# df$Clower <- df$cumHar - df$cumHarSD
# df$Cupper <- df$cumHar + df$cumHarSD
# df$AveHar = df$cumHar/200
# df$AveHarlower = df$Clower/200
# df$AveHarupper = df$Cupper/200
# ggplot(data = df, aes(x=Closure, y=AveHar)) + 
#   geom_point()+
#   geom_smooth(method="gam", se=TRUE) + 
#   labs(x="Closure", y="Cummulative Harvest") + 
#   theme(legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=15))
# ggsave("figures/cum-har-v-closure.png")
# #with 1 SD
# ggplot(data = df, aes(x=Closure, y=AveHar)) + 
#   geom_ribbon(aes(x=Closure, ymin=AveHarlower, ymax = AveHarupper), alpha=0.2, fill = "blue") +
#   geom_point()+
#   geom_smooth(method="gam", se=TRUE) + 
#   labs(x="Closure", y="Cummulative Harvest") + 
#   theme(legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=15))
# ggsave("figures/cum-har-v-closure-1SD.png")
# ##############################################
# ggplot(data = df, aes(x=Closure, y=pHunt)) + 
#   geom_point()+
#   geom_smooth(method="gam", se=TRUE)
# ##############################################
# ggplot(data = df, aes(x=Closure, y=mPop)) + 
#   geom_point()+
#   geom_smooth(method="gam", se=TRUE) + 
#   labs(x="Closure", y="Mean Population Size") + 
#   theme(legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=15))
# ggsave("figures/mPop-v-closure.png")
# #with SD
# df$Plower <- ifelse(df$mPop - df$mPopSD <0, 0, df$mPop - df$mPopSD)
# df$Pupper <- df$mPop + df$mPopSD
# ggplot(data = df, aes(x=Closure, y=mPop)) + 
#   geom_ribbon(aes(x=Closure, ymin=Plower, ymax = Pupper), alpha=0.2, fill = "blue") +
#   geom_point()+
#   geom_smooth(method="gam", se=TRUE) + 
#   labs(x="Closure", y="Mean Population Size") + 
#   theme(legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=15))
# ggsave("figures/mPop-v-closure-1SD.png")
# ##############################################
# #yield curve
# ggplot(data = df, aes(x=mPop, y=AveHar)) + 
#   geom_point()+
#   geom_smooth(method="gam", se=TRUE) +
#   labs(x="Population Size", y="Average Annual Harvest") + 
#   theme(legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=15))
# ggsave("figures/yield-curve.png")
# #with 1 SD
# ggplot(data = df, aes(x=mPop, y=AveHar)) + 
#   geom_ribbon(aes(x=mPop, ymin=AveHarlower, ymax = AveHarupper), alpha=0.2, fill = "blue") +
#   geom_point()+
#   geom_smooth(method="gam", se=TRUE) + 
#   labs(x="Population Size", y="Average Annual Harvest") + 
#   theme(legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=15))
# ggsave("figures/yield-curve-SD.png")
# ##############################################
# #find MSY closure level
# library(mgcv)
# fit <- gam(cumHar~s(Closure), data=df)
# pred <- predict(fit)
# df[which(pred == max(pred)),]
# #find % yield curve with min harvest
# slice_tail(df)$cumHar/df[which(pred == max(pred)),]$cumHar
# #94%!
# ################################################################################
# ## Utility curves
# ##Population
# x = seq(0, 1, by=0.001)
# u0 = x
# u1 = x^0.5
# u2 = x^2
# u3 = plogis(x, 0.5, 1)
# plot(x,u0, type="l")
# lines(x, u1)
# lines(x, u2)
# # from MBM
# # results from Elicitation of utility function from Eric Taylor
# # on 4-27-2016
# u1 <- list(N = c(0,50,60,80,100, 150, 175, 200), 
#            U = c(0,0, 0, 0.6, 0.75, 0.9, 1, 1))
# u2 <- list(N = c(0,60,80,100,125,155, 200), 
#            U = c(0,0.25,0.5,0.625, 0.75, 0.875, 1))
# 
# plot(u1$N, u1$U, pch=16, col=1, xlab="Population Size (x1000)", ylab="Utility", 
#      main = "Utility function from USFWS MBM (2016)", font.lab = 2)
# points(u2$N, u2$U, pch=16, col=2)
# legend("topleft", legend=c("First replicate", "Second replicate", "Best fit logistic function"), pch=c(16,16,NA), col=c(1,2,1), lty=c(NA, NA, 1))
# 
# u <- rbind(as.data.frame(u1),as.data.frame(u2))
# u$N2 <- u$N*u$N
# u$N3 <- u$N*u$N*u$N 
# 
# library(bbmle)
# nloglike <- function(m = 80, s=20, sigma=0.1){
#   -sum(dnorm(U, plogis(N, location=m, scale=s), sigma, log=TRUE))
# }
# 
# fit2 <- mle2(nloglike, start=list(m=80, s=20, sigma=0.1), data=u)
# summary(fit2)
# lines(1:200, plogis(1:200, location=coef(fit2)[1], scale=coef(fit2)[2]))
# #try ggplot
# df = data.frame(Utility = c(u1$U, u2$U), N = c(u1$N, u2$N), 
#                 Replicate = c(rep("First", length(u1$U)), rep("Second", length(u2$U))))
# df2 = data.frame(N = 1:200, Utility = plogis(1:200, location=coef(fit2)[1], scale=coef(fit2)[2]))
# ggplot(data = df, aes(x=N, y=Utility))+
#   geom_point(aes(color=Replicate), size=3)+
#   geom_path(data=df2, aes(x=N, y=Utility), size=1.5) + 
#   labs(x="Population Size", y="Utility") + 
#   theme(legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=15))
# ggsave("figures/utility.png")
# ##added linear utility
# ggsave("figures/utility.png")
# ggplot(data = df, aes(x=N, y=Utility))+
#   geom_point(aes(color=Replicate), size=3)+
#   #geom_path(data=df2, aes(x=N, y=Utility), size=1.5) + 
#   geom_segment(aes(x=0, y=0, xend=200, yend=1), size=1.5) +
#   #geom_smooth(method="lm") +
#   labs(x="Population Size", y="Utility") + 
#   theme(legend.text = element_text(size=20), 
#         legend.title = element_text(size=20),
#         axis.title = element_text(face="bold", size = 22), 
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=15))
# ggsave("figures/utility2.png")
#   
# #   ##ADFG utility:
# #   ##Fit curve to utility values elicited from ADFG representatives
# #   ## Mostly Dan Rosenberg, but Steve Fleischmen, Jason Schamber..., Bob Stehn, and Mike Petrula?  were present
# #   ##  Steve Bob and Jason gave feedback on utilities
# #   
# #   ##  Elicitation occurred on 5-20-2016
# #   # on 4-27-2016
# #   u1 <- list(N = c(0, 50, 75, 90, 150, 175, 200), 
# #              U = c(0, 0, 0.33, 0.8, 0.95, 1, 1))
# # u2 <- list(N = c(0, 65, 82, 112, 155, 180, 200), 
# #            U = c(0, 0.125, 0.25, 0.5, 0.75, 0.875, 1))
# # #pdf("plot1.pdf")
# # plot(u1$N, u1$U, pch=21, bg=1, xlab="Population Size (x1000)", ylab="Utility")
# # points(u2$N, u2$U, pch=22, bg=1)
# # u <- rbind(as.data.frame(u1),as.data.frame(u2))
# # u$N2 <- u$N*u$N
# # u$N3 <- u$N*u$N*u$N 
# # library(bbmle)
# # nloglike <- function(m = 80, s=20, sigma=0.1){
# #   -sum(dnorm(U, plogis(N, location=m, scale=s), sigma, log=TRUE))
# # }
# # fit.adfg <- mle2(nloglike, start=list(m=80, s=20, sigma=0.1), data=u)
# # lines(1:200, plogis(1:200, location=coef(fit.adfg)[1], scale=coef(fit.adfg)[2]), lwd=2, col=1)
# # # results from Elicitation of utiity function from Eric Tayor
# # # on 4-27-2016
# # u1 <- list(N = c(0,50,60,80,100, 150, 175, 200), 
# #            U = c(0,0, 0, 0.6, 0.75, 0.9, 1, 1))
# # u2 <- list(N = c(0,60,80,100,125,155, 200), 
# #            U = c(0,0.25,0.5,0.625, 0.75, 0.875, 1))
# # points(u1$N, u1$U, pch=21, bg="gray50")
# # points(u2$N, u2$U, pch=22, bg="gray50")
# # u <- rbind(as.data.frame(u1),as.data.frame(u2))
# # u$N2 <- u$N*u$N
# # u$N3 <- u$N*u$N*u$N 
# # nloglike <- function(m = 80, s=20, sigma=0.1){
# #   -sum(dnorm(U, plogis(N, location=m, scale=s), sigma, log=TRUE))
# # }
# # fit.usfws <- mle2(nloglike, start=list(m=80, s=20, sigma=0.1), data=u)
# # 
# # lines(1:200, plogis(1:200, location=coef(fit.usfws)[1], scale=coef(fit.usfws)[2]), col="gray50", lwd=2)
# # 
# # 
# # legend("topleft", legend=c("First replicate", "Second replicate", "Best fit logistic function", "ADFG", "USFWS"), 
# #        pch=c(21,22,NA, 22, 22), col=c(1,1,1,1,"gray50"), pt.bg=c(NA,NA,NA, 1, "gray50"), lty=c(NA, NA, 1, 1,1), lwd=c(NA, NA, 2,2,2))
# # 
# 
# 
