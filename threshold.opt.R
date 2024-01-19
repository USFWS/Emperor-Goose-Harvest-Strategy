# function to iterate through closure thresholds and find cumulative harvest and population size or utility
library(tidyverse)
#load data
out <- readRDS("data/out.RDS")
#load simulate function
source("simulate.pop.R")
#also need to run head of "theta.logistic.R to get fit object
#define threshold to search across; this is observed YKD index values to search over
t_red <- seq(0, 50000, by= 250)
Nsamples <- 1000 #number of samples of the posterior
TimeHor <- 100 #time over which to calculate return/yield; time horizon
discount <- 1 #time discount rate
df <- data.frame(Closure=t_red,
                 cumHar=rep(NA, length(t_red)),
                 pHunt=rep(NA, length(t_red)),
                 mPop = rep(NA, length(t_red)),
                 sdhar = rep(NA, length(t_red)),
                 GreenHar = rep(NA, length(t_red)),
                 pOpenClose = rep(NA, length(t_red)),
                 cumHarSD=rep(NA, length(t_red)),
                 pHuntSD=rep(NA, length(t_red)),
                 mPopSD = rep(NA, length(t_red)))

#set up loops
temp <- list(cumhar = numeric(), phunt = numeric(), mPop = numeric(), 
             sdhar = numeric(), GreenHar = numeric(), pOpenClose = numeric(), 
             pExt = numeric())
pick <- sample(1:length(out$sims.list$r.max), Nsamples)
for(t in 1:length(t_red)){
  for(i in 1:Nsamples){
    reward <- project.pop2(Tmax = TimeHor+1,
                           n1 = out$sims.list$N.tot[i,],
                           r = out$sims.list$r.max[i],
                           theta = out$sims.list$theta[i],
                           K = out$sims.list$CC[i],
                           Hgreen = out$sims.list$mu.green[i],
                           Hred = out$sims.list$m.har[i],
                           sdH = out$sims.list$sigma.har[i],
                           Hp = out$sims.list$m.har.p,
                           sdHp = out$sims.list$sd.har.p,
                           sdpop = out$sims.list$sigma.proc[i],
                           q = out$sims.list$q[i],
                           alpha1 = c(out$sims.list$alpha1[i, 5], 
                                      out$sims.list$alpha1[i, 10]),
                           sdesp = out$sims.list$sd.esp[i],
                           BETA = fit$coefficients[1], #linear model slope
                           SE = fit$coefficients[2], #linear model slope SE
                           SIGMA = fit$sigma, #linear model RMSE
                           DF = fit$df[2], #linear model degrees of freedom
                           crip = out$sims.list$c[i],
                           t_close = t_red[t],
                           total = TRUE)
    temp$cumhar[i] <- sum(reward$har[-(TimeHor+1)], na.rm = TRUE)
    temp$sdhar[i] <- sd(reward$har[-(TimeHor+1)], na.rm = TRUE)
    temp$GreenHar[i] <- sum(reward$har[-(TimeHor+1)]*reward$hunt[-(TimeHor+1)], na.rm = TRUE)
    temp$pOpenClose[i] <- sum(reward$hunt[-c(1,(TimeHor+1))] != reward$hunt[-c(TimeHor, (TimeHor+1))])/(length(reward$hunt[-(TimeHor+1)])-1)
    temp$phunt[i] <- sum(reward$hunt[-(TimeHor+1)], na.rm = TRUE)/length(reward$hunt[-(TimeHor+1)])
    temp$mPop[i] <- mean(reward$pop[-1], na.rm = TRUE)
    temp$pExt[i] <- ifelse(last(reward$pop) == 0, 1, 0)
  }
  #calculate summary/utility over parameter samples
  #Use mean as surrogate utility function for now
  #mean across parametric uncertain and time (stochastic variation)
  df$cumHar[t] <- mean(temp$cumhar)/TimeHor 
  df$pHunt[t] <- mean(temp$phunt)
  df$mPop[t] <- mean(temp$mPop)
  df$sdhar[t] <- mean(temp$sdhar) 
  df$GreenHar[t] <- mean(temp$GreenHar)/TimeHor
  df$pOpenClose[t] <- mean(temp$pOpenClose)
  df$pExt[t] <- mean(temp$pExt)
  #sd across parametric and time stochastic variation)
  df$cumHarSD[t] <- sd(temp$cumhar/TimeHor) 
  df$pHuntSD[t] <- sd(temp$phunt)
  df$mPopSD[t] <- sd(temp$mPop)

}

saveRDS(df, file = "data/optim1.RDS")
# df <- readRDS("data/optim1.RDS")
df$Clower <- df$cumHar - df$cumHarSD
df$Cupper <- df$cumHar + df$cumHarSD
ggplot(data = df, aes(x=Closure, y=cumHar)) + 
  geom_ribbon(aes(x=Closure, ymin=Clower, ymax = Cupper), alpha=0.2, fill = "blue") +
  geom_point()+
  geom_smooth(method="gam", se=TRUE)
ggplot(data = df, aes(x=Closure, y=pHunt)) + 
  geom_point()+
  geom_smooth(method="gam", se=TRUE)
ggplot(data = df, aes(x=Closure, y=mPop)) + 
  geom_point()+
  geom_smooth(method="gam", se=TRUE)
ggplot(data = df, aes(x=Closure, y=GreenHar)) + 
  #geom_ribbon(aes(x=Closure, ymin=Clower, ymax = Cupper), alpha=0.2, fill = "blue") +
  geom_point()+
  geom_smooth(method="gam", se=TRUE)
ggplot(data = df, aes(x=Closure, y=pOpenClose)) + 
  #geom_ribbon(aes(x=Closure, ymin=Clower, ymax = Cupper), alpha=0.2, fill = "blue") +
  geom_point()+
  geom_smooth(method="gam", se=TRUE)
ggplot(data = df, aes(x=Closure, y=pExt)) + 
  #geom_ribbon(aes(x=Closure, ymin=Clower, ymax = Cupper), alpha=0.2, fill = "blue") +
  geom_point()+
  geom_smooth(method="gam", se=TRUE)
ggplot(data = df, aes(x=Closure, y=cumHarSD)) + 
  #geom_ribbon(aes(x=Closure, ymin=Clower, ymax = Cupper), alpha=0.2, fill = "blue") +
  geom_point()+
  geom_smooth(method="gam", se=TRUE)
ggplot(data = df, aes(x=Closure, y=pHuntSD)) + 
  #geom_ribbon(aes(x=Closure, ymin=Clower, ymax = Cupper), alpha=0.2, fill = "blue") +
  geom_point()+
  geom_smooth(method="gam", se=TRUE)
ggplot(data = df, aes(x=Closure, y=mPopSD)) + 
  #geom_ribbon(aes(x=Closure, ymin=Clower, ymax = Cupper), alpha=0.2, fill = "blue") +
  geom_point()+
  geom_smooth(method="gam", se=TRUE)
#yield curve
ggplot(data = df, aes(x=mPop, y=cumHar)) + 
  geom_point()+
  geom_smooth(method="gam", se=TRUE)

#find MSY closure level
library(mgcv)
fit <- gam(cumHar~s(Closure), data=df)
pred <- predict(fit)
df[which(pred == max(pred)),]
#find % yield curve with min harvest
slice_tail(df)$cumHar/df[which(pred == max(pred)),]$cumHar

################################################################################
# out <- readRDS("out2.harOriginal.RDS")
# #define threshold to search across; this is observed YKD index values to search over
# t_red <- seq(0, 50000, by= 250)
# Nsamples <- 10000 #number of samples of the posterior
# TimeHor <- 200 #time over which to calculate return/yield; time horizon
# discount <- 1 #time discount rate
# df <- data.frame(Closure=t_red,
#                  cumHar=rep(NA, length(t_red)),
#                  pHunt=rep(NA, length(t_red)),
#                  mPop = rep(NA, length(t_red)),
#                  cumHarSD=rep(NA, length(t_red)),
#                  pHuntSD=rep(NA, length(t_red)),
#                  mPopSD = rep(NA, length(t_red)))
# 
# #set up loops
# temp <- list(cumhar = numeric(), phunt = numeric(), mPop = numeric())
# pick <- sample(1:length(out$sims.list$r.max), Nsamples)
# for(t in 1:length(t_red)){
#   for(i in 1:Nsamples){
#     reward <- project.pop2(Tmax = TimeHor+1,
#                            n1 = out$sims.list$N.tot[pick[i],], 
#                            r = out$sims.list$r.max[pick[i]], 
#                            theta = out$sims.list$theta[pick[i]],
#                            K = out$sims.list$CC[pick[i]], 
#                            Hgreen = out$sims.list$mu.green[pick[i]],
#                            Hred = out$sims.list$m.har[pick[i]],
#                            sdH = out$sims.list$sigma.har[pick[i]],
#                            sdpop = out$sims.list$sigma.proc[pick[i]],
#                            q = out$sims.list$q[pick[i]],
#                            t_close = t_red[t],
#                            total = TRUE)
#     temp$cumhar[i] <- sum(reward$har, na.rm = TRUE)
#     temp$phunt[i] <- sum(reward$hunt, na.rm = TRUE)/length(reward$hunt[!is.na(reward$hunt)])
#     temp$mPop[i] <- mean(reward$pop, na.rm = TRUE)
#   }
#   df$cumHar[t] <- mean(temp$cumhar) #Use mean as surrogate utility function for now
#   df$pHunt[t] <- mean(temp$phunt)
#   df$mPop[t] <- mean(temp$mPop)
#   df$cumHarSD[t] <- sd(temp$cumhar)
#   df$pHuntSD[t] <- sd(temp$phunt)
#   df$mPopSD[t] <- sd(temp$mPop)
# }
# 
# saveRDS(df, file = "optim2.harOriginal.RDS")
# # df <- readRDS("optim2.RDS")
# df$Clower <- df$cumHar - df$cumHarSD
# df$Cupper <- df$cumHar + df$cumHarSD
# ggplot(data = df, aes(x=Closure, y=cumHar)) + 
#   #geom_ribbon(aes(x=Closure, ymin=Clower, ymax = Cupper), alpha=0.2, fill = "blue") +
#   geom_point()+
#   geom_smooth(method="gam", se=TRUE)
# ggplot(data = df, aes(x=Closure, y=pHunt)) + 
#   geom_point()+
#   geom_smooth(method="gam", se=TRUE)
# ggplot(data = df, aes(x=Closure, y=mPop)) + 
#   geom_point()+
#   geom_smooth(method="gam", se=TRUE)
# #yield curve
# ggplot(data = df, aes(x=mPop, y=cumHar)) + 
#   geom_point()+
#   geom_smooth(method="gam", se=TRUE)
# #find MSY closure level
# library(mgcv)
# fit <- gam(cumHar~s(Closure), data=df)
# pred <- predict(fit)
# df[which(pred == max(pred)),]
# #find % yield curve with min harvest
# slice_tail(df)$cumHar/df[which(pred == max(pred)),]$cumHar
# #94%!
# ################################################################################
# ## Now for the 2016 posterior
# #load data
# out <- readRDS("out0.harOriginal.RDS")
# #define threshold to search across; this is observed YKD index values to search over
# t_red <- seq(0, 50000, by= 250)
# Nsamples <- 10000 #number of samples of the posterior
# TimeHor <- 200 #time over which to calculate return/yield; time horizon
# discount <- 1 #time discount rate
# df <- data.frame(Closure=t_red,
#                  cumHar=rep(NA, length(t_red)),
#                  pHunt=rep(NA, length(t_red)),
#                  mPop = rep(NA, length(t_red)),
#                  cumHarSD=rep(NA, length(t_red)),
#                  pHuntSD=rep(NA, length(t_red)),
#                  mPopSD = rep(NA, length(t_red)))
# 
# #set up loops
# temp <- list(cumhar = numeric(), phunt = numeric(), mPop = numeric())
# pick <- sample(1:length(out$sims.list$r.max), Nsamples)
# for(t in 1:length(t_red)){
#   for(i in 1:Nsamples){
#     reward <- project.pop2(Tmax = TimeHor+1,
#                            n1 = out$sims.list$N.tot[pick[i],], 
#                            r = out$sims.list$r.max[pick[i]], 
#                            theta = out$sims.list$theta[pick[i]],
#                            K = out$sims.list$CC[pick[i]], 
#                            Hgreen = out$sims.list$mu.green[pick[i]],
#                            Hred = out$sims.list$m.har[pick[i]],
#                            sdH = out$sims.list$sigma.har[pick[i]],
#                            sdpop = out$sims.list$sigma.proc[pick[i]],
#                            q = out$sims.list$q[pick[i]],
#                            t_close = t_red[t],
#                            total = TRUE)
#     temp$cumhar[i] <- sum(reward$har, na.rm = TRUE)
#     temp$phunt[i] <- sum(reward$hunt, na.rm = TRUE)/length(reward$hunt[!is.na(reward$hunt)])
#     temp$mPop[i] <- mean(reward$pop, na.rm = TRUE)
#   }
#   df$cumHar[t] <- mean(temp$cumhar) #Use mean as surrogate utility function for now
#   df$pHunt[t] <- mean(temp$phunt)
#   df$mPop[t] <- mean(temp$mPop)
#   df$cumHarSD[t] <- sd(temp$cumhar)
#   df$pHuntSD[t] <- sd(temp$phunt)
#   df$mPopSD[t] <- sd(temp$mPop)
# }
# 
# saveRDS(df, file = "optim0.harOriginal.RDS")
# df <- readRDS("optim0.harOriginal.RDS")
# df$Clower <- df$cumHar - df$cumHarSD
# df$Cupper <- df$cumHar + df$cumHarSD
# ggplot(data = df, aes(x=Closure, y=cumHar)) + 
#   geom_ribbon(aes(x=Closure, ymin=Clower, ymax = Cupper), alpha=0.2, fill = "blue") +
#   geom_point()+
#   geom_smooth(method="gam", se=TRUE)
# ggplot(data = df, aes(x=Closure, y=pHunt)) + 
#   geom_point()+
#   geom_smooth(method="gam", se=TRUE)
# ggplot(data = df, aes(x=Closure, y=mPop)) + 
#   geom_point()+
#   geom_smooth(method="gam", se=TRUE)
# #yield curve
# ggplot(data = df, aes(x=mPop, y=cumHar)) + 
#   geom_point()+
#   geom_smooth(method="gam", se=TRUE)
# #find closure level
# fit <- gam(cumHar~s(Closure), data=df)
# pred <- predict(fit)
# df[which(pred == max(pred)),]
