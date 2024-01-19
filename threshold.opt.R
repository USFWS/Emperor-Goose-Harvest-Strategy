# function to iterate through closure thresholds and find cumulative harvest and 
#  population size or utility
library(tidyverse)
#load data
out <- readRDS("data/out.RDS")
#load linear observation model
fit <- readRDS("data/fit.RDS")
#load simulate function
source("simulate.pop.R")
#also need to run head of "theta.logistic.R to get fit object
#define threshold to search across; this is observed YKD index values to search over
t_red <- seq(0, 50000, by= 250) #250
Nsamples <- 1000 #number of samples of the posterior 1000
TimeHor <- 100 #time over which to calculate return/yield; time horizon 100
discount <- 1 # 1 - time discount rate
df <- data.frame(Closure=t_red,
                 cumHar=rep(NA, length(t_red)),
                 pHunt=rep(NA, length(t_red)),
                 mPop = rep(NA, length(t_red)),
                 sdhar = rep(NA, length(t_red)),
                 GreenHar = rep(NA, length(t_red)),
                 pChange = rep(NA, length(t_red)),
                 cumHarSD=rep(NA, length(t_red)),
                 pHuntSD=rep(NA, length(t_red)),
                 mPopSD = rep(NA, length(t_red)))

#set up loops
temp <- list(cumhar = numeric(), phunt = numeric(), mPop = numeric(), 
             sdhar = numeric(), GreenHar = numeric(), pChange = numeric(), 
             pExt = numeric())
pick <- sample(1:length(out$sims.list$r.max), Nsamples)
for(t in 1:length(t_red)){
  for(i in 1:Nsamples){
    reward <- project_pop(Tmax = TimeHor,
                           n1 = out$sims.list$N.tot[i,],
                           r = out$sims.list$r.max[i],
                           theta = out$sims.list$theta[i],
                           K = out$sims.list$CC[i],
                           threshold = list(red = t_red[t], green = t_red[t]), 
                           Hlist = list(red = out$sims.list$m.har[i],
                                        yellow = out$sims.list$m.har[i],
                                        green = out$sims.list$mu.green[i]),
                           sdH = list(red = out$sims.list$sigma.har[i],
                                      yellow = out$sims.list$sigma.har[i],
                                      green = out$sims.list$sigma.har[i]),
                           Hp = list(green = out$sims.list$m.har.p[i]),
                           sdHp = list(green = out$sims.list$sd.har.p[i]),
                           sdpop = out$sims.list$sigma.proc[i],
                           q = out$sims.list$q[i],
                           alpha1 = c(out$sims.list$alpha1[i, 5], 
                                      out$sims.list$alpha1[i, 10]),
                           sdesp = out$sims.list$sd.esp[i],
                           BETA = fit$coefficients[1], #linear model slope
                           SE = fit$coefficients[2], #linear model slope SE
                           SIGMA = fit$sigma, #linear model RMSE
                           DF = fit$df[2], #linear model degrees of freedom
                           crip = out$sims.list$c[i])
    
    temp$cumhar[i] <- sum(reward$har, na.rm = TRUE)
    temp$sdhar[i] <- sd(reward$har, na.rm = TRUE)
    temp$GreenHar[i] <- sum(reward$har[reward$hunt == 3], na.rm = TRUE)
    temp$pChange[i] <- sum( (reward$hunt[-1] != reward$hunt[-TimeHor]) / 
                                length(reward$hunt[-1]))
    temp$phunt[i] <- sum(reward$hunt %in% c(2,3), na.rm = TRUE)/length(reward$hunt)
    temp$mPop[i] <- mean(reward$pop, na.rm = TRUE)
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
  df$pChange[t] <- mean(temp$pChange)
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
ggplot(data = df, aes(x=Closure, y=pChange)) + 
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
fit2 <- gam(cumHar~s(Closure), data=df)
pred <- predict(fit2)
df[which(pred == max(pred)),]
#find % yield curve with min harvest
slice_tail(df)$cumHar/df[which(pred == max(pred)),]$cumHar
