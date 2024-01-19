# function to iterate through closure thresholds and find cumulative harvest and 
#  population size or utility
# Add a "Yellow" season and explore optimum.
# Here the "Yellow" season is defined by a mean harvest mid way between "red" and "green" 
library(tidyverse)
#load data
out <- readRDS("data/out.RDS")
#load linear observation model
fit <- readRDS("data/fit.RDS")
#load simulate function
source("simulate.pop.R")
#also need to run head of "theta.logistic.R to get fit object
#define threshold to search across; this is observed YKD index values to search over
tstep <- 250
t_red <- seq(0, 50000, by= tstep) #250
Nsamples <- 1000 #number of samples of the posterior 1000
TimeHor <- 100 #time over which to calculate return/yield; time horizon 100
#discount <- 1 # 1 - time discount rate, not implemented
df <- list(Closure = numeric(),
                 Green = numeric(),
                 cumHar = numeric(),
                 pHunt = numeric(),
                 mPop = numeric(),
                 sdhar = numeric(),
                 GreenHar = numeric(),
                 pChange = numeric(),
                 cumHarSD = numeric(),
                 pHuntSD = numeric(),
                 mPopSD = numeric())

#set up loops
temp <- list(cumhar = numeric(), phunt = numeric(), mPop = numeric(), 
             sdhar = numeric(), GreenHar = numeric(), pChange = numeric(), 
             pExt = numeric())
pick <- sample(1:length(out$sims.list$r.max), Nsamples)
tt <- 1
for(t1 in 1:length(t_red)){
  t_green <- seq(t_red[t1], 50000, by = tstep)
  for(t2 in 1:length(t_green)){
    ii <- 1
    for(i in pick){
      reward <- project_pop(Tmax = TimeHor,
                            n1 = out$sims.list$N.tot[i,],
                            r = out$sims.list$r.max[i],
                            theta = out$sims.list$theta[i],
                            K = out$sims.list$CC[i],
                            threshold = list(red = t_red[t1], green = t_red[t2]), 
                            Hlist = list(red = out$sims.list$m.har[i],
                                         yellow = out$sims.list$m.har[i] + 
                                           (out$sims.list$mu.green[i] - 
                                              out$sims.list$m.har[i])/2,
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
      
      temp$cumhar[ii] <- sum(reward$har, na.rm = TRUE)
      temp$sdhar[ii] <- sd(reward$har, na.rm = TRUE)
      temp$GreenHar[ii] <- sum(reward$har[reward$hunt == 3], na.rm = TRUE)
      temp$pChange[ii] <- sum( (reward$hunt[-1] != reward$hunt[-TimeHor]) / 
                                 length(reward$hunt[-1]))
      temp$phunt[ii] <- sum(reward$hunt %in% c(2,3), na.rm = TRUE)/length(reward$hunt)
      temp$mPop[ii] <- mean(reward$pop, na.rm = TRUE)
      temp$pExt[ii] <- ifelse(last(reward$pop) == 0, 1, 0)
      ii <- ii + 1
    }
  #calculate summary/utility over parameter samples
  #Use mean as surrogate utility function for now
  #mean across parametric uncertain and time (stochastic variation)
  df$Closure[tt] <- t_red[t1]
  df$Green[tt] <- t_green[t2]
  df$cumHar[tt] <- mean(temp$cumhar)/TimeHor 
  df$pHunt[tt] <- mean(temp$phunt)
  df$mPop[tt] <- mean(temp$mPop)
  df$sdhar[tt] <- mean(temp$sdhar) 
  df$GreenHar[tt] <- mean(temp$GreenHar)/TimeHor
  df$pChange[tt] <- mean(temp$pChange)
  df$pExt[tt] <- mean(temp$pExt)
  #sd across parametric and time stochastic variation)
  df$cumHarSD[tt] <- sd(temp$cumhar/TimeHor) 
  df$pHuntSD[tt] <- sd(temp$phunt)
  df$mPopSD[tt] <- sd(temp$mPop)
  tt <- tt + 1
}}
df <- as.data.frame(df)
saveRDS(df, file = "data/optim2.RYG.RDS")
# df <- readRDS("data/optim1.RDS")
df$Clower <- df$cumHar - df$cumHarSD
df$Cupper <- df$cumHar + df$cumHarSD
dfmax <- filter(df, cumHar == max(cumHar))
library(viridis) 
ggplot() +  
  geom_tile(data=df, aes(x=Closure, y=Green, fill=cumHar)) + 
  scale_fill_viridis() + 
  geom_point(data = dfmax, aes(x = Closure, y = Green))
ggplot() +  
  geom_tile(data=df, aes(x=Closure, y=Green, fill=pHunt)) + 
  scale_fill_viridis()
ggplot() +  
  geom_tile(data=df, aes(x=Closure, y=Green, fill=mPop)) + 
  scale_fill_viridis()
ggplot() +  
  geom_tile(data=df, aes(x=Closure, y=Green, fill=pChange)) + 
  scale_fill_viridis()
ggplot() +  
  geom_tile(data=df, aes(x=Closure, y=Green, fill=pExt)) + 
  scale_fill_viridis()

#find MSY closure level
library(mgcv)
fit2 <- gam(cumHar~s(Closure, Green), data=df)
pred <- predict(fit2)
df[which(pred == max(pred)),]
#find % yield curve with min harvest
slice_tail(df)$cumHar/df[which(pred == max(pred)),]$cumHar
