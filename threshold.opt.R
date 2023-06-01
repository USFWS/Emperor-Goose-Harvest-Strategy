# function to iterate through closure thresholds and find cummulative harvest and population size or utility
#load data
# start with out1 posterior: posterior from 2016 original prior
out <- readRDS("out1.harOriginal.RDS")
#define threshold to search across; this is observed YKD index values to search over
t_red <- seq(0, 50000, by= 250)
Nsamples <- 10000 #number of samples of the posterior
TimeHor <- 200 #time over which to calculate return/yield; time horizon
discount <- 1 #time discount rate
df <- data.frame(Closure=t_red,
                 cumHar=rep(NA, length(t_red)),
                 pHunt=rep(NA, length(t_red)),
                 mPop = rep(NA, length(t_red)),
                 cumHarSD=rep(NA, length(t_red)),
                 pHuntSD=rep(NA, length(t_red)),
                 mPopSD = rep(NA, length(t_red)))

#set up loops
temp <- list(cumhar = numeric(), phunt = numeric(), mPop = numeric())
pick <- sample(1:length(out$sims.list$r.max), Nsamples)
for(t in 1:length(t_red)){
  for(i in 1:Nsamples){
    reward <- project.pop2(Tmax = TimeHor+1,
                           n1 = out$sims.list$N.tot[pick[i],], 
                           r = out$sims.list$r.max[pick[i]], 
                           theta = out$sims.list$theta[pick[i]],
                           K = out$sims.list$CC[pick[i]], 
                           Hgreen = out$sims.list$mu.green[pick[i]],
                           Hred = out$sims.list$m.har[pick[i]],
                           sdH = out$sims.list$sigma.har[pick[i]],
                           sdpop = out$sims.list$sigma.proc[pick[i]],
                           q = out$sims.list$q[pick[i]],
                           t_close = t_red[t],
                           total = TRUE)
    temp$cumhar[i] <- sum(reward$har, na.rm = TRUE)
    temp$phunt[i] <- sum(reward$hunt, na.rm = TRUE)/length(reward$hunt[!is.na(reward$hunt)])
    temp$mPop[i] <- mean(reward$pop, na.rm = TRUE)
  }
  df$cumHar[t] <- mean(temp$cumhar) #Use mean as surrogate utility function for now
  df$pHunt[t] <- mean(temp$phunt)
  df$mPop[t] <- mean(temp$mPop)
  df$cumHarSD[t] <- sd(temp$cumhar)
  df$pHuntSD[t] <- sd(temp$phunt)
  df$mPopSD[t] <- sd(temp$mPop)
}

saveRDS(df, file = "optim1.harOriginal.RDS")
# df <- readRDS("optim1..harOriginal.RDS")
df$Clower <- df$cumHar - df$cumHarSD
df$Cupper <- df$cumHar + df$cumHarSD
ggplot(data = df, aes(x=Closure, y=cumHar)) + 
  #geom_ribbon(aes(x=Closure, ymin=Clower, ymax = Cupper), alpha=0.2, fill = "blue") +
  geom_point()+
  geom_smooth(method="gam", se=TRUE)
ggplot(data = df, aes(x=Closure, y=pHunt)) + 
  geom_point()+
  geom_smooth(method="gam", se=TRUE)
ggplot(data = df, aes(x=Closure, y=mPop)) + 
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
#94%!
################################################################################
out <- readRDS("out2.harOriginal.RDS")
#define threshold to search across; this is observed YKD index values to search over
t_red <- seq(0, 50000, by= 250)
Nsamples <- 10000 #number of samples of the posterior
TimeHor <- 200 #time over which to calculate return/yield; time horizon
discount <- 1 #time discount rate
df <- data.frame(Closure=t_red,
                 cumHar=rep(NA, length(t_red)),
                 pHunt=rep(NA, length(t_red)),
                 mPop = rep(NA, length(t_red)),
                 cumHarSD=rep(NA, length(t_red)),
                 pHuntSD=rep(NA, length(t_red)),
                 mPopSD = rep(NA, length(t_red)))

#set up loops
temp <- list(cumhar = numeric(), phunt = numeric(), mPop = numeric())
pick <- sample(1:length(out$sims.list$r.max), Nsamples)
for(t in 1:length(t_red)){
  for(i in 1:Nsamples){
    reward <- project.pop2(Tmax = TimeHor+1,
                           n1 = out$sims.list$N.tot[pick[i],], 
                           r = out$sims.list$r.max[pick[i]], 
                           theta = out$sims.list$theta[pick[i]],
                           K = out$sims.list$CC[pick[i]], 
                           Hgreen = out$sims.list$mu.green[pick[i]],
                           Hred = out$sims.list$m.har[pick[i]],
                           sdH = out$sims.list$sigma.har[pick[i]],
                           sdpop = out$sims.list$sigma.proc[pick[i]],
                           q = out$sims.list$q[pick[i]],
                           t_close = t_red[t],
                           total = TRUE)
    temp$cumhar[i] <- sum(reward$har, na.rm = TRUE)
    temp$phunt[i] <- sum(reward$hunt, na.rm = TRUE)/length(reward$hunt[!is.na(reward$hunt)])
    temp$mPop[i] <- mean(reward$pop, na.rm = TRUE)
  }
  df$cumHar[t] <- mean(temp$cumhar) #Use mean as surrogate utility function for now
  df$pHunt[t] <- mean(temp$phunt)
  df$mPop[t] <- mean(temp$mPop)
  df$cumHarSD[t] <- sd(temp$cumhar)
  df$pHuntSD[t] <- sd(temp$phunt)
  df$mPopSD[t] <- sd(temp$mPop)
}

saveRDS(df, file = "optim2.harOriginal.RDS")
# df <- readRDS("optim2.RDS")
df$Clower <- df$cumHar - df$cumHarSD
df$Cupper <- df$cumHar + df$cumHarSD
ggplot(data = df, aes(x=Closure, y=cumHar)) + 
  #geom_ribbon(aes(x=Closure, ymin=Clower, ymax = Cupper), alpha=0.2, fill = "blue") +
  geom_point()+
  geom_smooth(method="gam", se=TRUE)
ggplot(data = df, aes(x=Closure, y=pHunt)) + 
  geom_point()+
  geom_smooth(method="gam", se=TRUE)
ggplot(data = df, aes(x=Closure, y=mPop)) + 
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
#94%!
################################################################################
## Now for the 2016 posterior
#load data
out <- readRDS("out0.harOriginal.RDS")
#define threshold to search across; this is observed YKD index values to search over
t_red <- seq(0, 50000, by= 250)
Nsamples <- 10000 #number of samples of the posterior
TimeHor <- 200 #time over which to calculate return/yield; time horizon
discount <- 1 #time discount rate
df <- data.frame(Closure=t_red,
                 cumHar=rep(NA, length(t_red)),
                 pHunt=rep(NA, length(t_red)),
                 mPop = rep(NA, length(t_red)),
                 cumHarSD=rep(NA, length(t_red)),
                 pHuntSD=rep(NA, length(t_red)),
                 mPopSD = rep(NA, length(t_red)))

#set up loops
temp <- list(cumhar = numeric(), phunt = numeric(), mPop = numeric())
pick <- sample(1:length(out$sims.list$r.max), Nsamples)
for(t in 1:length(t_red)){
  for(i in 1:Nsamples){
    reward <- project.pop2(Tmax = TimeHor+1,
                           n1 = out$sims.list$N.tot[pick[i],], 
                           r = out$sims.list$r.max[pick[i]], 
                           theta = out$sims.list$theta[pick[i]],
                           K = out$sims.list$CC[pick[i]], 
                           Hgreen = out$sims.list$mu.green[pick[i]],
                           Hred = out$sims.list$m.har[pick[i]],
                           sdH = out$sims.list$sigma.har[pick[i]],
                           sdpop = out$sims.list$sigma.proc[pick[i]],
                           q = out$sims.list$q[pick[i]],
                           t_close = t_red[t],
                           total = TRUE)
    temp$cumhar[i] <- sum(reward$har, na.rm = TRUE)
    temp$phunt[i] <- sum(reward$hunt, na.rm = TRUE)/length(reward$hunt[!is.na(reward$hunt)])
    temp$mPop[i] <- mean(reward$pop, na.rm = TRUE)
  }
  df$cumHar[t] <- mean(temp$cumhar) #Use mean as surrogate utility function for now
  df$pHunt[t] <- mean(temp$phunt)
  df$mPop[t] <- mean(temp$mPop)
  df$cumHarSD[t] <- sd(temp$cumhar)
  df$pHuntSD[t] <- sd(temp$phunt)
  df$mPopSD[t] <- sd(temp$mPop)
}

saveRDS(df, file = "optim0.harOriginal.RDS")
df <- readRDS("optim0.harOriginal.RDS")
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
#yield curve
ggplot(data = df, aes(x=mPop, y=cumHar)) + 
  geom_point()+
  geom_smooth(method="gam", se=TRUE)
#find closure level
fit <- gam(cumHar~s(Closure), data=df)
pred <- predict(fit)
df[which(pred == max(pred)),]
