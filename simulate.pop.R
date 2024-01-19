# simulate theta-logistic to find equilibrium at posterior mean
# Notes: (1) should population model be in terms of harvest rate? Functional response of hunters? --> yes, see below
#        (2) add a revised observation model to reduce process variance, 
#           average observer effects and annual detection effects.
# 20210916 (1) revise to do forward sims based on harvest rate. 
#              This allows harvest to decline with population size. 
# 20221116 (1) commented out dev models and worked on stochastic model
#          (2) revised to fully sample posterior initial population size 
#              (N.tot in 2022)
# 20230629 (1) revised prediction model to reflect new (2023) jags model; 
#              this include observer-specific observation model and the current 
#              decision index of the average estimate between observers 
# 2023???? (1) added permit harvest to model
# 20240118 worked to add multiple harvest strategies to population projection model,
#          specifically:
#          (1) deleted commented out dev code
#          (2) renamed function to project_pop
#          (3) moved functional response of harvest, (1 - exp(-0.0001*pop[t-1])), 
#              outside population projection step so that har[t] is realized harvest
#          (4) modified function to define a management strategy:
#                did this by making Hred and Hgreen a vector that defines harvest
#                and t_close another vector that defines the thresholds

project_pop <- function(
    #these inputs are scalar, 
  #  can use a sample of posterior or a summary value (posterior mean)
  #default is simulation at posterior mean
  #also assume that the observers are constant for observation model 
  # even when stochastic = FALSE, 
  Tmax = 100,
  n1 = out$mean$N.tot, 
  r = out$mean$r.max, 
  theta = out$mean$theta,
  K = out$mean$CC, 
  Hlist = list(red = out$mean$m.har, 
               yellow = out$mean$m.har + (out$mean$mu.green - out$mean$m.har)/2, 
               green = out$mean$mu.green), 
  sdH = list(red = 0, yellow = 0, green = 0),
  threshold = list(red = 23000, green = 28000), #closure threshold, 23000 is 2016 management plan default
  Hp = list(green = out$mean$m.har.p), #permit harvest
  sdHp = list(green = 0), #sd permit harvest
  sdpop = 0,
  q = out$mean$q, 
  alpha1 = c(out$mean$alpha1[5], out$mean$alpha1[10]), #default is HMW and MAS
  sdesp = out$mean$sd.esp, 
  BETA = NA, #linear model slope
  SE = NA, #linear model slope SE
  SIGMA = NA, #linear model RMSE
  DF = NA, #linear model degrees of freedom
  crip = out$mean$c, #crippling probability
  stochastic = TRUE){
  
  pop <- har <- hunt <- obs <- numeric(Tmax+1) #rep(0, Tmax)
  pop[1] <- n1[length(n1)] #starting population size
  s <- rchisq(1, DF) #linear model for observation error
  beta <- rnorm(1, BETA, SE)
  esp <- rnorm(length(alpha1), 0, sdesp)
  ############
  #determine harvest in year 1, har[1] #should turn into a function
  #assumune population is large and harvest is "green"
  #Note: har[t] is the harvest from april of t-1 to september of t; census date is June 
  har[1] <- rnorm(1, Hlist[[3]], sdH[[3]])+rnorm(1, Hp[[1]], sdHp[[1]])
  har[1] <- ifelse(har[1] < 0, 0, har[1]) #check that harvest is not negative
  ############
  
  for(t in 2:(Tmax+1)){
    #project population
    pop[t] <- pop[t-1] + pop[t-1]*r*(1-(pop[t-1]/K)^theta) - 
      har[t-1]/(1-crip)
    pop[t] <- max(pop[t], 1e-5)
    if(stochastic == TRUE){pop[t] <- rlnorm(1, log(pop[t]), sdpop)}
    if(pop[t] < 1){ 
      pop[t] <- 0
      break
    }
    #determine next years hunt
    mu <- q*pop[t-1] + alpha1 + esp
    sigma.obs <- rnorm(2, beta*mu, sqrt(((SIGMA^2)*DF)/s) )
    sigma.obs <- ifelse(sigma.obs < 0, 0, sigma.obs)
    obs[t-1] <- mean(rnorm(2, mu, sigma.obs))
    hunt[t] <- ifelse(obs[t-1] < threshold[[1]], 1,  #monitor hunting season type
                        ifelse(obs[t-1] < threshold[[2]], 2, 3))
    #Note: har[t] is the harvest from april of t-1 to september of t; census date is June 
    har[t] <- ifelse(obs[t-1] < threshold[[1]], rnorm(1, Hlist[[1]], sdH[[1]]), 
                       ifelse(obs[t-1] < threshold[[2]], 
                              rnorm(1, Hlist[[2]], sdH[[2]])+rnorm(1, Hp[[1]], sdHp[[1]]), 
                              rnorm(1, Hlist[[3]], sdH[[3]])+rnorm(1, Hp[[1]], sdHp[[1]]))
                       )
    har[t] <- ifelse(har[t] < 0, 0, har[t]) #check that harvest is not negative
    har[t] <- har[t]*(1 - exp(-0.0001*pop[t])) #functional response
  }
  #determine last obs
  mu <- q*pop[Tmax+1] + alpha1 + esp
  sigma.obs <- rnorm(2, beta*mu, sqrt(((SIGMA^2)*DF)/s) )
  sigma.obs <- ifelse(sigma.obs < 0, 0, sigma.obs)
  obs[Tmax+1] <- mean(rnorm(2, mu, sigma.obs))
  return(list(pop = pop[-1], obs = obs[-1], har = har[-1], hunt = hunt[-1]))
}
# #test it
# out <- readRDS("data/out.RDS")
# #load linear observation model
# fit <- readRDS("data/fit.RDS")
# #also need to run head of "theta.logistic.R to get fit object
# pop <- project_pop( BETA = fit$coefficients[1], #linear model slope
#                     SE = fit$coefficients[2], #linear model slope SE
#                     SIGMA = fit$sigma, #linear model RMSE
#                     DF = fit$df[2], 
#                     sdpop = out$mean$sigma.proc, 
#                     sdHp = list(green = out$mean$sd.har.p), 
#                     sdH = list(red = out$mean$sigma.har, 
#                                yellow = out$mean$sigma.har, 
#                                green = out$mean$sigma.har))
# pop <- as.data.frame(pop)
# pcols <- c("red", "orange", "green")
# plot(1:100, pop$pop, type = "b", pch = 16, col = pcols[pop$hunt])
# plot(1:100, pop$obs, type = "b", pch = 16, col = pcols[pop$hunt])
# abline(h=c(23000, 28000))
# plot(1:100, pop$har, type = "b", pch = 16, col = pcols[pop$hunt])
# #does function work if we set red and green threshold the same?
# pop <- project_pop( BETA = fit$coefficients[1], #linear model slope
#                     SE = fit$coefficients[2], #linear model slope SE
#                     SIGMA = fit$sigma, #linear model RMSE
#                     DF = fit$df[2], 
#                     sdpop = out$mean$sigma.proc, 
#                     threshold = list(red = 23000, green = 23000),
#                     sdHp = list(green = out$mean$sd.har.p), 
#                     sdH = list(red = out$mean$sigma.har, 
#                                yellow = out$mean$sigma.har, 
#                                green = out$mean$sigma.har))
# plot(1:100, pop$pop, type = "b", pch = 16, col = pcols[pop$hunt])
# plot(1:100, pop$obs, type = "b", pch = 16, col = pcols[pop$hunt])
# abline(h=c(23000, 28000))
# plot(1:100, pop$har, type = "b", pch = 16, col = pcols[pop$hunt])
# #yes!
