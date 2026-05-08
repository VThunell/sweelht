# DIASPARA WP2.2 LHT models - Viktor Thunell
# Eel length at silvering model and mcmc

library(nimbleHMC)

sweelsilv.code <- nimbleCode({
  
  # prior par. from Durif 2005 indiv. mean, https://doi.org/10.1111/j.0022-1112.2005.00662.x
  # sdsl <- sqrt(log(1 + (124/658)^2))  
  # msl <- log(658) - sdsl^2/2
  # sdbsl <- sqrt(log(1 + (23/393)^2))
  # mbsl <- log(393/658) - sdbsl^2/2
  # 
  # # overall "global" silver length
  # mu_gl ~ dnorm(msl,sdsl)
  # sd_gl ~ dgamma(shape = 1.5, rate = 1.5 / sdsl)
  # # group means 
  # mu_sig ~ dgamma(shape = 1.5, rate = 1.5 / sdsl) # group sd mean 
  # sd_sig ~ dgamma(shape = 1.5, rate = 1.5 / sdsl*0.75) # sd of the group sd means
  # 
  # # group basin estimates
  # for(j in 1:nmb){
  #   alpha[j] ~ dnorm(mu_gl, sd = sd_gl)
  #   sig_s[j] ~ dlnorm(meanlog = log(mu_sig) - sd_sig^2/2, sdlog = sd_sig)
  # }
  # 
  # bs ~ dnorm(mbsl, sd = sdbsl)
  alpha ~ dnorm(0, sd = 5)
  b0 ~ dnorm(0, sd = 5)
  # 
  # ps[1] ~ dbeta(1,1)  # prop of males (prob. to be 1)
  # ps[2] <- 1-ps[1]
  # 
  # Likelihood for heteroscedastic model
  for(i in 1:nobs){
    silver[i] ~ dbern(prob = p_s[i])
    p_s[i] <- 1 / (1 + exp(-z[i])) #logit link
    z[i] <- alpha + b0[mb[i]] * length_sc[i]
  
    L50 <- -alpha/b0
    #mu_s[i] <- alpha[mb[i]] + bs*sex[i] #+ bt*temp.sc[i]
    
    #sex[i]~dbern(ps[1])
  }
  
  # 
  # # first year alpha and bs
  # for(l in 1:ner){
  #   alpha[l,1] ~ dnorm(mu.a.gl[l], sd = sd.a.gl[l]) 
  #   bs[l] ~ dnorm(mbsl, sd = sdbsl)
  # }
  # 
  # # transition of pars from year i to i+1
  # for(i in 1:(nyear-1)){
  #   alpha[1:ner, (i+1)] ~ dmnorm(alpha[1:ner, i], prec = tau_p[1:ner,1:ner])
  # }
  # 
  # for(k in 1:ner){
  #   for(j in 1:ner){
  #     sigma_p[k,j] <- Rnew[k,j] * sd.a.gl[k] * sd.a.gl[j]  
  #   }
  # }
  # 
  # tau_p[1:ner, 1:ner]<- inverse(sigma_p[1:ner, 1:ner])
  
  # #Prior for correlation matrix Rnew (LKJ prior) on the cholesky of the Rnew, R. 
  # phi[1]  <- eta + (ner - 2)/2
  # corY[1] ~ dbeta(phi[1], phi[1])
  # r12   <- 2 * corY[1] - 1
  # ##
  # R[1,1]     <- 1
  # R[1,2]     <- r12
  # R[2,2]     <- sqrt(1 - r12^2)
  # 
  # R[2:ner,1]   <- 0
  # 
  # for (m in 2:(ner-1)) {
  #   ## Draw beta random variable
  #   phi[m] <- phi[(m-1)] - 0.5
  #   corY[m] ~ dbeta(m / 2, phi[m])
  #   ## Draw uniformly on a hypersphere
  #   for (jj in 1:m) {
  #     corZ[m, jj] ~ dnorm(0, 1)
  #   }
  #   scZ[m, 1:m] <- corZ[m, 1:m] / sqrt(inprod(corZ[m, 1:m], corZ[m, 1:m]))
  #   R[1:m,(m+1)] <- sqrt(corY[m]) * scZ[m,1:m] 
  #   R[(m+1),(m+1)] <- sqrt(1 - corY[m])
  #   for(jk in (m+1):ner){
  #     R[jk,m] <- 0 
  #   }
  # }  #m
  # 
  # Rnew[1:ner,1:ner] <- t(R[1:ner,1:ner]) %*% R[1:ner,1:ner]
  # 
  
  
})

initsilv = df.sweel2 %>% 
  mutate(sna = if_else(is.na(silver), 1, NA)) %>% pull(sna)

inits <- function() {
  list(
    alpha = rnorm(1, 0, 1),
    b0    = rnorm(1, 0, 1),
    p_s   =  runif(const$nobs,0.1,0.9),
    silver  = initsilv
  )
}

# initial values generating function for all nodes
# inits <- function() {
#   list(s.sigma = rexp(1,5),
#        alpha = matrix(rnorm(const$ner * const$nyear, log(658), 0.05),
#                       nrow = const$ner, ncol = const$nyear),
#        bs = rnorm(const$ner, log(393/658), 0.05),
#        bt = rnorm(1, 0, 0.1),
#        mu.a.gl  = rnorm(const$ner, log(658), 0.05),
#        sd.a.gl  = runif(const$ner, 0.2, 0.5),
#        corZ = matrix(rnorm((const$ner-1)*(const$ner-1), 0, 1),
#                      nrow = const$ner-1, ncol = const$ner-1),
#        corY = runif((const$ner-1), 0, 1)
#   )
# }

nobs <- nrow(df.sweel2)
#nmb <- length(unique(df.sweel$main_bas))
#nyear <- length(unique(df$year))

const <- list(nobs = nobs
              #eta = 2,
              #ner = ner,
              #nmb = nmb,
              #mb = df.sweel$mb
              #nyear = nyear,
              # temp.sc = data.silv$temp.sc,
              # year = data.silv$year,
              # er = data.silv$er
              )

# build model
sweelsilv.model <- nimbleModel(sweelsilv.code,
                               constants = const,
                               inits=inits(),
                               #data = df.sweel %>% select(length1,sex), buildDerivs = TRUE)
                               data = df.sweel2 %>% select(silver,length_sc), buildDerivs = TRUE)

sweelsilv.model$simulate()
sweelsilv.model$calculate()
sweelsilv.model$initializeInfo()

dataNodes <- sweelsilv.model$getNodeNames(dataOnly = TRUE)
parentNodes <- sweelsilv.model$getParents(dataNodes, stochOnly = TRUE) #all of these should be added to monitor below to recreate other model variables...
stnodes <- sweelsilv.model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
allvars <- sweelsilv.model$getVarNames(nodes = stnodes)
mvars <- allvars[!(grepl("lifted",allvars))]

# calculate vars to id NAs
vs <- mvars
for(i in 1:length(vs)){
  print(paste0(vs[i]," ",sweelsilv.model$calculate(vs[i]) ))
}

# compile model
sweelsilv.c <- compileNimble(sweelsilv.model)

monits = c(mvars,"p_s")

# configure and build mcmc and add hmc to alpha and sigma nodes
sweelsilv.confmcmc <- configureHMC(sweelsilv.c, monitors = monits, enableWAIC = TRUE)

sweelsilv.mcmc <- buildMCMC(sweelsilv.confmcmc, project = sweelsilv.model)

# compile mcmc
sweelsilv.mcmcc <- compileNimble(sweelsilv.mcmc, project = sweelsilv.model)

sweelsilv.samples <- runMCMC(sweelsilv.mcmcc, niter = 2000, nburnin = 1000, nchains = 1, thin = 1, WAIC=TRUE)

#saveRDS(sweelsilv.samples, file = paste0(home,"/models_eel/samples/sweelsilv.samples_",Sys.Date(),".RData"))
