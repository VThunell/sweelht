# DIASPARA WP2.2 LHT models - Viktor Thunell
# Eel length at silvering model and mcmc

library(nimbleHMC)

sweelsilv.code <- nimbleCode({
  
  # # prior par. Assuming L50 corresponsd to phase III in Durif 2005 for females indiv. mean, https://doi.org/10.1111/j.0022-1112.2005.00662.x
  # sdsl <- msl_sd
  # msl <- msl_p
  mbsl <- 393/658
  #sdbsl <- mbsl*0.5
  # 
  # overall "global" silver length
  mu_a_gl ~ dnorm(msa, sda)
  sd_a_gl ~ dgamma(shape = 1.5, rate = 1.5 / sda)
  mu_b_gl ~ dnorm(log(mb0), sd = log(sdb0))
  #mean(exp(rnorm(10,log(1), sd = log(2))))
  sd_b_gl ~ dgamma(shape = 1.5, rate = 1.5 / sdb0)
  
  for(j in 1:nmb){
    alpha[j] ~ dnorm(mu_a_gl, sd = sd_a_gl)
    log_b0[j] ~ dnorm(mu_b_gl, sd = sd_b_gl)
    b0[j] <- exp(log_b0[j])
  }
  
  #ps[1] ~ dbeta(1,1)  # prop of males (prob. to be 1)
  #ps[2] <- 1-ps[1]
  
  #bs ~ dnorm(mbsl, sd = mbsl*0.3)
  
  # Likelihood
  for(i in 1:nobs){
    silver[i] ~ dbern(prob = p[i])
    p[i] <- 1 / (1 + exp(-z[i])) #logit link
    z[i] <- alpha[mb[i]] + b0[mb[i]] * length_sc[i] #+ bs*sex[i] #+ bt*temp.sc[i]
  
    #sex[i] ~ dbern(ps[1])
  }
  
})

nobs <- nrow(df.sweel2)
nmb <- length(unique(df.sweel2$main_bas))

const <- list(nobs = nobs,
              msa = (658 - mean(df.sweel2$length1)) / sd(df.sweel2$length1),
              sda = 124/sd(df.sweel2$length1),
              mb0 = 1,
              sdb0 = 2,
              nmb = nmb,
              mb = df.sweel2$mb
              )

initsilv = df.sweel2 %>%
  mutate(sna = if_else(is.na(silver), 1, NA)) %>% pull(sna)
initsex = df.sweel2 %>% 
  mutate(sa = if_else(is.na(sex), 0, NA)) %>% pull(sa)

inits <- function() {
  list(
    alpha = rnorm(const$nmb, 0, 1),
    log_b0  = rnorm(const$nmb, 0, .1),
    #ps   = rnorm(2,0.5,.1),
    #bs = rnorm(1,0.5,.1),
    log_b0 = rep(log(const$mb0), const$nmb),
    mu_a_gl = rnorm(1, 500, 1),
    sd_a_gl = rnorm(1, 50, 1),
    mu_b_gl = rnorm(1, log(1), 1),
    sd_b_gl = rnorm(1, 1, 0.1),
    p   = runif(const$nobs,0.1,0.9),
    silver  = initsilv
    #sex = initsex
  )
}

# build model
sweelsilv.model <- nimbleModel(sweelsilv.code,
                               constants = const,
                               inits=inits(),
                               data = df.sweel2 %>% select(silver,length_sc), buildDerivs = TRUE)
#sweelsilv.model$silver[27811]
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

monits = c(mvars,"p","z")

# configure and build mcmc and add hmc to alpha and sigma nodes
sweelsilv.confmcmc <- configureHMC(sweelsilv.c, monitors = monits, enableWAIC = TRUE)

sweelsilv.mcmc <- buildMCMC(sweelsilv.confmcmc, project = sweelsilv.model)

# compile mcmc
sweelsilv.mcmcc <- compileNimble(sweelsilv.mcmc, project = sweelsilv.model)

sweelsilv.samples <- runMCMC(sweelsilv.mcmcc, niter = 2000, nburnin = 1000, nchains = 1, WAIC=TRUE)

saveRDS(sweelsilv.samples, file = paste0(home,"/data/samples/sweelsilv.samples_",Sys.Date(),".RData"))
