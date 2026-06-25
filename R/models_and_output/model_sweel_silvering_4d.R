# Swedish eel length at silvering model and mcmc 
# Viktor Thunell

library(nimbleHMC)
library(dplyr)

sweelsilv.code <- nimbleCode({
  
  # mean length of the data is 631: i.e. almost Durif 2005 FIII mean length. z = 0 with broad sd is ok
  mu_a_gl ~ dnorm(0, sd = 2)
  sd_a_gl ~ dgamma(shape = 2, rate = 6)
  mu_b_gl ~ dnorm(log(5), sd = 0.5)
  sd_b_gl ~ dgamma(shape = 2, rate = 6)
  
  for(j in 1:nsb){
    alpha_sb[j] ~ dnorm(mu_a_gl, sd = sd_a_gl)
    logb0_sb[j] ~ dnorm(mu_b_gl, sd = 0.5)
  }
  
  sd_a_su ~ dgamma(shape = 2, rate = 6)  # sd for su-level
  sd_b_su ~ dgamma(shape = 2, rate = 6)  # sd for su-level
  
  for(j in 1:nsu){
    alpha_raw[j] ~ dnorm(alpha_sb[susb[j]], sd = sd_a_su)  
    alpha[j] <- alpha_raw[j] * balanced[j] + alpha_sb[susb[j]] * (1 - balanced[j])
    logb0_raw[j] ~ dnorm(logb0_sb[susb[j]], sd = sd_b_su)  
    logb0[j] <- logb0_raw[j] * balanced[j] + logb0_sb[susb[j]] * (1 - balanced[j])
    b0[j] <- exp(logb0[j])
    }
  
  b1 ~ dnorm(0, sd = 0.5)
  b2 ~ dnorm(0, sd = 0.5)
  b3 ~ dnorm(0, sd = 0.5)
  for(j in 1:nsb){
     b4[j] ~ dnorm(0, sd = 0.5)
   }
  #b4 ~ dnorm(0, sd = 0.5)
  
  se_v ~ dbeta(1, 5) # p(classify as silver | is silver)
  sp_v ~ dbeta(1, 5) # p(classify as notsilver | is notsilver)
  se_d ~ dbeta(1, 10) 
  sp_d ~ dbeta(1, 10) 
  
  # fish base: mean a = 0.0010, mean b = 3.15, SD log10(W) = 0, SD log10(a) = 0.3030 SD b = 0.1811
  # sd of a on log10 scale sd is 0.303, 
  a_mean <- -0.7 # log(10) + log10(.001)
  b_mean <- 3.15
  a_sd <- 2 # log(10) + log10(.303)
  b_sd <- 0.1811*1.5 
  ab[1] ~ dnorm(a_mean, sd = a_sd)
  ab[2] ~ dnorm(b_mean, sd = b_sd)
  
  w_prec <- 1/w_sd^2
  w_sd ~ dexp(5)
  w_nu ~ dexp(0.1)
  
  # Likelihood
  for(i in 1:nobs){
    
    age_sc[i] ~ dnorm(0, 1)
    #length_sc[i] ~ dnorm(0, 1)
    
    z[i] <- alpha[su[i]] + b0[su[i]]*length_sc[i] + b1*age_sc[i] + b2*le_cren[i] + b3*lat_sc[i] + b4[sb[i]]*habitat[i]
    p[i] <- 1 / (1 + exp(-z[i]))
    silver[i] ~ dbern(p[i])
    
    # probability of each classifier to return 1, given true stage
    p_vis[i] <- silver[i]*se_v + (1 - silver[i])*(1 - sp_v)
    p_dur[i] <- silver[i]*se_d + (1 - silver[i])*(1 - sp_d)
    
    # classifiers observing the latent true stage
    silver_vis[i] ~ dbern(p_vis[i])
    silver_dur[i] ~ dbern(p_dur[i])
    
    # posterior predictive nodes
    silver_rep[i] ~ dbern(p[i])
    visual_rep[i]   ~ dbern(silver_rep[i] * se_v + (1 - silver_rep[i]) * (1 - sp_v))
    durif_rep[i] ~ dbern(silver_rep[i] * se_d + (1 - silver_rep[i]) * (1 - sp_d))
    
    #Condition le cren
    logweight[i] ~ dt(mu = logw_mu[i], tau = w_prec, df = w_nu)
    logw_mu[i] <- ab[1] + ab[2]*log(length[i])

    # calculate le cren (observed/expected weight)
    le_cren[i] <- exp(logweight[i]) / exp(logw_mu[i])
  }
  
})

nobs <- nrow(df.sweel3)
nsu <- length(unique(df.sweel3$su))
nsb <- length(unique(df.sweel3$sb))
susb <- df.sweel3 %>% distinct(su,sb) %>% arrange(su) %>% pull(sb)
balanced <- df.sweel3 %>% distinct(su,balance) %>% arrange(su) %>% pull(balance)
const <- list(nobs = nobs,
              sb = df.sweel3$sb,
              nsu = nsu,
              nsb = nsb,
              su = df.sweel3$su,
              susb = susb,
              balanced = balanced
              )

initdur = df.sweel3 %>%
  mutate(dna = if_else(is.na(silver_dur), 1, NA)) %>% pull(dna)
initvis = df.sweel3 %>%
  mutate(vna = if_else(is.na(silver_vis), 1, NA)) %>% pull(vna)
initlen = df.sweel3 %>% 
  mutate(lna = if_else(is.na(length_sc), 0, NA)) %>% pull(lna)
initage = df.sweel3 %>% 
  mutate(ana = if_else(is.na(age_sc), 0, NA)) %>% pull(ana)
initlogwe = df.sweel3 %>% 
  mutate(wna = if_else(is.na(logweight), 5, NA)) %>% pull(wna)

inits <- function() {list(
  alpha_raw = rnorm(const$nsu, 0, 1),
  logb0_raw  = rnorm(const$nsu, 2, .1),
  mu_a_gl = rnorm(1, 0, 1),
  sd_a_gl = rnorm(1, 0.1, 1),
  mu_b_gl = rnorm(1, log(5), 1),
  sd_b_gl = rnorm(1, 0.1, 0.1),
  sd_a_su = rnorm(1, 0.1, 0.1),
  sd_b_su = rnorm(1, 0.1, 0.1),
  silver = rbinom(const$nobs,size = 1, 0.5),
  silver_vis = initvis,
  silver_dur = initdur,
  length_sc = initlen,
  age_sc = initage,
  logweight = initlogwe,
  b1 = rnorm(1, 0, .1),
  b2 = rnorm(1, 0, .1),
  b3 = rnorm(1, 0, .1),
  b4 = rnorm(nsb, 0, .1),
  ab = rnorm(2, 0, .1)
)}

# build model
sweelsilv.model <- nimbleModel(sweelsilv.code,
                               constants = const,
                               inits=inits(),
                               data = df.sweel3 %>% select(silver_vis,silver_dur,length,length_sc,age_sc,lat_sc,habitat,logweight), 
                               buildDerivs = TRUE,
                               calculate = FALSE)
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
sweelsilv.c <- compileNimble(sweelsilv.model, resetFunctions = TRUE )

monits = c(mvars,"p","z","alpha","logb0","le_cren","logProb_silver_vis","logProb_silver_dur")
# configure and build mcmc and add hmc to alpha and sigma nodes
t <- Sys.time()
sweelsilv.confmcmc <- configureHMC(sweelsilv.c, monitors = monits, enableWAIC = TRUE)

sweelsilv.mcmc <- buildMCMC(sweelsilv.confmcmc, project = sweelsilv.model)
# compile mcmc
sweelsilv.mcmcc <- compileNimble(sweelsilv.mcmc, project = sweelsilv.model)
Sys.time() - t

t <- Sys.time()
sweelsilv.samples <- runMCMC(sweelsilv.mcmcc, niter = 3000, nburnin = 1500, nchains = 1, WAIC = TRUE)
Sys.time() - t

saveRDS(sweelsilv.samples, file = paste0(home,"/data/samples/sweelsilv.samples_d",Sys.Date(),".RData"))
