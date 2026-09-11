# Swedish eel length at silvering model and mcmc 
# Viktor Thunell

library(nimbleHMC)
library(dplyr)

sweelsilv.code <- nimbleCode({
  
  # mean length of the data is 631: i.e. almost Durif 2005 FIII mean length. z = 0 with broad sd is ok
  mu_a_gl ~ dnorm(0, sd = 2)
  mu_b_gl ~ dnorm(log(3), sd = 0.5)
  
  sd_a_su ~ dgamma(shape = 2, rate = 6)  # sd for su-level
  sd_b_su ~ dgamma(shape = 2, rate = 6)  # sd for su-level
  
  for(j in 1:nsu){
    alpha_raw[j] ~ dnorm(mu_a_gl, sd = sd_a_su)
    alpha[j] <- alpha_raw[j] * balanced[j] + mu_a_gl * (1 - balanced[j])
    logb0_raw[j] ~ dnorm(mu_b_gl, sd = sd_b_su)
    logb0[j] <- logb0_raw[j] * balanced[j] + mu_b_gl * (1 - balanced[j])
    b0[j] <- exp(logb0[j])
    }
  
  bA ~ dnorm(1, sd = 2)
  bC ~ dnorm(-1, sd = 2)
  bE ~ dnorm(-1, sd = 2)
  bP ~ dnorm(-1, sd = 2)
  bLa ~ dnorm(-1, sd = 2)
  for(j in 1:nsb){
     bH[j] ~ dnorm(0, sd = 2)
   }
  
  se_v ~ dbeta(10, 2) # p(classify as silver | is silver)
  sp_v ~ dbeta(10, 2) # p(classify as notsilver | is notsilver)
  
  # Likelihood
  for(i in 1:nobs){
    
    # for imputation
    age_sc[i] ~ dnorm(0,1)
    con_index_sc[i] ~ dnorm(0,1)
    eye_index_sc[i] ~ dnorm(0,1)
    pec_index_sc[i] ~ dnorm(0,1)
    
    # logistic
    z[i] <- alpha[su[i]] + b0[su[i]]*length_sc[i] + bA*age_sc[i] + 
      bC*con_index_sc[i] + bE*eye_index_sc[i] + 
      bP*pec_index_sc[i] + bLa*lat_sc[i] + bH[sb[i]]*habitat[i]
    
    p[i] <- 1 / (1 + exp(-z[i]))
    silver[i] ~ dbern(p[i])
    
    # probability of visual classifier to return 1, given true stage
    p_vis[i] <- silver[i]*se_v + (1 - silver[i])*(1 - sp_v)
    
    # classifier observing the latent true stage
    silver_vis[i] ~ dbern(p_vis[i])
    
    # posterior predictive nodes
    silver_pp[i] ~ dbern(p[i])
    visual_pp[i] ~ dbern(silver_pp[i] * se_v + (1 - silver_pp[i]) * (1 - sp_v))
    
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
              #susb = susb,
              balanced = balanced
              )

initvis = df.sweel3 %>%
   mutate(vna = if_else(is.na(silver_vis), 1, NA)) %>% pull(vna)
initage = df.sweel3 %>% 
  mutate(ana = if_else(is.na(age_sc), 0, NA)) %>% pull(ana)
initcon = df.sweel3 %>% 
  mutate(cna = if_else(is.na(con_index_sc), 1, NA)) %>% pull(cna)
initpect = df.sweel3 %>% 
  mutate(pna = if_else(is.na(pec_index_sc), 1, NA)) %>% pull(pna)
initeye = df.sweel3 %>% 
  mutate(ena = if_else(is.na(eye_index_sc), 1, NA)) %>% pull(ena)


inits <- function() {list(
  alpha_raw = rnorm(const$nsu, 0, 1),
  logb0_raw  = rnorm(const$nsu, 2, .1),
  mu_a_gl = rnorm(1, 0, 1),
  mu_b_gl = rnorm(1, log(5), 1),
  sd_a_su = rnorm(1, 0.1, 0.1),
  sd_b_su = rnorm(1, 0.1, 0.1),
  silver = rbinom(const$nobs,size = 1, 0.5),
  silver_vis = initvis,
  age_sc = initage,
  con_index_sc = initcon,
  pec_index_sc = initpect,
  eye_index_sc = initeye,
  bA = rnorm(1, 0, .1),
  bC = rnorm(1, 0, .1),
  bLa = rnorm(1, 0, .1),
  bP = rnorm(1, 0, .1),
  bE = rnorm(1, 0, .1),
  bH = rnorm(nsb, 0, .1)
)}

# build model
sweelsilv.model <- nimbleModel(sweelsilv.code,
                               constants = const,
                               inits=inits(),
                               data = df.sweel3 %>% select(silver_vis,length_sc,age_sc,con_index_sc,lat_sc,habitat,eye_index_sc,pec_index_sc), 
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

monits = c(mvars,"p","z","alpha","logb0","logProb_silver_vis","p_vis")
# configure and build mcmc and add hmc to alpha and sigma nodes
t <- Sys.time()
sweelsilv.confmcmc <- configureHMC(sweelsilv.c, monitors = monits, enableWAIC = TRUE)

sweelsilv.mcmc <- buildMCMC(sweelsilv.confmcmc, project = sweelsilv.model)
# compile mcmc
sweelsilv.mcmcc <- compileNimble(sweelsilv.mcmc, project = sweelsilv.model)
Sys.time() - t

t <- Sys.time()
sweelsilv.samples <- runMCMC(sweelsilv.mcmcc, niter = 2000, nburnin = 1000, nchains = 1, WAIC = TRUE)
Sys.time() - t

# test current_filename()
saveRDS(sweelsilv.samples, file = paste0(home,"/data/samples/sweelsilv.samples_i",Sys.Date(),".RData"))
