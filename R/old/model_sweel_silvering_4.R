# Swedish eel length at silvering model and mcmc 
# Viktor Thunell


library(nimbleHMC)
library(dplyr)

sweelsilv.code <- nimbleCode({
  
  # prior par. Assuming L50 corresponsd to phase IV in Durif 2005 for females indiv. mean, https://doi.org/10.1111/j.0022-1112.2005.00662.x
  # overall "global" silver length
  mu_a_gl ~ dnorm(msa, sda)
  sd_a_gl ~ dgamma(shape = 2, rate = 4)
  #mu_b_gl ~ dnorm(log(mb0), sd = log(sdb0))
  mu_b_gl ~ dnorm(log(mb0), sd = 0.5)
  #sd_b_gl ~ dgamma(shape = 1.5, rate = 1.5 / sdb0) # to weak 
  sd_b_gl ~ dgamma(shape = 2, rate = 6) 
  
  sd_a_su ~ dgamma(shape = 2, rate = 20)   # sd for su-level
  sd_b_su ~ dgamma(shape = 2, rate = 20)   # sd for su-level
  
  sb_alpha[1] ~ dnorm(mu_a_gl, sd = sd_a_gl)
  sb_alpha[2] ~ dnorm(mu_a_gl, sd = sd_a_gl)
  
  sb_logb0[1] ~ dnorm(mu_b_gl, sd = sd_b_gl)
  sb_logb0[2] ~ dnorm(mu_b_gl, sd = sd_b_gl)
  
  for(j in 1:nsu){
    alpha[j] ~ dnorm(sb_alpha[susb[j]], sd = sd_a_su)
    #alpha[j] ~ dnorm(mu_a_gl, sd = sd_a_gl) # for non-nested alpha
    log_b0[j] ~ dnorm(sb_logb0[susb[j]], sd = sd_b_su)
    #log_b0[j] ~ dnorm(mu_b_gl, sd = sd_b_gl) # for non-nested alpha
    b0[j] <- exp(log_b0[j])
  }
  
  b1 ~ dnorm(0, sd = 0.5)
  
  se_v ~ dbeta(5, 2) # p(classify as silver | is silver)
  sp_v ~ dbeta(5, 2) # p(classify as notsilver | is notsilver)
  se_d ~ dbeta(10, 1) 
  sp_d ~ dbeta(10, 1) 
  
  # Likelihood
  for(i in 1:nobs){
    
    age_sc[i] ~ dnorm(0, 1)
    length_sc[i] ~ dnorm(0, 1)
    
    z[i] <- alpha[su[i]] + b0[su[i]]*length_sc[i] + b1*age_sc[i]
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
    
    # Condition le cren
    # logweight[i] ~ dt(mu = logw.mu[i], tau = w.prec, df = w.nu) 
    # logw.mu[i] <- ab[mb[i],year[i],1] + ab[mb[i],year[i],2]*log(length_sc[i])  
    # 
    # # calculate le cren
    # le.cren[i] <- weight[i] / exp(logw.mu[i])
  }
  
  # w.prec <- 1/w.sd^2
  # w.sd ~ dexp(7)
  # w.nu ~ dexp(0.1)
})

nobs <- nrow(df.sweel3)
nsu <- length(unique(df.sweel3$su))
susb <- df.sweel3 %>% distinct(su,sb) %>% arrange(su) %>% pull(sb)

const <- list(nobs = nobs,
              msa = (746 - mean(df.sweel3$length, na.rm = TRUE)) / sd(df.sweel3$length, na.rm = TRUE),
              sda = 150/sd(df.sweel3$length, na.rm = TRUE),
              mb0 = 5,
              #sdb0 = 1,
              nsu = nsu,
              su = df.sweel3$su,
              susb = susb
)

initdur = df.sweel3 %>%
  mutate(dna = if_else(is.na(silver_dur), 1, NA)) %>% pull(dna)
initvis = df.sweel3 %>%
  mutate(vna = if_else(is.na(silver_vis), 1, NA)) %>% pull(vna)
initlen = df.sweel3 %>% 
  mutate(lna = if_else(is.na(length_sc), 0, NA)) %>% pull(lna)
initage = df.sweel3 %>% 
  mutate(ana = if_else(is.na(age_sc), 0, NA)) %>% pull(ana)

inits <- function() {list(
  alpha = rnorm(const$nsu, 0, 1),
  log_b0  = rnorm(const$nsu, 0, .1),
  log_b0 = rep(log(const$mb0), const$nsu),
  mu_a_gl = rnorm(1, 500, 1),
  sd_a_gl = rnorm(1, 50, 1),
  mu_b_gl = rnorm(1, log(1), 1),
  sd_b_gl = rnorm(1, 1, 0.1),
  p = runif(const$nobs,0.1,0.9),
  silver_vis = initvis,
  silver_dur = initdur,
  length_sc = initlen,
  age_sc = initage
)}

# build model
sweelsilv.model <- nimbleModel(sweelsilv.code,
                               constants = const,
                               inits=inits(),
                               data = df.sweel3 %>% select(silver_vis,silver_dur,length_sc,age_sc), buildDerivs = TRUE)
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
t <- Sys.time()
sweelsilv.confmcmc <- configureHMC(sweelsilv.c, monitors = monits, enableWAIC = TRUE)

sweelsilv.mcmc <- buildMCMC(sweelsilv.confmcmc, project = sweelsilv.model)
# compile mcmc
sweelsilv.mcmcc <- compileNimble(sweelsilv.mcmc, project = sweelsilv.model)
Sys.time() - t

t <- Sys.time()
sweelsilv.samples <- runMCMC(sweelsilv.mcmcc, niter = 3000, nburnin = 1500, nchains = 1, WAIC = TRUE)
Sys.time() - t

saveRDS(sweelsilv.samples, file = paste0(home,"/data/samples/sweelsilv.samples_",Sys.Date(),".RData"))
