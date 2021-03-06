######################################################################
#    Amplified cyclicality in mast seeding dynamics positively       #
#    influences the dynamics of a seed consumer species              #
######################################################################


### 2. SHORT TERM SIMULATIONS


# Autocorrelation 
# 2 environmental states
# 3 x 3 matrix
# 4 frequencies
# 11 autocorrelation coefficients (including 0!)
# 1,000 simulations of 150 time steps
# 50 first iterations discarded


library(popbio)

rm(list=ls(all=TRUE))

frequency <- c(0.35,0.45,0.55,0.65) # frequency of good years
autocorrelation <- c(-0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50)          # coefficient of autocorrelation
trun = 150                 # run each simulation for trun years
env = c("good","bad")       # the different environments 
env.n = 1:2
rep.simul = 1000          # number of repetition for each condition (frequency x autocorrelation)

time_steps = trun - 50

rows = 3
cols = 3
nstage = 3

f = frequency
v1 = autocorrelation

env_at_t = 1                                # state of the environment


### Calculation of the growth rates 

lambda_t <- array(0,c(length(frequency), length(autocorrelation),time_steps,rep.simul)) # store growth rate at each time step

for (z in 1:length(frequency)) {
  for (y in 1:length(autocorrelation)) {
    
    q = f[z]*(1-v1[y])
    p = v1[y]+q
    
    P = matrix(c(p,1-p,
                 q,1-q),nrow=2,ncol=2,byrow=F)
    
    colnames(P) <- c("1","2")
    row.names(P) <- c("1","2")
    
    for(rep in 1:rep.simul) {
      
      
      #### 1. STORE MATRICES FOR EACH SIMULATION 
      A_all <- array(NA, c(3, 3, trun)) # store all matrices 
      n <- t(rep((1/nstage),nstage))
      
      simul.n = numeric(trun)
      simul.n[1] = 1                   # define environmental state at time 1
      
      # initiate A_all with Matrix A
      Ssf = 1 - rnorm(1, mean = 0.022, sd = 0.002)
      Smf = 1 - rnorm(1, mean = 0.080, sd = 0.007)
      Slf = 1 - rnorm(1, mean = 0.083, sd = 0.006)
      
      pSS = rnorm(1, mean = 0.137, sd = 0.014)
      pSM = rnorm(1, mean = 0.284, sd = 0.025)
      pSL = rnorm(1, mean = 0.555, sd = 0.028)
      pML = rnorm(1, mean = 0.490, sd = 0.036)
      
      BPs_A = rnorm(1, mean = 0.178, sd = 0.016)
      BPm_A = rnorm(1, mean = 0.680, sd = 0.018)
      BPl_A = rnorm(1, mean = 0.842, sd = 0.019)
      
      Fs = rnorm(1, mean = 3.917, sd = 0.250)
      Fm = rnorm(1, mean = 4.708, sd = 0.105)
      Fl = rnorm(1, mean = 6.095, sd = 0.121)
      
      h = 0.425
      Spn = 0.750
      piOs = 0.600
      
      Matrix_A <- matrix(c(
        BPs_A * Fs * 0.5 * Spn * piOs * Ssf * (1 - h) + pSS * Ssf * (1 - h), BPm_A * Fm * 0.5 * Spn * piOs * Ssf * (1 - h), BPl_A * Fl * 0.5 * Spn * piOs * Ssf * (1 - h), 
        BPs_A * Fs * 0.5 * Spn * (1 - piOs) * Smf * (1 - h) + pSM * Smf * (1 - h), BPm_A * Fm * 0.5 * Spn * (1 - piOs) * Smf * (1 - h) + (1 - pML) * Smf * (1 - h), BPl_A * Fl * 0.5 * Spn * (1 - piOs) * Smf * (1 - h),  
        pSL * Slf * (1 - h), pML * Slf * (1 - h), Slf * (1 - h)),nrow = 3, byrow = T)      
      
      A_all[,,1] <- Matrix_A
      
      for (t in 2:trun) {
        
        simul.n[t] = env_at_t = sample(ncol(P),1,pr = P[,env_at_t]) # sample environmental state
        states = simul.n[t]
        
        # sample demographic parameters
        Ssf = 1 - rnorm(1, mean = 0.022, sd = 0.002)
        Smf = 1 - rnorm(1, mean = 0.080, sd = 0.007)
        Slf = 1 - rnorm(1, mean = 0.083, sd = 0.006)
        
        pSS = rnorm(1, mean = 0.137, sd = 0.014)
        pSM = rnorm(1, mean = 0.284, sd = 0.025)
        pSL = rnorm(1, mean = 0.555, sd = 0.028)
        pML = rnorm(1, mean = 0.490, sd = 0.036)
        
        BPs_A = rnorm(1, mean = 0.178, sd = 0.016)
        BPm_A = rnorm(1, mean = 0.680, sd = 0.018)
        BPl_A = rnorm(1, mean = 0.842, sd = 0.019)
        
        BPs_N = rnorm(1, mean = 0.155, sd = 0.012)
        BPm_N = rnorm(1, mean = 0.501, sd = 0.016)
        BPl_N = rnorm(1, mean = 0.602, sd = 0.019)
        
        Fs = rnorm(1, mean = 3.917, sd = 0.250)
        Fm = rnorm(1, mean = 4.708, sd = 0.105)
        Fl = rnorm(1, mean = 6.095, sd = 0.121)
        
        h = 0.425
        Spn = 0.750
        piOs = 0.600
        
        # compute matrices
        Matrix_A <- matrix(c(
          BPs_A * Fs * 0.5 * Spn * piOs * Ssf * (1 - h) + pSS * Ssf * (1 - h), BPm_A * Fm * 0.5 * Spn * piOs * Ssf * (1 - h), BPl_A * Fl * 0.5 * Spn * piOs * Ssf * (1 - h), 
          BPs_A * Fs * 0.5 * Spn * (1 - piOs) * Smf * (1 - h) + pSM * Smf * (1 - h), BPm_A * Fm * 0.5 * Spn * (1 - piOs) * Smf * (1 - h) + (1 - pML) * Smf * (1 - h), BPl_A * Fl * 0.5 * Spn * (1 - piOs) * Smf * (1 - h),  
          pSL * Slf * (1 - h), pML * Slf * (1 - h), Slf * (1 - h)),nrow = 3, byrow = T)      
        
        Matrix_N <- matrix(c(
          BPs_N * Fs * 0.5 * Spn * piOs * Ssf * (1 - h) + pSS * Ssf * (1 - h), BPm_N * Fm * 0.5 * Spn * piOs * Ssf * (1 - h), BPl_N * Fl * 0.5 * Spn * piOs * Ssf * (1 - h), 
          BPs_N * Fs * 0.5 * Spn * (1 - piOs) * Smf * (1 - h) + pSM * Smf * (1 - h), BPm_N * Fm * 0.5 * Spn * (1 - piOs) * Smf * (1 - h) + (1 - pML) * Smf * (1 - h), BPl_N * Fl * 0.5 * Spn * (1 - piOs) * Smf * (1 - h),  
          pSL * Slf * (1 - h), pML * Slf * (1 - h), Slf * (1 - h)),nrow = 3, byrow = T)
        
        if(states == 1){mat1 <- Matrix_A}  else {mat1 <- Matrix_N} # call the right projection matrix
        
        A_all[,,t] <- mat1
        
      } # from the stored matrices, compute reproductive values stable stage distribution 
      
      
      ### 2. COMPUTE POPULATION GROWTH RATES 
      
      w <- matrix(0,rows,(trun + 1))
      wvec <- (matrix(1,rows,1))/rows
      w[,1] <- wvec
      
      growth <- matrix(0,trun,1) # storage matrix for realized growth rate
      
      for (j in 1:trun) {
        
        A <- A_all[,,j]
        wvec <- A %*% wvec
        growth[j] <- sum(wvec) # per time step population growth rate
        wvec <- wvec/growth[j]
        w[,(j+1)] <- wvec
        
      }
      
      w <- w[,-c(1:50)]
      growth <- growth[-c(1:50)]
      
      lambda_t[z,y,,rep] = growth
      
    }
  }
}


# 3. COMPUTING POPULATION SIZE AT EACH TIME STEP

pop_t <- array(0,c(length(frequency), length(autocorrelation), time, rep.simul)) # store nb. of individuals at each time step
pop_t[,,1,] = 600 # starting with a population of 600 individuals

for (z in 1:length(frequency)) {
   for (y in 1:length(autocorrelation)) {
     for (rep in 1:rep.simul) {
       for (t in 2:time) {

         pop_t[z,y,t,rep] = pop_t[z,y,t-1,rep] * lambda_t[z,y,t,rep]

       }
     }
   }
}


# 4. COMPUTE DOUBLING TIME

doubling_time <- array(0,c(length(frequency), length(autocorrelation), rep.simul))

for (i in 1:length(frequency)) {
  for (j in 1:length(autocorrelation)) {
    for (k in 1:rep.simul) {
      
      doubling_time[i,j,k] = log(2)/(log(mean(lambda_t[i,j,,k])))
      
    }
  }
}


