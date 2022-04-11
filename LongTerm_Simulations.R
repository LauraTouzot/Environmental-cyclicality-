######################################################################
#    Amplified cyclicality in mast seeding dynamics positively       #
#    influences the dynamics of a seed consumer species              #
######################################################################


### 1. LONG TERM SIMULATIONS


# Autocorrelation 
# 2 environmental states
# 3 x 3 matrix
# 4 frequencies
# 11 autocorrelation coefficients (including 0!)
# 1,000 simulations of 100,000 time steps
# 1,000 first iterations discarded


# EXAMPLE FOR AUTOCORRELATION COEFFICIENT v1 = -0.50
library(popbio)

rm(list=ls(all=TRUE))

frequency <- c(0.35,0.45, 0.55, 0.65) # frequency of good years
autocorrelation <- c(-0.50)          # coefficient of autocorrelation
trun = 600                           # run each simulation for trun years
env = c("good","bad")                       # the different environments 
env.n = 1:2
rep.simul = 2                           # number of repetition for each condition (frequency x autocorrelation)

rows = 3
cols = 3
nstage = 3

f = frequency
v1 = autocorrelation

env_at_t = 1                                # state of the environment


### Calculation of the stochastic growth rates and elasticities

lambda_stoch <- array(0,c(length(frequency), length(autocorrelation),rep.simul)) # store stochastic growth rate
elas_mean <- array(0,c(rows, cols, length(frequency), length(autocorrelation),rep.simul))    # store stochastic elasticity - mean
elas_var <- array(0,c(rows, cols, length(frequency), length(autocorrelation),rep.simul))     # store stochastic elasticity - variance
elas_all <- array(0,c(rows, cols, length(frequency), length(autocorrelation),rep.simul))     # store stochastic elasticity - mean + variance


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
        
      } # from the stored matrices, compute reproductive values stable stage distribution and elasticities
      
      
      ### 2. COMPUTE STOCHASTIC GROWTH RATES & STOCHASTIC ELASTICITIES
      
      v <- matrix(0,rows,(trun + 1))
      w <- matrix(0,rows,(trun + 1))
      
      vvec <- (matrix(1,1,cols))/rows
      v[,(trun + 1)] <- t(vvec)
      
      wvec <- (matrix(1,rows,1))/rows
      w[,1] <- wvec
      
      for (i in trun:1) {
        
        A <- A_all[,,i]   # call the corresponding projection matrix
        vvec <- vvec %*% A
        vvec <- vvec/sum(vvec)
        v[,i] <- t(vvec) # store the sequence of reproductive value vectors in reverse time
        
      } 
      
      v <- v[,-c(1:1000)]
      
      growth <- matrix(0,trun,1) # storage matrix for realized growth rate
      
      for (j in 1:trun) {
        
        A <- A_all[,,j]
        wvec <- A %*% wvec
        growth[j] <- sum(wvec) # per time step population growth rate
        wvec <- wvec/growth[j]
        w[,(j+1)] <- wvec
        
      }
      
      w <- w[,-c(1:1000)]
      
      growth <- growth[-c(1:1000)]
      A_all <- A_all[,,-c(1:1000)]
      
      A_mean <- matrix(0,rows,cols) # compute mean matrix
      
      for (k in 1:(trun - 1000)) {
        A_mean <- A_mean + A_all[,,k]
      }
      
      A_mean <- A_mean/(trun - 1000)
      
      elasticity_matrix_all <- matrix(0,rows,cols)
      elasticity_matrix_mean <- matrix(0,rows,cols)
      elasticity_matrix_var <- matrix(0,rows,cols)
      
      
      
      for (l in 1:(trun - 1000)) {
        
        A <- A_all[,,l]
        elasticity_matrix_all <- elasticity_matrix_all + (v[,(l+1)] %*% t(w[,l]) * A_all[,,l])/as.numeric(growth[l] * t(v[,(l+1)]) %*% w[,(l+1)])
        elasticity_matrix_mean <- elasticity_matrix_mean + (v[,(l+1)] %*% t(w[,l]) * A_mean)/as.numeric(growth[l] * t(v[,(l+1)]) %*% w[,(l+1)])
        elasticity_matrix_var <- elasticity_matrix_var + (v[,(l+1)] %*% t(w[,l]) * (A_all[,,l] - A_mean))/as.numeric(growth[l] * t(v[,(l+1)]) %*% w[,(l+1)])
        
      }
      
      
      lambda_stoch[z,y,rep] = mean(growth) 
      elas_all[,,z,y,rep] = round(elasticity_matrix_all/(trun - 1000),4)
      elas_mean[,,z,y,rep] = round(elasticity_matrix_mean/(trun - 1000),4)
      elas_var[,,z,y,rep] = round(elasticity_matrix_var/(trun - 1000),4)
    }
  }
}

