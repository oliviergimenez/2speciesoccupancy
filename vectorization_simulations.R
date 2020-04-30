library(mipfp)
library(unmarked)

# occupancy probabilities
psi01 <- 81/175
psi10 <- 36/175
psi11 <- 4/175
psi00 <- 1 - (psi01 + psi10 + psi11) # 54/175
psiS1 <- psi10 + psi11
psiS2 <- psi01 + psi11
or <- matrix(c(1, (psiS1*(1-psiS2))/(psiS2*(1-psiS1)), 
               (psiS2*(1-psiS1))/(psiS1*(1-psiS2)), 1), nrow = 2, ncol = 2, byrow = TRUE)
rownames(or) <- colnames(or) <- c("sp1", "sp2")
marg.probs <- c(psiS1, psiS2)
p.joint <- ObtainMultBinaryDist(odds = or, marg.probs = marg.probs)

# detection probabilities
ps <- c(0.5,0.9)

# function for simulating data - brute force
occu_sim <- function(
  N = 50, # number of sites
  J = 5, # number of surveys
  S = 2, # nb species 
  n_sim = 500){ # number of simulations
  
  # simulate n_sim datasets
  for (j in 1:n_sim){
    z <- RMultBinary(n = N, mult.bin.dist = p.joint)$binary.sequences 
    y <- list()
    for (i in 1:S){
      y[[i]] <- matrix(NA,N,J)
      for (j in 1:N){
        for (k in 1:J){
          y[[i]][j,k] <- rbinom(1,1,z[j,i]*ps[i])
        }
      }
    }
    names(y) <- c('sp1','sp2')
  }
}

# function for simulating data - with vectorization
occu_sim_vect <- function(
  N = 50, # number of sites
  J = 5, # number of surveys
  S = 2, # nb species 
  n_sim = 500){ # number of simulations
  
  # simulate n_sim datasets
  for (j in 1:n_sim){
    z <- RMultBinary(n = N, mult.bin.dist = p.joint)$binary.sequences 
    z_rep <- rep(z, each = J)
    ps_rep <- rep(ps, each = N*J)
    y <- array(rbinom(N*J*S,1,z_rep*ps_rep),
               dim = c(J,N,S))
    y <- aperm(y, c(2,1,3))
  }
}

# benchmark w/ 250 sites, 20 visits and 500 simulations
system.time(occu_sim(N = 250, J = 20, S = 2, n_sim = 500))
system.time(occu_sim_vect(N = 250, J = 20, S = 2, n_sim = 500))

