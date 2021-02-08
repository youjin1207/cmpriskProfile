library(MASS)
library(survival)
source("auxfunctions.R")

J = 100 # the number of centers 
nJ = 200 # the number of subjects per each center
N = J*nJ # total number of subjects in the simulated data

## multivariate normal
set.seed(1234)
U = mvrnorm(100, mu = c(0,0), Sigma = matrix(c(1, -0.2, -0.2, 1) , 2, 2))
phi1 = 0.00 + 0.03*plogis(U[,1])
phi2 = 0.10 + 0.10*plogis(U[,2])


# generate covariates under scenario (iii)
X1 = X2 = X3 = X4 = X5 = matrix(0, nrow = J, ncol = nJ)
for(j in 1:J){
  #prob = rep(0.1, 10) # scenario (i)
  #prob = c(rep(0.1-0.001*(j-1), 3), rep(0.1, 4), rep(0.1+0.001*(j-1), 3)) # scenario (ii)
  prob = c(rep(0.1-0.001*(j-1), 5), rep(0.1+0.001*(j-1), 5))
  X1[j,] = apply(rmultinom(nJ, 1, prob = prob), 2, function(x) which(x == 1)) / 10
  X2[j,] = apply(rmultinom(nJ, 1, prob = prob), 2, function(x) which(x == 1)) / 10
  X3[j,] = apply(rmultinom(nJ, 1, prob = prob), 2, function(x) which(x == 1)) / 10
  X4[j,] = apply(rmultinom(nJ, 1, prob = prob), 2, function(x) which(x == 1)) / 10
  X5[j,] = apply(rmultinom(nJ, 1, prob = prob), 2, function(x) which(x == 1)) / 10
}

T0 = Delta = C0 = Censor = Y0 = matrix(0, nrow = J, ncol = nJ);
beta1 = c(1, -1)*2
beta2 = c(1, -1)*2
beta3 = c(0.5, -1)*2
beta4 = c(-1, 0.5)*2
beta5 = c(0.5, 0.5)*2
gammas = rep(-0.5, 5)*2
tau = 3/2 # shape parameter
# 
for(j in 1:J){
  phis[1] = phi1[j]
  phis[2] = phi2[j]
  W1 = runif(nJ, 0, 1)
  T0[j,] = (-log(W1) / (exp(beta1[1]*X1[j,] + beta2[1]*X2[j,] + 
                              beta3[1]*X3[j,] + beta4[1]*X4[j,] + beta5[1]*X5[j,])*(phis[1])^(tau) + 
                          exp(beta1[2]*X1[j,] + beta2[2]*X2[j,] + 
                                beta3[2]*X3[j,] + beta4[2]*X4[j,] + beta5[2]*X5[j,])*(phis[2])^(tau)))^(1 / tau)
  csh1 = (T0[j,])^(tau - 1)*exp(beta1[1]*X1[j,] + beta2[1]*X2[j,] + 
                                  beta3[1]*X3[j,] + beta4[1]*X4[j,] + beta5[1]*X5[j,])*(phis[1]^(tau))*tau
  csh2 = (T0[j,])^(tau - 1)*exp(beta1[2]*X1[j,] + beta2[2]*X2[j,] + 
                                  beta3[2]*X3[j,] + beta4[2]*X4[j,] + beta5[2]*X5[j,])*(phis[2]^(tau))*tau
  Delta[j,] = rbinom(nJ, 1, csh1 / (csh1 + csh2))
  
  W2 = runif(nJ, 0, 1)
  C0[j,] = (-log(W2) / (0.5*exp(gammas[1]*X1[j,] + gammas[2]*X2[j,] + 
                                  gammas[3]*X3[j,] + gammas[4]*X4[j,] + gammas[5]*X5[j,])))
  
  Censor[j,] = T0[j,] <= C0[j,]
  
  Y0[j,] = pmin(T0[j,], C0[j,])
}

sim1 = data.frame(timeto = as.numeric(t(Y0)), X1 = as.numeric(t(X1)), 
                  X2 = as.numeric(t(X2)), X3 = as.numeric(t(X3)),
                  X4 = as.numeric(t(X4)), X5 = as.numeric(t(X5)),
                  center = rep(1:J, each = nJ),
                  id = 1:N, Delta = as.numeric(t(Delta)),
                  Censor = as.numeric(t(Censor)))
# misspecified covariates W
sim1$W1 = exp(sim1$X1/2)
sim1$W2 = sim1$X2 / (1+exp(sim1$X1)) + 10
sim1$W3 = (sim1$X1*sim1$X3/25 + 0.6 )^3
sim1$W4 = (sim1$X2 + sim1$X4 + 20)^2
sim1$W5 = sim1$X5

# fit cause-specific hazards to derive prognostic scores for each cause
sim1$cause1 = sim1$Delta == 1 & sim1$Censor == 1
sim1$cause2 = sim1$Delta == 0 & sim1$Censor == 1
fit1 = coxph( Surv(timeto, cause1) ~ X1 + X2 + X3 + X4 + X5 + strata(center), data = sim1) # correctly specified case
fit2 = coxph( Surv(timeto, cause2) ~ X1 + X2 + X3 + X4 + X5 + strata(center), data = sim1) # correctly specified case
#fit1 = coxph( Surv(timeto, cause1) ~ W1 + W2 + W3 + W4 + W5 + strata(center), data = sim1) # misspecified covariates case
#fit2 = coxph( Surv(timeto, cause2) ~ W1 + W2 + W3 + W4 + W5 + strata(center), data = sim1) # misspecified covariates case
#first.scores.whole = additive.score.sim.competing(dat = sim1, cause = 1) # misspecified prognostic score model
#second.scores.whole = additive.score.sim.competing(dat = sim1, cause = 2) # misspecified prognostic score model
first.scores.whole = predict(fit1, reference = "sample")
second.scores.whole = predict(fit2, reference = "sample")
sim1$first.scores.whole = first.scores.whole
sim1$second.scores.whole = second.scores.whole
# thresholds by quintiles 
first.scores = c(quantile(first.scores.whole , seq(0.2, 1, 0.2))[1:4], max(first.scores.whole) + 0.1)
second.scores = c(quantile(second.scores.whole, seq(0.2, 1, 0.2))[1:4], max(second.scores.whole) + 0.1)

# population level risk class probability
dist.mat = matrix(0, nrow(sim1), 2)
dist.mat[,1] = apply(as.matrix(sim1$first.scores.whole), 1, function(x) sum(x > first.scores))
dist.mat[,2] = apply(as.matrix(sim1$second.scores.whole), 1, function(x) sum(x > second.scores))
table.mat = table(dist.mat[,1], dist.mat[,2])
table.mat = as.matrix(table.mat) / nrow(sim1)

time.list = seq(0, 30, 1)
count1 = count2 = atrisk =  Lambda1 = Lambda2 = tmp.S = tmp.F1 = tmp.F2 = matrix(0, J, nJ)
delta1 = delta2 = S = F1 = F2 = matrix(0, J, length(time.list))
F1.all = F2.all = rep(0, length(time.list))
S[,1] = 1
sim1$ww = rep(NA, nrow(sim1)) # individual-level weights
# iterate across J centers:
for(j in 1:J){
  print(j)
  tmp.center = sim1[sim1$center == j,]
  order.center = tmp.center[order(tmp.center$timeto),]
  tmp.dist.mat = matrix(0, nrow(order.center), 2)
  tmp.dist.mat = matrix(0, nrow(order.center), 2)
  tmp.dist.mat[,1] = apply(as.matrix(order.center$first.scores.whole), 1, function(x) sum(x > first.scores))
  tmp.dist.mat[,2] = apply(as.matrix(order.center$second.scores.whole), 1, function(x) sum(x > second.scores))
  
  # center-level risk class probability
  result.mat = print.55table(table(tmp.dist.mat[,1], tmp.dist.mat[,2]))
  weights = table.mat/result.mat
  weights = ifelse(is.infinite(weights), 0, weights)
  weights = ifelse(is.na(weights), 0, weights)
  weights = weights*nJ/sum(table.mat[result.mat > 0])
  order.center$ww = rep(NA, nrow(order.center))
  
  # iterate across observed failure times within center:
  for(t in 1:nrow(order.center)){
    if(order.center$cause1[t] == 1 & order.center$Censor[t] == 1){
      n1.mat = print.55table(table(tmp.dist.mat[t,1], tmp.dist.mat[t, 2]))
    }else{
      n1.mat = 0
    }
    
    if(order.center$cause2[t] == 1 & order.center$Censor[t] == 1){
      n2.mat = print.55table(table(tmp.dist.mat[t,1], tmp.dist.mat[t, 2]))
    }else{
      n2.mat = 0
    }
    count1[j,t] = sum(weights*n1.mat)
    count2[j,t] = sum(weights*n2.mat)
    
    y.index = which((order.center$timeto >= order.center$timeto[t]))
    y.mat = print.55table(table(tmp.dist.mat[y.index,1], tmp.dist.mat[y.index, 2]))      
    atrisk[j,t] = sum(weights*y.mat)
    
    order.center$ww[t] = sum(weights* print.55table(table(tmp.dist.mat[t,1], tmp.dist.mat[t, 2])))
    
    # cumulative cause-specific hazards
    if(t==1){
      Lambda1[j,t] = 0 + (count1[j,t]/atrisk[j,t])
      Lambda2[j,t] = 0  + (count2[j,t]/atrisk[j,t])
    }else{
      Lambda1[j,t] = Lambda1[j,(t-1)] + (count1[j,t]/atrisk[j,t])
      Lambda2[j,t] = Lambda2[j,(t-1)] + (count2[j,t]/atrisk[j,t])
    }
    
    # survival function
    tmp.S[j,t] = exp(-Lambda1[j,t]-Lambda2[j,t])
    tmp.S[j,t] = exp(-Lambda1[j,t]-Lambda2[j,t])
    # CIF for each cause
    if(t==1){
      tmp.F1[j,t] = 0 + ((1 +tmp.S[j,t])/2)*(count1[j,t]/atrisk[j,t])
      tmp.F2[j,t] = 0 + 1*(count2[j,t]/atrisk[j,t])
    }else{
      tmp.F1[j,t] = tmp.F1[j,(t-1)] + ((tmp.S[j,(t-1)] +tmp.S[j,t])/2)*(count1[j,t]/atrisk[j,t])
      tmp.F2[j,t] = tmp.F2[j,(t-1)] + ((tmp.S[j,(t-1)] +tmp.S[j,t])/2)*(count2[j,t]/atrisk[j,t])
    }
    sim1$ww[as.integer(rownames(order.center))] = order.center$ww # save individual-level weights
  }
  # approximate the time points to specified times points at time.list 
  for(k in 1:length(time.list)){
    S[j,k] = tmp.S[j,which.min(abs(time.list[k] - order.center$timeto))]
    F1[j,k] = tmp.F1[j,which.min(abs(time.list[k] - order.center$timeto))]
    F2[j,k] = tmp.F2[j,which.min(abs(time.list[k] - order.center$timeto))]
  }
}

# derive excess CIF across time points
for(t in 1:length(time.list)){
  F1.all[t] = mean(F1[,t])
  F2.all[t] = mean(F2[,t])
  for(j in 1:J){
    delta1[j,t] = F1[j,t] - F1.all[t]
    delta2[j,t] = F2[j,t] - F2.all[t]
  }
}

# variance estimation at t=20
sim1$cause1 = sim1$Delta == 1 & sim1$Censor == 1
sim1$cause2 = sim1$Delta == 0 & sim1$Censor == 1
sim1$obs.Y = sim1$timeto
sim1$center = sim1$center
sim1$Delta = ifelse(sim1$cause1 == 1, 1, ifelse(sim1$cause2 == 1, 2, 0))
dat = sim1
F1var.if = F2var.if = c()
for(k in 1:100){
  if.result = IF.cif.weight(dat[dat$center == k, ], 20, 2, weight = TRUE)
  F1var.if[k] = sum((if.result[,1])^2)/(nrow(dat[dat$center == k, ]))^2 
  F2var.if[k] = sum((if.result[,2])^2)/(nrow(dat[dat$center == k, ]))^2
}

# results 
print(list(F1 = F1, F2 = F2,  # estimated CIF
           delta1 = delta1, delta2 = delta2,  # estimated center effect
           weights = sim1$ww,  # estimated individual-level weights 
          F1var.if = F1var.if, F2var.if = F2var.if  # estimated variance of the CIF estmates
          ))

