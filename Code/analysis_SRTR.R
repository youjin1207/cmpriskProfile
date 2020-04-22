library(MASS)
library(survival)
source("auxfunctions.R")

#### exact matching ####
matchingdat = original
# exact matching with three key covariates: age, diabetes, bloodtype
age.cate = ifelse(matchingdat$age + 50 < 40, 1, 
                  ifelse(matchingdat$age + 50 < 50, 2, 
                         ifelse(matchingdat$age + 50 < 60, 3, 4)))
diabetes.cate = ifelse(matchingdat$diabetes == 0, 1, 2)
bloodtype.cate = ifelse(matchingdat$blood_a == 1, 1, 
                        ifelse(matchingdat$blood_b == 1, 2, 
                               ifelse(matchingdat$blood_ab == 1, 3, 4)))
matchingdat$age.cate = age.cate
matchingdat$diabetes.cate = diabetes.cate
matchingdat$bloodtype.cate = bloodtype.cate
rownames(matchingdat) = as.character(c(1:nrow(matchingdat)))

# population level risk class probability
table.mat = ftable(matchingdat$age.cate, matchingdat$diabetes.cate, matchingdat$bloodtype.cate)
template.mat = as.matrix(table.mat) 
prop.mat = template.mat / nrow(matchingdat)

center.ind = names(table(original$center))
J = length(center.ind)
time.list = seq(0, 2000, 100)
delta1 = delta2 = S = F1 = F2 = matrix(NA, J, length(time.list))
F1.all = F2.all = rep(NA, length(time.list))
S[,1] = 1
matchingdat$ww = rep(NA, nrow(matchingdat))
for(j in 1:length(center.ind)){
  print(j)
  
  tmp.center = matchingdat[matchingdat$center == center.ind[j],]
  
  if(nrow(tmp.center) > 0){
    count1 = count2 = atrisk =  Lambda1 = Lambda2 = 
      tmp.S = tmp.F1 = tmp.F2 = rep(0, nrow(tmp.center))
    
    # center level risk class probability
    order.center = tmp.center[order(tmp.center$timeto),]
    tmp.table.mat =  ftable(order.center$age.cate, order.center$diabetes.cate, order.center$bloodtype.cate)
    result.mat = real.matching(tmp.table.mat, template.mat)
    
    weights = prop.mat/result.mat
    weights = ifelse(is.infinite(weights), 0, weights)
    weights = weights*nrow(order.center)/sum(prop.mat[result.mat > 0])
    order.center$ww = rep(NA, nrow(order.center))
    
    for(t in 1:nrow(order.center)){
      if(order.center$cause1[t] == 1){
        tmp.single.table = ftable(order.center$age.cate[t], order.center$diabetes.cate[t], order.center$bloodtype.cate[t])
        n1.mat = real.matching(tmp.single.table, template.mat)
      }else{
        n1.mat = 0
      }
      
      if(order.center$cause2[t] == 1){
        tmp.single.table = ftable(order.center$age.cate[t], order.center$diabetes.cate[t], order.center$bloodtype.cate[t])
        n2.mat = real.matching(tmp.single.table, template.mat)
      }else{
        n2.mat = 0
      }
      count1[t] = sum(weights*n1.mat)
      count2[t] = sum(weights*n2.mat)
      
      y.index = which((order.center$timeto >= order.center$timeto[t]))
      tmp.riskset.table = ftable(order.center$age.cate[y.index], order.center$diabetes.cate[y.index], order.center$bloodtype.cate[y.index])
      y.mat = real.matching(tmp.riskset.table, template.mat)    
      atrisk[t] = sum(weights*y.mat)
      
      tmp.single.table = ftable(order.center$age.cate[t], order.center$diabetes.cate[t], order.center$bloodtype.cate[t])
      order.center$ww[t] = sum(weights*real.matching(tmp.single.table, template.mat))
      
      if(t==1){
        Lambda1[t] = 0 + (count1[t]/atrisk[t])
        Lambda2[t] = 0  + (count2[t]/atrisk[t])
      }else{
        Lambda1[t] = Lambda1[(t-1)] + (count1[t]/atrisk[t])
        Lambda2[t] = Lambda2[(t-1)] + (count2[t]/atrisk[t])
      }
      
      
      tmp.S[t] = exp(-Lambda1[t]-Lambda2[t])
      if(t==1){
        tmp.F1[t] = 0 + 1*(count1[t]/atrisk[t])
        tmp.F2[t] = 0 + 1*(count2[t]/atrisk[t])
      }else{
        tmp.F1[t] = tmp.F1[(t-1)] + tmp.S[(t-1)]*(count1[t]/atrisk[t])
        tmp.F2[t] = tmp.F2[(t-1)] + tmp.S[(t-1)]*(count2[t]/atrisk[t])
      }
      matchingdat$ww[as.integer(rownames(order.center))] = order.center$ww  
    }
    
    for(k in 1:length(time.list)){
      S[j,k] = tmp.S[which.min(abs(time.list[k] - order.center$timeto))]
      F1[j,k] = tmp.F1[which.min(abs(time.list[k] - order.center$timeto))]
      F2[j,k] = tmp.F2[which.min(abs(time.list[k] - order.center$timeto))]
    }
  }
}

for(t in 1:length(time.list)){
  F1.all[t] = mean(F1[,t], na.rm = TRUE)
  F2.all[t] = mean(F2[,t], na.rm = TRUE)
  for(j in 1:J){
    delta1[j,t] = F1[j,t] - F1.all[t]
    delta2[j,t] = F2[j,t] - F2.all[t]
  }
}
# variance estimation
matchingdat$obs.Y = matchingdat$timeto
matchingdat$Delta = ifelse(matchingdat$cause1 == 1, 1, ifelse(matchingdat$cause2 == 1, 2, 0))
matchingdat$Censor = matchingdat$censor
F1var.if = F2var.if = rep(NA, length(center.ind))
for(j in 1:length(center.ind)){
  print(j)
  tmp.center = matchingdat[matchingdat$center == center.ind[j],]
  if(nrow(tmp.center) > 0){
    if.result = IF.cif.weight(tmp.center, 1800, 2, weight = TRUE)
    F1var.if[j] = sum((if.result[,1])^2)/(nrow(tmp.center))^2
    F2var.if[j] = sum((if.result[,2])^2)/(nrow(tmp.center))^2
  }
}
# save the result
matching_onetime = list(F1 = F1, F2 = F2, delta1 = delta1, delta2 = delta2, 
                        F1var.if = F1var.if, F2var.if = F2var.if)

#### multivariate prognostic scores ####
scoredat = original
rownames(scoredat) = as.character(c(1:nrow(scoredat)))
# fit a cause-specific hazards
fit1 = coxph(Surv(timeto, cause1) ~  age + bmi + female +     
               race_black + race_hispanic + race_asian +  
               race_other + blood_a + blood_ab +     
               blood_b + diag_poly + diag_diab +     
               diag_hyper + diag_other + copd +        
               hyperten + malig + pvd +       
               diabetes + strata(center), 
              data = scoredat)
fit2 = coxph(Surv(timeto, cause2) ~  age + bmi + female +     
               race_black + race_hispanic + race_asian +  
               race_other + blood_a + blood_ab +     
               blood_b + diag_poly + diag_diab +     
               diag_hyper + diag_other + copd +        
               hyperten + malig + pvd +       
               diabetes + strata(center), 
             data = scoredat)

first.scores.whole = predict(fit1, reference = "sample")
second.scores.whole = predict(fit2, reference = "sample")
scoredat$first.scores.whole = first.scores.whole
scoredat$second.scores.whole = second.scores.whole
first.scores.raw = c(quantile(scoredat$first.scores.whole, seq(0.2, 1, 0.2))[1:4], max(scoredat$first.scores.whole) + 0.1)
second.scores.raw = c(quantile(scoredat$second.scores.whole, seq(0.2, 1, 0.2))[1:4], max(scoredat$second.scores.whole) + 0.1)  

# population-level risk class probability
dist.mat = matrix(0, nrow(scoredat), 2)
score1 = scoredat$first.scores.whole
score2 = scoredat$second.scores.whole
dist.mat[,1] = apply(as.matrix(score1), 1, function(x) sum(x > first.scores.raw))
dist.mat[,2] = apply(as.matrix(score2), 1, function(x) sum(x > second.scores.raw))
table.raw = table(dist.mat[,1], dist.mat[,2])
table.raw = as.matrix(table.raw) / nrow(scoredat)

center.ind = names(table(original$center))
J = length(center.ind)
time.list = seq(0, 2000, 100)
delta1 = delta2 = S = F1 = F2 = matrix(NA, J, length(time.list))
F1.all = F2.all = rep(NA, length(time.list))
S[,1] = 1
scoredat$ww = rep(NA, nrow(scoredat))
for(j in 1:length(center.ind)){
  print(j)
  
  tmp.center = scoredat[scoredat$center == center.ind[j],]
  
  if(nrow(tmp.center) > 0){
    count1 = count2 = atrisk =  Lambda1 = Lambda2 = 
      tmp.S = tmp.F1 = tmp.F2 = rep(0, nrow(tmp.center))
    
    order.center = tmp.center[order(tmp.center$timeto),]
    tmp.dist.mat = matrix(0, nrow(order.center), 2)
    score1 = scoredat$first.scores.whole[scoredat$center==center.ind[j]]
    score2 = scoredat$second.scores.whole[scoredat$center==center.ind[j]]
    
    # center level risk class probability
    tmp.dist.mat = matrix(0, nrow(tmp.center), 2)
    tmp.dist.mat[,1] = apply(as.matrix(score1), 1, function(x) sum(x > first.scores.raw))
    tmp.dist.mat[,2] = apply(as.matrix(score2), 1, function(x) sum(x > second.scores.raw))
    result.mat = print.55table(table(tmp.dist.mat[,1], tmp.dist.mat[,2]))
    weights = table.raw/result.mat
    weights = ifelse(is.infinite(weights), 0, weights)
    weights = weights*nrow(order.center)/sum(table.raw[result.mat > 0])
    
    order.center$ww = rep(NA, nrow(order.center))
    
    for(t in 1:nrow(order.center)){
      if(order.center$cause1[t] == 1){
        n1.mat = print.55table(table(tmp.dist.mat[t,1], tmp.dist.mat[t, 2]))
      }else{
        n1.mat = 0
      }
      
      if(order.center$cause2[t] == 1){
        n2.mat = print.55table(table(tmp.dist.mat[t,1], tmp.dist.mat[t, 2]))
      }else{
        n2.mat = 0
      }
      count1[t] = sum(weights*n1.mat)
      count2[t] = sum(weights*n2.mat)
      
      y.index = which((order.center$timeto >= order.center$timeto[t]))
      y.mat = print.55table(table(tmp.dist.mat[y.index,1], tmp.dist.mat[y.index, 2]))      
      atrisk[t] = sum(weights*y.mat)
      
      order.center$ww[t] = sum(weights*print.55table(table(tmp.dist.mat[t,1], tmp.dist.mat[t, 2])))
      
      if(t==1){
        Lambda1[t] = 0 + (count1[t]/atrisk[t])
        Lambda2[t] = 0  + (count2[t]/atrisk[t])
      }else{
        Lambda1[t] = Lambda1[(t-1)] + (count1[t]/atrisk[t])
        Lambda2[t] = Lambda2[(t-1)] + (count2[t]/atrisk[t])
      }
      
      
      tmp.S[t] = exp(-Lambda1[t]-Lambda2[t])
      if(t==1){
        tmp.F1[t] = 0 + 1*(count1[t]/atrisk[t])
        tmp.F2[t] = 0 + 1*(count2[t]/atrisk[t])
      }else{
        tmp.F1[t] = tmp.F1[(t-1)] + tmp.S[(t-1)]*(count1[t]/atrisk[t])
        tmp.F2[t] = tmp.F2[(t-1)] + tmp.S[(t-1)]*(count2[t]/atrisk[t])
      }
      scoredat$ww[as.integer(rownames(order.center))] = order.center$ww  
    }
    
    for(k in 1:length(time.list)){
      S[j,k] = tmp.S[which.min(abs(time.list[k] - order.center$timeto))]
      F1[j,k] = tmp.F1[which.min(abs(time.list[k] - order.center$timeto))]
      F2[j,k] = tmp.F2[which.min(abs(time.list[k] - order.center$timeto))]
    }
  }
}

# calculate the excess CIF 
for(t in 1:length(time.list)){
  F1.all[t] = mean(F1[,t], na.rm = TRUE)
  F2.all[t] = mean(F2[,t], na.rm = TRUE)
  for(j in 1:J){
    delta1[j,t] = F1[j,t] - F1.all[t]
    delta2[j,t] = F2[j,t] - F2.all[t]
  }
}

# variance estimation
scoredat$obs.Y = scoredat$timeto
scoredat$Delta = ifelse(scoredat$cause1 == 1, 1, ifelse(scoredat$cause2 == 1, 2, 0))
scoredat$Censor = scoredat$censor
F1var.if = F2var.if = rep(NA, length(center.ind))
for(j in 1:length(center.ind)){
  tmp.center = scoredat[scoredat$center == center.ind[j],]
  if(nrow(tmp.center) > 0){
    if.result = IF.cif.weight(tmp.center, 1800, 2, weight = TRUE)
    F1var.if[j] = sum((if.result[,1])^2)/(nrow(tmp.center))^2
    F2var.if[j] = sum((if.result[,2])^2)/(nrow(tmp.center))^2
  }
}
# save the result
score_onetime = list(F1 = F1, F2 = F2, delta1 = delta1, delta2 = delta2,  
                     F1var.if = F1var.if, F2var.if = F2var.if)


### naive method (without any weighting) ###
noweightdat = original
rownames(noweightdat) = as.character(c(1:nrow(noweightdat)))

center.ind = names(table(original$center))
J = length(center.ind)
time.list = seq(0, 2000, 100)
delta1 = delta2 = S = F1 = F2 = matrix(NA, J, length(time.list))
F1.all = F2.all = rep(NA, length(time.list))
S[,1] = 1
for(j in 1:length(center.ind)){
  print(j)
  tmp.center = noweightdat[noweightdat$center == center.ind[j],]
  order.center = tmp.center[order(tmp.center$timeto),]
  
  if(nrow(tmp.center) > 0){
    count1 = count2 = atrisk =  Lambda1 = Lambda2 = 
      tmp.S = tmp.F1 = tmp.F2 = rep(0, nrow(tmp.center))
    order.center$ww = rep(NA, nrow(order.center))
    for(t in 1:nrow(order.center)){
      if(order.center$cause1[t] == 1){
        n1.mat = 1
      }else{
        n1.mat = 0
      }
      
      if(order.center$cause2[t] == 1){
        n2.mat = 1
      }else{
        n2.mat = 0
      }
      count1[t] = sum(weights*n1.mat)
      count2[t] = sum(weights*n2.mat)
      
      y.index = which((order.center$timeto >= order.center$timeto[t]))
      atrisk[t] = length(y.index)
      
      if(t==1){
        Lambda1[t] = 0 + (count1[t]/atrisk[t])
        Lambda2[t] = 0  + (count2[t]/atrisk[t])
      }else{
        Lambda1[t] = Lambda1[(t-1)] + (count1[t]/atrisk[t])
        Lambda2[t] = Lambda2[(t-1)] + (count2[t]/atrisk[t])
      }
      
      
      tmp.S[t] = exp(-Lambda1[t]-Lambda2[t])
      if(t==1){
        tmp.F1[t] = 0 + 1*(count1[t]/atrisk[t])
        tmp.F2[t] = 0 + 1*(count2[t]/atrisk[t])
      }else{
        tmp.F1[t] = tmp.F1[(t-1)] + tmp.S[(t-1)]*(count1[t]/atrisk[t])
        tmp.F2[t] = tmp.F2[(t-1)] + tmp.S[(t-1)]*(count2[t]/atrisk[t])
      }
    }
    
    for(k in 1:length(time.list)){
      S[j,k] = tmp.S[which.min(abs(time.list[k] - order.center$timeto))]
      F1[j,k] = tmp.F1[which.min(abs(time.list[k] - order.center$timeto))]
      F2[j,k] = tmp.F2[which.min(abs(time.list[k] - order.center$timeto))]
    }
  }
}

# excess CIF
for(t in 1:length(time.list)){
  F1.all[t] = mean(F1[,t], na.rm = TRUE)
  F2.all[t] = mean(F2[,t], na.rm = TRUE)
  for(j in 1:J){
    delta1[j,t] = F1[j,t] - F1.all[t]
    delta2[j,t] = F2[j,t] - F2.all[t]
  }
}

# variance estimation
noweightdat$ww = rep(1, nrow(noweightdat))
noweightdat$obs.Y = noweightdat$timeto
noweightdat$Delta = ifelse(noweightdat$cause1 == 1, 1, ifelse(noweightdat$cause2 == 1, 2, 0))
noweightdat$Censor = noweightdat$censor
F1var.if = F2var.if = rep(NA, length(center.ind))
for(j in 1:length(center.ind)){
  print(j)
  tmp.center = noweightdat[noweightdat$center == center.ind[j],]
  if(nrow(tmp.center) > 0){
    result =  var.function.weight(tmp.center, 1800, 2, weight = TRUE)
    if.result = IF.cif.weight(tmp.center, 1800, 2, weight = TRUE)
    F1var.est[j] = result[1]
    F2var.est[j] = result[2]
    F1var.if[j] = sum((if.result[,1])^2)/(nrow(tmp.center))^2
    F2var.if[j] = sum((if.result[,2])^2)/(nrow(tmp.center))^2
  }
}
# save the result
noweight_onetime = list(F1 = F1, F2 = F2, delta1 = delta1, delta2 = delta2, F1var.est = F1var.est,
                        F2var.est = F2var.est, F1var.if = F1var.if, F2var.if = F2var.if)


