library(MASS)
library(survival)
source("auxfunctions.R")
scoredat = partdat # read the original data
rownames(scoredat) = as.character(c(1:nrow(scoredat)))
# fit a cause-specific hazards model to obtain prognostic scores
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

# center effect measure: excess CIF
for(t in 1:length(time.list)){
  F1.all[t] = mean(F1[,t], na.rm = TRUE)
  F2.all[t] = mean(F2[,t], na.rm = TRUE)
  for(j in 1:J){
    delta1[j,t] = F1[j,t] - F1.all[t]
    delta2[j,t] = F2[j,t] - F2.all[t]
  }
}


# variance estimation at ten different time points
scoredat$obs.Y = scoredat$timeto
scoredat$Delta = ifelse(scoredat$cause1 == 1, 1, ifelse(scoredat$cause2 == 1, 2, 0))
scoredat$Censor = scoredat$censor
F1var.est = F2var.est = F1var.if = F2var.if = matrix(NA, length(center.ind), 10)
F1var.est.uw = F2var.est.uw = F1var.if.uw = F2var.if.uw = matrix(NA, length(center.ind), 10)
tentimes = seq(200, 2000, 200)
for(j in 1:length(center.ind)){
  for(t in 1:length(time.list)){
    print(j)
    tmp.center = matchingdat[matchingdat$center == center.ind[j],]
    if(nrow(tmp.center) > 0){
      if.result = IF.cif.weight(tmp.center, tentimes[t], 2, weight = TRUE)
      F1var.if[j,t] = sum((if.result[,1])^2)/(nrow(tmp.center))^2
      F2var.if[j,t] = sum((if.result[,2])^2)/(nrow(tmp.center))^2
    }
  }
}


score_onetime_timevarying = list(F1 = F1, F2 = F2, delta1 = delta1, delta2 = delta2, F1var.est = F1var.est,
                                 F2var.est = F2var.est, F1var.if = F1var.if, F2var.if = F2var.if)


# read the estimates at ten time points
noweight.delta1.est = noweight_onetime_timevarying$delta1[,seq(3, 21, 2)]
noweight.delta2.est = noweight_onetime_timevarying$delta2[,seq(3, 21, 2)]
score.delta1.est = score_onetime_timevarying$delta1[,seq(3, 21, 2)]
score.delta2.est = score_onetime_timevarying$delta2[,seq(3, 21, 2)]
matching.delta1.est = matching_onetime_timevarying$delta1[,seq(3, 21, 2)]
matching.delta2.est = matching_onetime_timevarying$delta2[,seq(3, 21, 2)]

time.list = seq(200, 2000, 200)
noweight.delta1.low = noweight.delta1.high = score.delta1.low = score.delta1.high = 
  matching.delta1.low = matching.delta1.high = matrix(NA, nrow(matching.delta2.est),length(time.list))
noweight.delta2.low = noweight.delta2.high = score.delta2.low = score.delta2.high = 
  matching.delta2.low = matching.delta2.high = matrix(NA, nrow(matching.delta2.est),length(time.list))
# read the lower and upper bound of the estimates at ten time points
for(t in 1:length(time.list)){
  noweight.delta1.low[,t] = noweight.delta1.est[,t] - 1.96*sqrt((210/212)*noweight_onetime_timevarying$F1var.if[,t] + sum((1/212)^2*noweight_onetime_timevarying$F1var.if[!is.na(noweight_onetime$F1var.if), t]))
  noweight.delta1.high[,t] = noweight.delta1.est[,t] + 1.96*sqrt((210/212)*noweight_onetime_timevarying$F1var.if[,t] + sum((1/212)^2*noweight_onetime_timevarying$F1var.if[!is.na(noweight_onetime$F1var.if),t]))
  score.delta1.low[,t] = score.delta1.est[,t] - 1.96*sqrt((210/212)*score_onetime_timevarying$F1var.if[,t] + sum((1/212)^2*score_onetime_timevarying$F1var.if[!is.na(noweight_onetime$F1var.if),t]))
  score.delta1.high[,t] = score.delta1.est[,t] + 1.96*sqrt((210/212)*score_onetime_timevarying$F1var.if[,t] + sum((1/212)^2*score_onetime_timevarying$F1var.if[!is.na(noweight_onetime$F1var.if),t]))
  matching.delta1.low[,t] = matching.delta1.est[,t] - 1.96*sqrt((210/212)*matching_onetime_timevarying$F1var.if[,t] + sum((1/212)^2*matching_onetime_timevarying$F1var.if[!is.na(noweight_onetime$F1var.if),t]))
  matching.delta1.high[,t] = matching.delta1.est[,t] + 1.96*sqrt((210/212)*matching_onetime_timevarying$F1var.if[,t] + sum((1/212)^2*matching_onetime_timevarying$F1var.if[!is.na(noweight_onetime$F1var.if),t]))
  
  
  noweight.delta2.low[,t] = noweight.delta2.est[,t] - 1.96*sqrt((210/212)*noweight_onetime_timevarying$F2var.if[,t] + sum((1/212)^2*noweight_onetime_timevarying$F2var.if[!is.na(noweight_onetime$F2var.if),t]))
  noweight.delta2.high[,t] = noweight.delta2.est[,t] + 1.96*sqrt((210/212)*noweight_onetime_timevarying$F2var.if[,t] + sum((1/212)^2*noweight_onetime_timevarying$F2var.if[!is.na(noweight_onetime$F2var.if),t]))
  score.delta2.low[,t] = score.delta2.est[,t] - 1.96*sqrt((210/212)*score_onetime_timevarying$F2var.if[,t] + sum((1/212)^2*score_onetime_timevarying$F2var.if[!is.na(noweight_onetime$F2var.if),t]))
  score.delta2.high[,t] = score.delta2.est[,t] + 1.96*sqrt((210/212)*score_onetime_timevarying$F2var.if[,t] + sum((1/212)^2*score_onetime_timevarying$F2var.if[!is.na(noweight_onetime$F2var.if),t]))
  matching.delta2.low[,t] = matching.delta2.est[,t] - 1.96*sqrt((210/212)*matching_onetime_timevarying$F2var.if[,t] + sum((1/212)^2*matching_onetime_timevarying$F2var.if[!is.na(noweight_onetime$F2var.if),t]))
  matching.delta2.high[,t] = matching.delta2.est[,t] + 1.96*sqrt((210/212)*matching_onetime_timevarying$F2var.if[,t] + sum((1/212)^2*matching_onetime_timevarying$F2var.if[!is.na(noweight_onetime$F2var.if),t]))
}


time.list = seq(0, 2000, 100)
inds = order(score.delta2.est[,10])[c(12, 54, 107, 159, 202)] # (5th, 25th, 50th, 75th, 95th)-percentiles
ind1 = inds[1]; ind2 = inds[2]; ind3 = inds[3]; ind4 = inds[4]; ind5 = inds[5]
pdf("Figure/delta_nonrandom_five.pdf", width = 18, height = 6)
par(mfrow = c(1,2), cex.lab = 2, cex.axis = 1.5, cex.main = 2,
    mar=c(4,6,3,1), tcl = 0.5, oma = c(2, 2, 0, 3), xpd = FALSE)
plot(time.list[seq(3,21,2)], 
     score.delta1.est[ind1,1:10] , col = "dodgerblue", 
     cex = 1, ylim = c(-0.2, 0.2),  
     ylab = expression(paste("Estimated ", hat(delta)[j1], "(t)")), 
     xlab = "Time",
     main = "Cause 1 (Death)", type = "b", lwd = 2, pch = 20)
lines(time.list[seq(3,21,2)], score.delta1.low[ind1,1:10], col = "dodgerblue", lwd = 0.5, lty = 2, type = "l", pch = 20)
lines(time.list[seq(3,21,2)], score.delta1.high[ind1,1:10], col = "dodgerblue", lwd = 0.5, lty = 2, type = "l", pch = 20)

lines(time.list[seq(3,21,2)], score.delta1.est[ind2,1:10], col = "firebrick1", lwd = 2, type = "b", pch = 20)
lines(time.list[seq(3,21,2)], score.delta1.low[ind2,1:10] , col = "firebrick1", lwd = 0.5, lty = 2, type = "l", pch = 20)
lines(time.list[seq(3,21,2)], score.delta1.high[ind2,1:10] , col = "firebrick1", lwd = 0.5, lty = 2, type = "l", pch = 20)


lines(time.list[seq(3,21,2)], score.delta1.est[ind3,1:10] , col = "lightgreen", lwd = 2, type = "b", pch = 20)
lines(time.list[seq(3,21,2)], score.delta1.low[ind3,1:10] , col = "lightgreen", lwd = 0.5, lty = 2, type = "l", pch = 20)
lines(time.list[seq(3,21,2)], score.delta1.high[ind3,1:10] , col = "lightgreen", lwd = 0.5, lty = 2, type = "l", pch = 20)

lines(time.list[seq(3,21,2)], score.delta1.est[ind4,1:10] , col = "orange", lwd = 2, type = "b", pch = 20)
lines(time.list[seq(3,21,2)], score.delta1.low[ind4,1:10] , col = "orange", lwd = 0.5, lty = 2, type = "l", pch = 20)
lines(time.list[seq(3,21,2)], score.delta1.high[ind4,1:10] , col = "orange", lwd = 0.5, lty = 2, type = "l", pch = 20)


lines(time.list[seq(3,21,2)], score.delta1.est[ind5,1:10] , col = "purple", lwd = 2, type = "b", pch = 20)
lines(time.list[seq(3,21,2)], score.delta1.low[ind5,1:10] , col = "purple", lwd = 0.5, lty = 2, type = "l", pch = 20)
lines(time.list[seq(3,21,2)], score.delta1.high[ind5,1:10] , col = "purple", lwd = 0.5, lty = 2, type = "l", pch = 20)
legend("topleft", c("Center 64 (5th)", "Center 404 (25th)", "Center 813 (50th)", "Center 866 (75th)", "Center 499 (95th)"),
       col = c("dodgerblue", "firebrick1", "lightgreen", "orange", "purple"), bty = "n", lwd = 2,
       cex = 1.5)


abline(h = 0, col = "grey")
plot(time.list[seq(3,21,2)], 
     score.delta2.est[ind1,1:10], col = "dodgerblue", 
     cex = 1, ylim = c(-0.4, 0.4),  
     ylab = expression(paste("Estimated ", hat(delta)[j2], "(t)")), 
     xlab = "Time",
     main = "Cause 2 (Transplant)", type = "b", lwd = 2, pch = 20)
lines(time.list[seq(3,21,2)], score.delta2.low[ind1,1:10], col = "dodgerblue", lwd = 0.5, lty = 2, type = "l", pch = 20)
lines(time.list[seq(3,21,2)], score.delta2.high[ind1,1:10], col = "dodgerblue", lwd = 0.5, lty = 2, type = "l", pch = 20)

lines(time.list[seq(3,21,2)], score.delta2.est[ind2,1:10], col = "firebrick1", lwd = 2, type = "b", pch = 20)
lines(time.list[seq(3,21,2)], score.delta2.low[ind2,1:10], col = "firebrick1", lwd = 0.5, lty = 2, type = "l", pch = 20)
lines(time.list[seq(3,21,2)], score.delta2.high[ind2,1:10], col = "firebrick1", lwd = 0.5, lty = 2, type = "l", pch = 20)

lines(time.list[seq(3,21,2)], score.delta2.est[ind3,1:10], col = "lightgreen", lwd = 2, type = "b", pch = 20)
lines(time.list[seq(3,21,2)], score.delta2.low[ind3,1:10], col = "lightgreen", lwd = 0.5, lty = 2, type = "l", pch = 20)
lines(time.list[seq(3,21,2)], score.delta2.high[ind3,1:10], col = "lightgreen", lwd = 0.5, lty = 2, type = "l", pch = 20)

lines(time.list[seq(3,21,2)], score.delta2.est[ind4,1:10], col = "orange", lwd = 2, type = "b", pch = 20)
lines(time.list[seq(3,21,2)], score.delta2.low[ind4,1:10], col = "orange", lwd = 0.5, lty = 2, type = "l", pch = 20)
lines(time.list[seq(3,21,2)], score.delta2.high[ind4,1:10], col = "orange", lwd = 0.5, lty = 2, type = "l", pch = 20)

lines(time.list[seq(3,21,2)], score.delta2.est[ind5,1:10], col = "purple", lwd = 2, type = "b", pch = 20)
lines(time.list[seq(3,21,2)], score.delta2.low[ind5,1:10], col = "purple", lwd = 0.5, lty = 2, type = "l", pch = 20)
lines(time.list[seq(3,21,2)], score.delta2.high[ind5,1:10], col = "purple", lwd = 0.5, lty = 2, type = "l", pch = 20)

abline(h = 0, col = "grey")
dev.off()
