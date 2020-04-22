## a function to print out 5 x 5 risk classes
print.55table = function(mat){
  result.mat = matrix(0, 5, 5); 
  for(i in 1:5){
    for(j in 1:5){
      q = mat[rownames(mat) == as.character(i-1),
              colnames(mat) == as.character(j-1)]
      if(length(q) > 0){
        result.mat[i,j] = mat[rownames(mat) == as.character(i-1),
                              colnames(mat) == as.character(j-1)]
      }
    }
  }
  return(result.mat)
}

## a function to match the observed risk classes to the population (template) risk classes
real.matching = function(mat, template.mat){
  mat = as.matrix(mat)
  result.mat = matrix(0, nrow = nrow(template.mat), ncol = ncol(template.mat))
  rownames(result.mat) = rownames(template.mat); colnames(result.mat) = colnames(template.mat)
  
  mat2 = result.mat
  if(ncol(mat) < ncol(template.mat)){
    for(j in 1:ncol(template.mat)){
      if(colnames(template.mat)[j] %in% colnames(mat)){
        mat2[which(rownames(mat2) %in% rownames(mat) ) ,j] = mat[,which(colnames(mat) %in% colnames(template.mat)[j])]
      }else{
        mat2[,j] = rep(0, nrow(mat2))
      }
    }
  }else{
    mat2 = mat
  }
  
  for(i in 1:nrow(result.mat)){
    if(rownames(result.mat)[i] %in% rownames(mat2)){
      result.mat[i,] = mat2[which(rownames(mat2) %in% rownames(template.mat)[i]),]
    }
  }
  return(result.mat)
}

## a function to return the influence function-based variance estimate for (weighted) CIF estimator
IF.cif.weight = function(dat, time.point, K, weight = FALSE){
  # dat : obs.Y, center, Delta
  # K : the number of competing risks
  if(weight == TRUE){
    ww = dat$ww
    km_fit = survfit(Surv(timeto, Censor) ~ 1, data = dat, weights = ww)
  }else{
    ww = rep(1, nrow(dat))
    km_fit = survfit(Surv(timeto, Censor) ~ 1, data = dat, weights = rep(1, nrow(dat)))
  }
  
  Shat = km_fit$surv
  time.list = km_fit$time
  #dat = dat[order(dat$obs.Y),]
  # define aux, dM 
  aux = Lambda_l = matrix(0, length(time.list), K)
  risk.mean = Lambda = rep(0, length(time.list))
  for(t in 1:length(time.list)){
    risk.mean[t] = mean(ww*(dat$obs.Y >= time.list[t]))
    if(t == 1){
      Lambda[t] = sum(ww*(dat$obs.Y == time.list[t] & dat$Censor == 1))/
        sum(ww*(dat$obs.Y >= time.list[t]))
    }else{
      Lambda[t] = Lambda[t-1] + sum(ww*(dat$obs.Y == time.list[t] & dat$Censor == 1))/
        sum(ww*(dat$obs.Y >= time.list[t]))
    }
    #######################
    ## for each center l ##
    for(l in 1:K){
      if(t == 1){
        Lambda_l[t,l] =  sum(ww*(dat$obs.Y == time.list[t] & dat$Censor == 1 & dat$Delta == l))/
          sum(ww*(dat$obs.Y >= time.list[t]))
      }else{
        Lambda_l[t,l] =  Lambda_l[t-1,l] + sum(ww*(dat$obs.Y == time.list[t] & dat$Censor == 1 & dat$Delta == l))/
          sum(ww*(dat$obs.Y >= time.list[t]))
      }
    }
  }
  time.m = which.min(abs(time.point - time.list))
  ## for cause k
  results = matrix(NA, nrow(dat), K)
  aux = matrix(0, nrow(dat), time.m)
  for(l in 1:K){
    A = B =  0
    for(t in 1:time.m){
      if(t == 1){
        aux[,t] = (ww*(dat$obs.Y == time.list[t]) - ww*(dat$obs.Y >= time.list[t])*(Lambda[t] - 0))/(risk.mean[t])
        A = A + 1*(ww*(dat$obs.Y == time.list[t] & dat$Delta == l) - ww*(dat$obs.Y >= time.list[t])*(Lambda_l[t,l] - 0))/(risk.mean[t])
        B = B + 1*0*(Lambda_l[t,l]-0)
      }else{
        aux[,t] = aux[,t-1] + (ww*(dat$obs.Y == time.list[t]) - ww*(dat$obs.Y >= time.list[t])*(Lambda[t] - Lambda[t-1]))/(risk.mean[t])
        A = A + Shat[(t-1)]*(ww*(dat$obs.Y == time.list[t]& dat$Delta == l) - ww*(dat$obs.Y >= time.list[t])*(Lambda_l[t,l] - Lambda_l[t-1,l]))/(risk.mean[t])
        B = B + Shat[(t-1)]*aux[,t-1]*(Lambda_l[t,l]-Lambda_l[t-1,l])
      }
    }
    results[,l] = A-B
  }
  
  return(results)
}