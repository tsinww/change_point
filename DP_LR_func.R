library(gglasso)
library(glmnet)
data.generate2 = function(sparsity,true_change_point,p,n,sigma){
  ncp = 2
  seed = 10
  X = matrix(rnorm(p*n,0,sigma^2),p,n)
  y = matrix(0,n,1)
  d0 = sparsity
  nonzero.element.loc = c(1:d0)
  cp = c(0,true_change_point,n)
  beta1 = matrix(0,p,1)
  beta1[nonzero.element.loc,] = 1
  y[(1+cp[1]):cp[2],] = rnorm(cp[2] - cp[1],t(X[,(1+cp[1]):cp[2]])%*%beta1,sigma)
  beta2 = matrix(0,p,1)
  beta2[nonzero.element.loc,] = 2
  y[(1+cp[2]):cp[3],] = rnorm(cp[3] - cp[2],t(X[,(1+cp[2]):cp[3]])%*%beta2,sigma)
  beta3 = matrix(0,p,1)
  beta3[nonzero.element.loc,] = 3
  y[(1+cp[3]):cp[4],] = rnorm(cp[4] - cp[3],t(X[,(1+cp[3]):cp[4]])%*%beta3,sigma)
  List = list(true.change.point = true_change_point, X = X, y = y)
  return(List)
}

data.generate.k = function(sparsity,ncp,p,n,sigma){
  seed = 10
  cp = seq(0,n,by = n/(ncp+1))
  true_change_point = cp[2:(length(cp)-1)]
  X = matrix(rnorm(p*n,0,1),p,n)
  y = matrix(0,n,1)
  d0 = sparsity
  nonzero.element.loc = c(1:d0)#sort(sample.int(p, d0))
  #cp = c(0,true_change_point,n)
  beta = matrix(0,p,ncp+1)
  betafullmat = matrix(0,p,n)
  for (i in 1:(ncp+1)){
    beta[nonzero.element.loc,i] = i
    y[(1+cp[i]):cp[i+1],] = rnorm(cp[i+1] - cp[i],t(X[,(1+cp[i]):cp[i+1]])%*%beta[,i],sigma)
    for (j in (1+cp[i]):cp[i+1]){
      betafullmat[,j] = beta[,i] 
    }
  }
  
  List = list(true.change.point = true_change_point, X = X, y = y, betafullmat = betafullmat)
  return(List)
}

data.generate.k.nequal = function(d0,change.point,p,n,sigma){
  seed = 10
  #cp = seq(0,n,by = n/(ncp+1))
  #true_change_point = cp[2:(length(cp)-1)]
  true_change_point = change.point
  ncp = length(change.point)
  X = matrix(rnorm(p*n,0,1),p,n)
  y = matrix(0,n,1)
  nonzero.element.loc = c(1:d0)#sort(sample.int(p, d0))
  cp = c(0,true_change_point,n)
  beta = matrix(0,p,ncp+1)
  betafullmat = matrix(0,p,n)
  for (i in 1:(ncp+1)){
    beta[nonzero.element.loc,i] = i*2#*5#/sqrt(d0)#*5
    y[(1+cp[i]):cp[i+1],] = rnorm(cp[i+1] - cp[i],t(X[,(1+cp[i]):cp[i+1]])%*%beta[,i],sigma)
    for (j in (1+cp[i]):cp[i+1]){
      betafullmat[,j] = beta[,i] 
    }
  }
  
  List = list(true.change.point = true_change_point, X = X, y = y, betafullmat = betafullmat)
  return(List)
}

data.generate.k.fluctuate = function(d0,change.point,p,n,sigma){
  seed = 10
  #cp = seq(0,n,by = n/(ncp+1))
  #true_change_point = cp[2:(length(cp)-1)]
  true_change_point = change.point
  ncp = length(change.point)
  X = matrix(rnorm(p*n,0,1),p,n)
  y = matrix(0,n,1)
  nonzero.element.loc = c(1:d0)#sort(sample.int(p, d0))
  cp = c(0,true_change_point,n)
  beta = matrix(0,p,ncp+1)
  betafullmat = matrix(0,p,n)
  for (i in 1:(ncp+1)){
    if (i%%2 == 1){
      beta[nonzero.element.loc,i] = 1
    }
    else{
      beta[nonzero.element.loc,i] = 2
    }
    #beta[nonzero.element.loc,i] = i*10
    y[(1+cp[i]):cp[i+1],] = rnorm(cp[i+1] - cp[i],t(X[,(1+cp[i]):cp[i+1]])%*%beta[,i],sigma)
    for (j in (1+cp[i]):cp[i+1]){
      betafullmat[,j] = beta[,i] 
    }
  }
  
  List = list(true.change.point = true_change_point, X = X, y = y, betafullmat = betafullmat)
  return(List)
}

data.generate.k.fluctuate2 = function(d0,change.point,p,n,sigma,kappa){
  seed = 10
  #cp = seq(0,n,by = n/(ncp+1))
  #true_change_point = cp[2:(length(cp)-1)]
  true_change_point = change.point
  ncp = length(change.point)
  X = matrix(rnorm(p*n,0,1),p,n)
  y = matrix(0,n,1)
  nonzero.element.loc = c(1:d0)#sort(sample.int(p, d0))
  cp = c(0,true_change_point,n)
  beta = matrix(0,p,ncp+1)
  betafullmat = matrix(0,p,n)
  for (i in 1:(ncp+1)){
    if (i%%2 == 1){
      beta[nonzero.element.loc,i] = kappa/(2*sqrt(d0))
    }
    else{
      beta[nonzero.element.loc,i] = -kappa/(2*sqrt(d0))
    }
    #beta[nonzero.element.loc,i] = i*10
    y[(1+cp[i]):cp[i+1],] = rnorm(cp[i+1] - cp[i],t(X[,(1+cp[i]):cp[i+1]])%*%beta[,i],sigma)
    for (j in (1+cp[i]):cp[i+1]){
      betafullmat[,j] = beta[,i] 
    }
  }
  
  List = list(true.change.point = true_change_point, X = X, y = y, betafullmat = betafullmat)
  return(List)
}

DPR = function(gamma,lambda,dataX, datay,gridsize,sigma,ncp,sparsity){
  ## N must be the multiple of gridsize
  N = ncol(dataX)
  p = nrow(dataX)
  Bestvalue = c()
  partition = c()
  Bestvalue = c(Bestvalue,-gamma*(ncp+1)* (sigma^2) * (sparsity^2) *log(max(N,p)))
  #gridsize = 1## gridsize cannot control error, will cause larger error
  for (r in 1:N){
    #b = c()
    b = Inf
    if (r %% gridsize == 0){
      #b = c()
      b = sapply(1:r,function(l)inner.func(Bestvalue,gamma,lambda,dataX,datay,sigma,sparsity,l,r,ncp,gridsize))
      #print(b)
      #for (l in 1:r){
      #btemp = Inf
      #if (l %% gridsize == 0){
      #btemp = Bestvalue[l] + gamma*(ncp+1)* (sigma^2) * (sparsity^2) *log(max(N,p)) + distanceR(l,r,dataX,datay,lambda,sigma,sparsity)
      #}
      ##btemp = Bestvalue[l] + gamma + distanceR(l, r,dataX,datay,lambda)
      #b = c(b, btemp)
      #}
    }
    Bestvalue = c(Bestvalue,min(b))
    partition = c(partition, which.min(b) - 1)
  }
  R = N
  L = c()
  while (R > 0) {
    L = c(partition[R], L)
    R = partition[R]
  }
  L = L[-1]
  return(L)
}

DPR1 = function(gamma,lambda,dataX, datay,gridsize,sigma,ncp,sparsity){
  ## N must be the multiple of gridsize
  N = ncol(dataX)
  p = nrow(dataX)
  Bestvalue = c()
  partition = c()
  Bestvalue = c(Bestvalue,-gamma*(ncp+1)* (sigma^2) * (sparsity^2) *log(max(N,p)))
  #gridsize = 1## gridsize cannot control error, will cause larger error
  for (r in 1:N){
    Best.value[r+1]=Inf
    #b = Inf
    if (r %% gridsize == 0){
      #b = c()
      b = sapply(1:r,function(l)inner.func(Bestvalue,gamma,lambda,dataX,datay,sigma,sparsity,l,r,ncp,gridsize))
      #print(b)
      #for (l in 1:r){
      #btemp = Inf
      #if (l %% gridsize == 0){
      #btemp = Bestvalue[l] + gamma*(ncp+1)* (sigma^2) * (sparsity^2) *log(max(N,p)) + distanceR(l,r,dataX,datay,lambda,sigma,sparsity)
      #}
      ##btemp = Bestvalue[l] + gamma + distanceR(l, r,dataX,datay,lambda)
      #b = c(b, btemp)
      #}
    }
    #Bestvalue = c(Bestvalue,min(b))
    #partition = c(partition, which.min(b) - 1)
    if (min(b) <Best.value[r+1]){
      Best.value[r+1]=min(b)
      partition[r+1]=which.min(b)-1
    } 
  }
  R = N
  L = c()
  while (R > 0) {
    L = c(partition[R], L)
    R = partition[R]
  }
  L = L[-1]
  return(L)
}                 

inner.func = function(Bestvalue,gamma,lambda,dataX,datay,sigma,sparsity,l,r,ncp,gridsize){
  btemp = Inf
  N = ncol(dataX)
  p = nrow(dataX)
  if (l %% gridsize == 0){
    btemp = Bestvalue[l] + gamma*(ncp+1)* (sigma^2) * (sparsity^2) *log(max(N,p)) + distanceR(l,r,dataX,datay,lambda,sigma,sparsity)
  }
  #btemp = Bestvalue[l] + gamma + distanceR(l, r,dataX,datay,lambda)
  #b = c(b, btemp)
  return(btemp)
}
# residual
distanceR = function(s, e, dataX,datay,lambda,sigma,sparsity){
  n = ncol(dataX)
  if (abs(s-e) >0){
    fit = glmnet(x=t(dataX[, s:e]), y=datay[s:e,], family=c("gaussian"),
                 alpha = 1, lambda=lambda*sigma*sqrt(max(log(n),e-s))*sqrt(sparsity*log(n))/(e-s),intercept=F)#lambda=lambda/sqrt(max(log(p),e-s))#lambda*sigma*sqrt(sparsity*log(max(e-s+1, p)))
    coef_est = t(t(as.vector(fit$beta)))
    yhat = t(dataX[,s:e])%*%coef_est
    d = norm(datay[s:e] - yhat, type = "2")
  }
  else{
    d = 0
  }
  return(d^2)
}
## divide X into [s,eta] and [eta+1,e]
convert.design.matrix.one.change=function(X,eta,s_ceiling){
  ee=ncol(X)
  xx1=t(X)
  t = eta - s_ceiling +1
  xx1[ (t+1):ee,]=0##
  xx2=t(X)
  xx2[1:t,]=0
  #xx = cbind(xx1,xx2)
  xx=cbind(xx1/sqrt(t-1),xx2/sqrt(ee-t+1))
  return(xx)
}

LocalRefine = function(init_cp,X,y,w,zeta){
  partition = c(0)
  n = ncol(X)
  #w = 18/19# 5/7
  #lambda.LR = 0.03#*sqrt(log(nrow(X)))
  init_cp_ = c(0,init_cp,n)
  K = length(init_cp_) - 2
  for (k in 1:K){
    #s = init_cp_[k]
    #e = init_cp_[k+2]
    #s = w*init_cp_[k] + (1-w)*init_cp_[k+1]#w = 0.9
    s = w*partition[k] + (1-w)*init_cp_[k+1]
    e = (1-w)*init_cp_[k+1] + w*init_cp_[k+2]#w = 0.9  
    b = c()
    lower = ceiling(s)+1
    upper = floor(e)-1
    #print(lower)
    #print(upper)
    b = sapply(lower:upper,function(eta)inner.func.lr(s,e,eta,X,y,zeta))
    #for (eta in (ceiling(s)+1):(floor(e)-1)){
    #group = rep(1:nrow(X),2)
    #convertX = convert.design.matrix.one.change(X[,(ceiling(s)):(floor(e))],eta,ceiling(s))
    #y_ = y[(ceiling(s)):(floor(e)),]
    #lambda.LR = cv.gglasso(x=convertX,y=y_,group=group, loss="ls",pred.loss="L1")$lambda.min
    ##lambda.LR = 0.001*sqrt(log(300))
    #auxfit = gglasso(x=convertX,y=y_,group=group, loss="ls",
    #lambda=lambda.LR,intercept = FALSE,eps =
    #0.001)#*sqrt(log(max(p,ncol(X))))/(floor(e)-ceiling(s)+1)#lambda=lambda.LR/ncol(X)#lambda.LR*sqrt(log(max(p,floor(e)-ceiling(s)+1)))
    #coef = as.vector(auxfit$beta)
    #btemp = norm(y_ - convertX %*% coef, type = "2")
    #b = c(b, btemp)
    #}
    partition = c(partition, ceiling(s) + which.min(b) )
    #init_cp_[k+1] = ceiling(s) + which.min(b) 
  }
  return(partition[-1])
}

inner.func.lr = function(s,e,eta,X,y,zeta){
  p = nrow(X)
  group = rep(1:nrow(X),2)
  convertX = convert.design.matrix.one.change(X[,(ceiling(s)):(floor(e))],eta,ceiling(s))
  y_ = y[(ceiling(s)):(floor(e)),]
  #lambda.LR = cv.gglasso(x=convertX,y=y_,group=group, loss="ls",pred.loss="L1")$lambda.min
  lambda.LR = zeta*sqrt(log(max(ncol(X),nrow(X))))
  auxfit = gglasso(x=convertX,y=y_,group=group, loss="ls",
                   lambda=lambda.LR/(floor(e)-ceiling(s)+1),intercept = FALSE,eps =
                     0.001)#*sqrt(log(max(p,ncol(X))))/(floor(e)-ceiling(s)+1)#lambda=lambda.LR/ncol(X)#lambda.LR*sqrt(log(max(p,floor(e)-ceiling(s)+1)))
  coef = as.vector(auxfit$beta)
  coef1 = coef[1:p]
  coef2 = coef[(p+1):(2*p)]
  #print(p)
  #print(lambda.LR*sum(sqrt(coef1^2 + coef2^2)))
  btemp = norm(y_ - convertX %*% coef, type = "2")^2 + lambda.LR*sum(sqrt(coef1^2 + coef2^2))
  #print(btemp)
  return(btemp)
}


measurement = function(est,true){
  dist = matrix(0,nrow = length(est),ncol = length(true))
  for (i in 1:nrow(dist)){
    for (j in 1:ncol(dist)){
      dist[i,j] = abs(est[i] - true[j])
    }
  }
  
  error = mean(apply(dist, 2, function(x) min(x) ))
  return(error)
}
Hausdorff = function(vec1,vec2){
  dist = matrix(0,nrow = length(vec1),ncol = length(vec2))
  for (i in 1:nrow(dist)){
    for (j in 1:ncol(dist)){
      dist[i,j] = abs(vec1[i] - vec2[j])
    }
  }
  
  dH = max(max(apply(dist, 2, function(x) min(x) )), max(apply(dist, 1, function(x) min(x) )))
  return(dH)
}



#training
distance = function(s, e, X,y, dataX,datay,lambda,sigma,d0){
  n = ncol(dataX)
  p = nrow(dataX)
  #even_indexes = seq(2,ncol(dataX),2)
  #odd_indexes = seq(1,(ncol(dataX)-1),2)
  lower = s-1#odd_indexes[s]-1
  upper = e#odd_indexes[e]
  #print(s)
  #print(e)
  #print()
  #print(lower)
  #print(upper)
  #print(ncol(X))
  print(lambda*sigma*sqrt(max(log(max(n,p)),upper-lower))*sqrt(d0*log(max(n,p)))/(upper-lower))
  fit = glmnet(x=t(dataX[, lower:upper]), y=datay[lower:upper,], family=c("gaussian"),alpha = 1, 
               lambda=lambda*sigma*sqrt(max(log(max(n,p)),upper-lower))*sqrt(d0*log(max(n,p)))/(upper-lower),intercept=F)
  #fit = glmnet(x=t(dataX[, s:e]), y=datay[s:e,], family=c("gaussian"),
  #             alpha = 1, lambda=lambda,intercept=F)#lambda*sigma*sqrt(max(log(max(n,p)),e-s))*sqrt(d0*log(max(n,p)))
  
  coef_est = t(t(as.vector(fit$beta)))
  yhat = t(dataX[,s:e])%*%coef_est
  d = norm(datay[s:e,] - yhat, type = "2")
  result = list("mse" = d^2, "beta" = coef_est)
  return(result)
}

residual = function(lower, upper,dataX,datay,beta_est){
  #beta_est = betamat[,index]
  res = norm(datay[lower:upper] - t(dataX[,lower:upper])%*%beta_est, type = "2")^2
  return(res)
} 


cv.grid.dp = function(lambda,gamma,X,y,d0,grid,sigma,ncp,betafullmat,GT){
  n = ncol(X)
  even_indexes = seq(2,n,2)
  odd_indexes = seq(1,(n-1),2)
  train.X = X[,odd_indexes]
  train.y = y[odd_indexes,]
  train.y = as.matrix(train.y)
  validation.X = X[,even_indexes]
  validation.y = y[even_indexes,]
  betamatfull = betafullmat[,odd_indexes]
  colnames(train.X) = c()
  row.names(train.y) = c()
  init_cp_train = DPR(gamma,lambda,train.X, train.y,grid,sigma,ncp,d0)#c(75,135)
  init_cp = c()
  for (i in init_cp_train){
    init_cp = c(init_cp,odd_indexes[i])
  }
  print("error")
  print(Hausdorff(init_cp,GT))
  print(init_cp)
  print(c(lambda,gamma))
  #how to estimate beta????
  len = length(init_cp)
  init_cp_long = c(init_cp_train,n/2)
  interval = matrix(0,nrow = len + 1,ncol = 2)
  interval[1,] = c(1,init_cp_long[1])
  if (len > 0){
    for (j in 2:(1+len)){
      interval[j,] = c(init_cp_long[j-1]+1,init_cp_long[j])
    }
  }
  p = nrow(train.X)
  trainmat = sapply(1:(len+1),function(index) distance(interval[index,1], interval[index,2], X, y, train.X,train.y,lambda,sigma,d0))
  #trainmat = sapply(1:(len+1),function(index) distance(interval[index,1], interval[index,2], train.X,train.y,lambda,sigma,d0))
  betamat = matrix(0,nrow = p,ncol = len+1)
  training_loss = matrix(0,nrow = 1,ncol = len+1)                
  for (col in 1:(len+1)){
    betamat[,col] = as.numeric(trainmat[2,col]$beta)
    training_loss[,col] = as.numeric(trainmat[1,col]$mse)
  }
  betafullmatest = matrix(0,nrow = p,ncol = n/2)
  for(index in 1:(len+1)){
    betafullmatest[,interval[index,1]:interval[index,2]] = betamat[,index]
  } 
  deltabeta = betafullmatest - betamatfull
  betaerror = 0
  betaerrorvec = c()                  
  for (column in 1:ncol(deltabeta)){
    betaerror = betaerror + norm(deltabeta[,column],type = "2")^2
    betaerrorvec = c(betaerrorvec,norm(deltabeta[,column],type = "2")^2)
  }
  #print(unique(betaerrorvec))                  
  #print(betamat)       
  validationmat = sapply(1:(len+1),function(index) residual(interval[index,1], interval[index,2], validation.X,validation.y,betamat[,index]))
  #training_loss = sapply(1:(len+1),function(index) distance(interval[index,1], interval[index,2], validation.X,validation.y,betamat[,index]))
  result = list("cp" = init_cp, "K" = len, "res" = sum(validationmat),"trloss" = sum(training_loss),"betaerror" = betaerror)                       
  return(result)
}

cv.grid.search.dp = function(lambda.set,gamma.set,X,y,d0,grid,sigma,ncp,betafullmat,GT){
  output = sapply(1:length(lambda.set), function(i) sapply(1:length(gamma.set), 
                                                           function(j) cv.grid.dp(lambda.set[i],gamma.set[j],X,y,d0,grid,sigma,ncp,betafullmat,GT)))                                                 
  #cp[1:length(init_cp),j+length(gamma.set)*(i-1)] = init_cp
  print(output)
  cp = output[seq(1,5*length(gamma.set),5),]
  len = output[seq(2,5*length(gamma.set),5),]
  ply = output[seq(3,5*length(gamma.set),5),]
  trloss = output[seq(4,5*length(gamma.set),5),] 
  betaerror = output[seq(5,5*length(gamma.set),5),]                                                         
  result = list("changepoint" = cp, "K" = len, "res" = ply, "trloss" = trloss, "betaerror" = betaerror)
  return(result)
}                           


# CV DP
distance = function(s, e, X,y, dataX,datay,zeta,delta){
  #dataX is training X not X
  n = ncol(dataX)
  p = nrow(dataX)
  #even_indexes = seq(2,ncol(dataX),2)
  #odd_indexes = seq(1,(ncol(dataX)-1),2)
  lower = s-1#odd_indexes[s]-1
  upper = e#odd_indexes[e]
  #print(s)
  #print(e)
  #print()
  #print(lower)
  #print(upper)
  group = c(1:nrow(dataX))#rep(1:nrow(dataX),2)
  convertX = convert.matrix.1.group(dataX[,(ceiling(s)):(floor(e))])
  y_ = datay[(ceiling(s)):(floor(e)),]
  #lambda.LR = cv.gglasso(x=convertX,y=y_,group=group, loss="ls",pred.loss="L1")$lambda.min
  lambda.LR = zeta*sqrt(log(max(n,p)))
  #auxfit = gglasso(x=convertX,y=y_,group=group, loss="ls",
  #                   lambda=lambda.LR/(floor(e)-ceiling(s)+1),intercept = FALSE,eps =
  #                     0.001)#*sqrt(log(max(p,ncol(X))))/(floor(e)-ceiling(s)+1)#lambda=lambda.LR/ncol(X)#lambda.LR*sqrt(log(max(p,floor(e)-ceiling(s)+1)))
  fit = glmnet(x=t(dataX[,lower:upper]), y=datay[lower:upper,], family=c("gaussian"),
               alpha = 1,lambda = lambda.LR,intercept=F)#lambda=lambda*sigma*sqrt(max(log(max(n,p)),e-s))*sqrt(d0*log(max(n,p)))
  coef_est = t(t(as.vector(fit$beta)))
  #print(coef_est)
  #coef_est = as.vector(auxfit$beta)
  #coef1 = coef[1:p]
  #coef2 = coef[(p+1):(2*p)]
  #d = norm(y_ - convertX %*% coef_est, type = "2")^2 + lambda.LR*sum(sqrt(coef_est^2))
  #print(ncol(X))
  #print(lambda*sigma*sqrt(max(log(max(n,p)),upper-lower))*sqrt(d0*log(max(n,p)))/(upper-lower))
  #fit = glmnet(x=t(dataX[, lower:upper]), y=datay[lower:upper,], family=c("gaussian"),alpha = 1, 
  #lambda=lambda*sigma*sqrt(max(log(max(n,p)),upper-lower))*sqrt(d0*log(max(n,p)))/(upper-lower),intercept=F)
  #fit = glmnet(x=t(dataX[, s:e]), y=datay[s:e,], family=c("gaussian"),
  #             alpha = 1, lambda=lambda,intercept=F)#lambda*sigma*sqrt(max(log(max(n,p)),e-s))*sqrt(d0*log(max(n,p)))
  
  #coef_est = t(t(as.vector(fit$beta)))
  #yhat = t(dataX[,s:e])%*%coef_est
  #d = norm(datay[s:e,] - yhat, type = "2")
  yhat = t(dataX[,s:e])%*%coef_est
  d = norm(datay[s:e,] - yhat, type = "2")
  result = list("mse" = d^2, "beta" = coef_est)
  return(result)
}

convert.matrix.1.group=function(X){
  ee=ncol(X)
  xx1=t(X)
  #t = eta - s_ceiling +1
  #xx1[ (t+1):ee,]=0##
  #xx2=t(X)
  #xx2[1:t,]=0
  xx=cbind(xx1/sqrt(t-1),xx2/sqrt(ee-t+1))
  return(xx)
}

residual = function(lower, upper,dataX,datay,beta_est){
  #beta_est = betamat[,index]
  res = norm(datay[lower:upper] - t(dataX[,lower:upper])%*%beta_est, type = "2")^2
  return(res)
} 

## CV DP.LR

cv.grid.pb.lr = function(lambda,gamma,zeta,X,y,betafullmat,w){
  n = ncol(X)
  even_indexes = seq(2,n,2)
  odd_indexes = seq(1,(n-1),2)
  train.X = X[,odd_indexes]
  
  train.y = y[odd_indexes,]
  train.y = as.matrix(train.y)
  validation.X = X[,even_indexes]
  validation.y = y[even_indexes,]
  colnames(train.X) = c()
  row.names(train.y) = c()
  betamatfull = betafullmat[,odd_indexes]
  init_cp_train.pb = DPR(gamma,lambda,train.X, train.y,grid,sigma,ncp,d0)#c(75,135)
  #init_cp_train = LocalRefine(init_cp_train.pb, train.X, train.y,w,zeta)
  init_cp.pb = c()
  for (i in init_cp_train.pb){
    init_cp.pb = c(init_cp.pb,odd_indexes[i])
  }
  init_cp = LocalRefine(init_cp.pb, X, y,w,zeta)
  print("error")
  print(Hausdorff(init_cp.pb,c(120,220,350,450)))
  print(init_cp.pb)
  print(Hausdorff(init_cp,c(120,220,350,450)))
  print(init_cp)
  print(c(lambda,gamma))
  #how to estimate beta????
  len = length(init_cp)
  init_cp_train = (1+init_cp)/2
  init_cp_long = c(init_cp_train,n/2)
  interval = matrix(0,nrow = len + 1,ncol = 2)
  interval[1,] = c(1,init_cp_long[1])
  if (len > 0){
    for (j in 2:(1+len)){
      interval[j,] = c(init_cp_long[j-1]+1,init_cp_long[j])
    }
  }
  p = nrow(train.X)
  trainmat = sapply(1:(len+1),function(index) distance(interval[index,1], interval[index,2], X,y,train.X,train.y,zeta,delta))
  betamat = matrix(0,nrow = p,ncol = len+1)
  training_loss = matrix(0,nrow = 1,ncol = len+1)                
  for (col in 1:(len+1)){
    betamat[,col] = as.numeric(trainmat[2,col]$beta)
    training_loss[,col] = as.numeric(trainmat[1,col]$mse)
  }
  betafullmatest = matrix(0,nrow = p,ncol = n/2)
  for(index in 1:(len+1)){
    betafullmatest[,interval[index,1]:interval[index,2]] = betamat[,index]
  } 
  deltabeta = betafullmatest - betamatfull
  betaerror = 0
  betaerrorvec = c()                  
  for (column in 1:ncol(deltabeta)){
    betaerror = betaerror + norm(deltabeta[,column],type = "2")^2
    betaerrorvec = c(betaerrorvec,norm(deltabeta[,column],type = "2")^2)
  }
  #print(unique(betaerrorvec))                  
  #print(betamat)       
  validationmat = sapply(1:(len+1),function(index) residual(interval[index,1], interval[index,2], validation.X,validation.y,betamat[,index]))
  #training_loss = sapply(1:(len+1),function(index) distance(interval[index,1], interval[index,2], validation.X,validation.y,betamat[,index]))
  result = list("cp" = init_cp, "K" = len, "res" = sum(validationmat),"trloss" = sum(training_loss),"betaerror" = betaerror)                       
  return(result)
}

cv.grid.search.pb.lr.lg = function(lambda.set,gamma.set,zeta,X,y,betafullmat,w){
  output = sapply(1:length(lambda.set), function(i) sapply(1:length(gamma.set), 
                                                           function(j) cv.grid.pb.lr(lambda.set[i],gamma.set[j],zeta, X,y,betafullmat,w)))                                                 
  #cp[1:length(init_cp),j+length(gamma.set)*(i-1)] = init_cp
  print(output)
  cp = output[seq(1,5*length(gamma.set),5),]
  len = output[seq(2,5*length(gamma.set),5),]
  ply = output[seq(3,5*length(gamma.set),5),]
  trloss = output[seq(4,5*length(gamma.set),5),] 
  betaerror = output[seq(5,5*length(gamma.set),5),]                                                         
  result = list("changepoint" = cp, "K" = len, "res" = ply, "trloss" = trloss, "betaerror" = betaerror)
  return(result)
}                           

cv.grid.search.pb.lr = function(lambda.set,gamma.set,zeta.set,X,y,delta,betafullmat,w){
  output.2 = sapply(1:length(zeta.set), function(q) cv.grid.search.pb.lr.lg(lambda.set,gamma.set,zeta.set[q], X,y,betafullmat,w))
  print("output with zeta")
  print(output.2) 
  cp = output.2[1,]#[seq(1,5*length(gamma.set),5),]
  len = output.2[2,]#[seq(2,5*length(gamma.set),5),]
  ply = output.2[3,]#[seq(3,5*length(gamma.set),5),]
  trloss = output.2[4,]#[seq(4,5*length(gamma.set),5),] 
  betaerror = output.2[5,]#[seq(5,5*length(gamma.set),5),]                    
  #return(list("changepoint" = cp, "K" = len, "res" = ply, "trloss" = trloss, "betaerror" = betaerror))
  return(output.2)                  
}                   
