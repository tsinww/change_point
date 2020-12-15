#### DP

start.all = Sys.time()
print('starttime')
print(start.all)
#lambda.dp.set = c(0.002,0.004,0.006,0.008,0.01)
#gamma.dp.set = c(0.02,0.04,0.06,0.08,0.1)
lambda.dp.set = c(seq(0.001,0.01,0.002),seq(0.011,0.1,0.02),seq(0.11,0.45,0.1))#seq(0.001,0.25,0.001)#c(0.001,0.002)#seq(0.01,0.6,0.1)#c(0.12,0.4)#c(0.1,0.2)#c(0.02,0.04,0.06)#c(0.01,0.02,0.03)#c(0.4,0.5,0.6,0.7,0.8)#c(0.01,0.015,0.02,0.03,0.04)#seq(0.0001,0.001,0.0002)#c(0.0001,0.0002,0.0003)#c(0.01,0.02,0.03,0.05)#0.0001#0.0002#c(0.0006,0.0007,0.0008,0.0009)#c(0.0012,0.0014,0.0016,0.0018)#c(0.016,0.02)
gamma.dp.set = seq(0.05,0.15,0.02)#seq(1,15,2)#c(3,5)#c(0.0025,0.005,0.01)#c(0.04,0.08,0.12,0.16,0.2)#c(0.08)#seq(0.02,0.2,0.04)#c(0.04,0.05,0.06,0.07)#seq(0.01,0.1,0.01)#c(0.03,0.04,0.05,0.06,0.07,0.08)
grid = 1
RR = 1
K = matrix(0,nrow = length(gamma.dp.set),ncol = length(lambda.dp.set))
resmat = matrix(0,nrow = length(gamma.dp.set),ncol = length(lambda.dp.set))
#cp = matrix(0,nrow = length(gamma.dp.set),ncol = length(lambda.dp.set))
init.change.point = matrix(0,20,RR)
len.est = matrix(0,RR,1)
error = matrix(0,RR,1)
nrowmat = length(gamma.dp.set)
ncolmat = length(lambda.dp.set)
X = read.csv("substract_na_dataX_18_365_by_avg.csv")
y = read.csv("substract_na_datay_18_365_by_avg.csv")
X = X[-nrow(X),-1]/sd(as.matrix(y[,2]))
colnames(X) = c()
y = as.matrix(y[,2])/sd(as.matrix(y[,2]))
for (time in 1:RR){
  #data = data.generate.k.fluctuate2(d0,change.point,p,n,sigma,kappa)#data.generate.k.nequal
  #X = data$X
  #y = data$y
  #true = data$true.change.point
  #betafullmat = data$betafullmat
  result = cv.grid.search.dp(lambda.dp.set,gamma.dp.set,X,y,grid)
  K = matrix(as.numeric(result$K),nrowmat,ncolmat)
  resmat = as.numeric(result$res)
  cp = result$changepoint
  #cp = cbind(cp,result$cp)
  beta = result$betaerror
  haus.index = which(matrix(resmat,nrowmat,ncolmat) == min(matrix(resmat,nrowmat,ncolmat)), arr.ind=TRUE)
  #K.avg = abs(K - ncp)
  K.avg = K
  K.from.pair = c()
  #K.from.pair = K.avg[haus.index[i,1],haus.index[i,2]]
  for(i in 1:nrow(haus.index)){
    K.from.pair = c(K.from.pair, K.avg[haus.index[i,1],haus.index[i,2]])
  }
  I = which.min(K.from.pair)
  #rowlabel = haus.index[I,1]; veclabel = haus.index[I,2]
  haus.pair = c(gamma.dp.set[haus.index[I,1]],lambda.dp.set[haus.index[I,2]])#gamma,lambda
  #haus.pair
  #print(cp)
  selected.change.point = as.numeric(unlist(cp[haus.index[I,1],haus.index[I,2]]))#zenmexuanabuzhidao!!
  init.change.point[1:length(selected.change.point), time] = selected.change.point
  #write.csv(init.change.point, file = "kappa_4_dp_d_0_10.csv")
  #error[time,1] = Hausdorff(selected.change.point,change.point)
  #write.csv(error, file = "kappa_4_dp_error_d_0_10.csv")
  len.est[time,1] = K[haus.index[I,1],haus.index[I,2]]
  #write.csv(len.est, file = "kappa_4_dp_ncp_d_0_10.csv")
}
init.change.point
error
len.est
end.all = Sys.time()
print('duration')
print(difftime(end.all, start.all, units = "secs")[[1]])


########## DPLR

start.all = Sys.time()
print('starttime')
print(start.all)
#lambda.dp.set = c(0.002,0.004,0.006,0.008,0.01)
#gamma.dp.set = c(0.02,0.04,0.06,0.08,0.1)
lambda.dp.set = c(seq(0.001,0.01,0.002),seq(0.011,0.1,0.02),seq(0.11,0.45,0.1))#seq(0.001,0.25,0.001)#c(0.001,0.002)#seq(0.01,0.6,0.1)#c(0.12,0.4)#c(0.1,0.2)#c(0.02,0.04,0.06)#c(0.01,0.02,0.03)#c(0.4,0.5,0.6,0.7,0.8)#c(0.01,0.015,0.02,0.03,0.04)#seq(0.0001,0.001,0.0002)#c(0.0001,0.0002,0.0003)#c(0.01,0.02,0.03,0.05)#0.0001#0.0002#c(0.0006,0.0007,0.0008,0.0009)#c(0.0012,0.0014,0.0016,0.0018)#c(0.016,0.02)
gamma.dp.set = seq(0.05,0.15,0.02)#seq(0.05,0.15,0.02)#seq(10,60,5)#seq(0.1,2,0.2)#c(seq(0.01,1,0.1),seq(1.1,2,0.5))#c(0.04,0.08,0.12,0.16,0.2)#c(0.08)#seq(0.02,0.2,0.04)#c(0.04,0.05,0.06,0.07)#seq(0.01,0.1,0.01)#c(0.03,0.04,0.05,0.06,0.07,0.08)
delta = 1
grid = 1
RR = 100
w = 0.9
zeta = seq(0.001,0.01,0.002)#c(0.006)#seq(0.09,0.15,0.02)
K = matrix(0,nrow = length(gamma.dp.set),ncol = length(lambda.dp.set))
resmat = matrix(0,nrow = length(gamma.dp.set),ncol = length(lambda.dp.set))
#cp = matrix(0,nrow = length(gamma.dp.set),ncol = length(lambda.dp.set))
init.change.point = matrix(0,20,RR)
len.est = matrix(0,RR,1)
error = matrix(0,RR,1)
nrowmat = length(gamma.dp.set)
ncolmat = length(lambda.dp.set)
X = read.csv("substract_na_dataX_18_365_by_avg.csv")
y = read.csv("substract_na_datay_18_365_by_avg.csv")
X = X[-nrow(X),-1]/sd(as.matrix(y[,2]))
colnames(X) = c()
y = as.matrix(y[,2])/sd(as.matrix(y[,2]))
for (time in 1:RR){
  result = cv.grid.search.dp.lr(lambda.dp.set,gamma.dp.set,zeta,X,y,w)
  end.all = Sys.time()
  print('duration')
  print(difftime(end.all, start.all, units = "secs")[[1]])
  K = matrix(as.numeric(result[2,]$K),nrowmat,ncolmat)
  resmat = matrix(as.numeric(result[3,]$res),nrowmat,ncolmat)
  cp = result[1,]$changepoint
  #cp = cbind(cp,result$cp)
  #beta = result$betaerror
  haus.index = which(matrix(resmat,nrowmat,ncolmat) == min(matrix(resmat,nrowmat,ncolmat)), arr.ind=TRUE)
  K.avg = K
  K.from.pair = c()
  #K.from.pair = K.avg[haus.index[i,1],haus.index[i,2]]
  for(i in 1:nrow(haus.index)){
    K.from.pair = c(K.from.pair, K.avg[haus.index[i,1],haus.index[i,2]])
  }
  I = which.min(K.from.pair)
  #rowlabel = haus.index[I,1]; veclabel = haus.index[I,2]
  haus.pair = c(gamma.dp.set[haus.index[I,1]],lambda.dp.set[haus.index[I,2]])#gamma,lambda
  #haus.pair
  #print(cp)
  selected.change.point = as.numeric(unlist(cp[haus.index[I,1],haus.index[I,2]]))#zenmexuanabuzhidao!!
  #init.change.point[1:length(selected.change.point), time] = selected.change.point
  #write.csv(init.change.point, file = "1_dplr_d_0_10.csv")
  #error[time,1] = Hausdorff(selected.change.point,change.point)
  #write.csv(error, file = "1_dplr_error_d_0_10.csv")
  len.est[time,1] = K[haus.index[I,1],haus.index[I,2]]

}

L = 20
cpmat = matrix(0,L*nrow(cp),ncol(cp))
for(row in 1:nrow(cp)){
  for (col in 1:ncol(cp)){
    changepoint = unlist(cp[row,col]) 
    cpmat[L*(row-1)+c(1:length(changepoint)),col] =  changepoint
  }
}

selected.change.point
len.est
end.all = Sys.time()
print('duration')
print(difftime(end.all, start.all, units = "secs")[[1]])

