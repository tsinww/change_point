#main part
# remember to change d0,lambda,gamma and the names of csv files
n = 600
p = 200
start.all = Sys.time()
print('starttime')
print(start.all)
m = 1
d0 = 15##remember to change
RR = 25
w = 0.9
grid = 1
sigma = 1
ncp = 4
kappa = 5
change.point = c(120,220,350,450)
gamma.dp = c(0.001,0.005,0.01)
lambda.dp = c(0.05,0.1,0.5,1)
zeta.dp = c(0.01,0.05,0.1)#0.05
start.dp = Sys.time()
ply.dp = matrix(0,1,RR)
ply.lr = matrix(0,1,RR)
time.dp = matrix(0,RR,1)
time.lr = matrix(0,RR,1)
ply.dp.1 = matrix(0,1,RR)
ply.lr.1 = matrix(0,1,RR)
len.est = matrix(0,RR,1)
init.change.point = matrix(0,20,RR)
lr.change.point = matrix(0,20,RR)
sigma = 1
d0 = 20
ncp = 4
grid = 1
RR = 10
p = 200 
n = 600#300*3
kappa = 6
change.point = c(120,220,350,450)
K = matrix(0,nrow = length(gamma.dp.set),ncol = length(lambda.dp.set))
resmat = matrix(0,nrow = length(gamma.dp.set),ncol = length(lambda.dp.set))
#cp = matrix(0,nrow = length(gamma.dp.set),ncol = length(lambda.dp.set))
init.change.point = matrix(0,20,RR)
len.est = matrix(0,RR,1)
error = matrix(0,RR,1)
nrowmat = length(gamma.dp.set)
ncolmat = length(lambda.dp.set)
for (time in 1:RR){
  data = data.generate.k.fluctuate2(d0,change.point,p,n,sigma,kappa)#data.generate.k.nequal
  X = data$X
  y = data$y
  true = data$true.change.point
  betafullmat = data$betafullmat
  result = cv.grid.search.dp(lambda.dp.set,gamma.dp.set,X,y,d0,grid,sigma,ncp,betafullmat,change.point)
  K = matrix(as.numeric(result$K),nrowmat,ncolmat)
  resmat = as.numeric(result$res)
  cp = result$changepoint
  #cp = cbind(cp,result$cp)
  beta = result$betaerror
  haus.index = which(matrix(resmat,nrowmat,ncolmat) == min(matrix(resmat,nrowmat,ncolmat)), arr.ind=TRUE)
  K.avg = abs(K - ncp)
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
  write.csv(init.change.point, file = "kappa_6_dp_d_0_20.csv")
  error[time,1] = Hausdorff(selected.change.point,change.point)
  write.csv(error, file = "kappa_6_dp_error_d_0_20.csv")
  #Hausdorff(selected.change.point,change.point)
  len.est[time,1] = K[haus.index[I,1],haus.index[I,2]]
  write.csv(len.est, file = "kappa_6_dp_ncp_d_0_20.csv")
}
init.change.point
error
len.est
end.all = Sys.time()
print('duration')
print(difftime(end.all, start.all, units = "secs")[[1]])


## DP.LR

# first try of cross validation
start.all = Sys.time()
print('starttime')
print(start.all)
#lambda.dp.set = c(0.002,0.004,0.006,0.008,0.01)
#gamma.dp.set = c(0.02,0.04,0.06,0.08,0.1)
#gamma.dp = c(0.001,0.005,0.01)#0.008
#lambda.dp = c(0.05,0.1,0.5,1)#0.5#c(0.04,0.08,0.12,0.16,0.2)#c(0.08)#seq(0.02,0.2,0.04)#c(0.04,0.05,0.06,0.07)#seq(0.01,0.1,0.01)#c(0.03,0.04,0.05,0.06,0.07,0.08)
sigma = 1
d0 = 10
ncp = 4
grid = 1
RR = 40
p = 200 
n = 600#300*3
kappa = 5
w = 0.9
zeta = c(0.05)
change.point = c(120,220,350,450)
K = matrix(0,nrow = length(gamma.dp.set),ncol = length(lambda.dp.set))
resmat = matrix(0,nrow = length(gamma.dp.set),ncol = length(lambda.dp.set))
#cp = matrix(0,nrow = length(gamma.dp.set),ncol = length(lambda.dp.set))
init.change.point = matrix(0,20,RR)
len.est = matrix(0,RR,1)
error = matrix(0,RR,1)
nrowmat = length(gamma.dp.set)
ncolmat = length(lambda.dp.set)
for (time in 1:RR){
  data = data.generate.k.fluctuate2(d0,change.point,p,n,sigma,kappa)#data.generate.k.nequal
  X = data$X
  y = data$y
  true = data$true.change.point
  betafullmat = data$betafullmat#lambda.set,gamma.set,zeta.set,X,y,delta,betafullmat,w
  result = cv.grid.search.pb.lr(lambda.dp.set,gamma.dp.set,zeta,X,y,betafullmat,w)
  end.all = Sys.time()
  print('duration')
  print(difftime(end.all, start.all, units = "secs")[[1]])
  K = matrix(as.numeric(result[2,]$K),nrowmat,ncolmat)
  resmat = matrix(as.numeric(result[3,]$res),nrowmat,ncolmat)
  cp = result[1,]$changepoint
  #cp = cbind(cp,result$cp)
  #beta = result$betaerror
  haus.index = which(matrix(resmat,nrowmat,ncolmat) == min(matrix(resmat,nrowmat,ncolmat)), arr.ind=TRUE)
  K.avg = abs(K - ncp)
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
  selected.change.point = as.numeric(unlist(cp[haus.index[I,1],haus.index[I,2]]))
  init.change.point[1:length(selected.change.point), time] = selected.change.point
  write.csv(init.change.point, file = "1_dplr_d_0_10.csv")
  error[time,1] = Hausdorff(selected.change.point,change.point)
  write.csv(error, file = "1_dplr_error_d_0_10.csv")
  #Hausdorff(selected.change.point,change.point)
  len.est[time,1] = K[haus.index[I,1],haus.index[I,2]]
  write.csv(len.est, file = "1_dplr_ncp_d_0_10.csv")
}
end.all = Sys.time()
print('duration')
print(difftime(end.all, start.all, units = "secs")[[1]])