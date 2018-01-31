set.seed(123)
####################################################################################################
path= getwd()
workpath=paste0(path,"/data/Thomas data1.Rdata")
load(workpath)
#load("C:/Users/CZhao/Dropbox/Research/PairCorrelationFunction/simulation/SimulationFinal/simulation with selection/setting1/rho1_sigma0.025/rho=1.0/data/Thomas data1.Rdata")
m=length(process)



# first get all process length
process_len = sapply(process, length)
# combine all data points into one big vector 
events_process = unlist(process)
# generate 50 simulation permuation for all data points 

i = 1
while (i <= m){
  
}


























#########################################################################
library(snowfall)
library(snow)
# read data into r
R=0.18
n.h = 8
h.min = 0.002
TT=30
h.max = 0.016
h.list = seq(h.min,h.max,by= 0.002)
lcEstimator<- function(sim){
  #######################
  load(paste0(workpath,"/data/Thomas data",sim,".Rdata"))
  ptm<-proc.time()
  #TT=30
  R.max = R+2*h.max
  #################################################################################
  N  =length(process)
  ptm<-proc.time()
  event_num_N = sapply(process,length)
  event_process = unlist(process)
  # 
  l_event_process = length(event_process)
  dp=list(l_event_process)
  sp=list(l_event_process)
  loc_all = 1:l_event_process
  for(i in 1:N){
    N_i = event_num_N[i]
    sum_N_i = sum(event_num_N[1:i])
    start= sum_N_i-N_i+1
    loc_N_i =start:sum_N_i
    for(j in loc_N_i){
      start_point = event_process[j]
      end_loc_sp = loc_N_i[loc_N_i!=j]
      events_ij_sp = abs(start_point-event_process[ end_loc_sp])
      retain_logi_sp= events_ij_sp<=R.max
      dist_retain_sp=events_ij_sp[retain_logi_sp]
      sp[[j]]=dist_retain_sp
      # different process distance
      events_ij_dp = abs(start_point-event_process[-(start:sum_N_i)])
      retain_logi_dp= events_ij_dp<=R.max
      dist_retain_dp=events_ij_dp[retain_logi_dp] 
      dp[[j]]=dist_retain_dp
    }
  }
  distVS=sort(unlist(sp))
  distVS=distVS[duplicated(distVS)]
  distVD=sort(unlist(dp))
  distVD=distVD[duplicated(distVD)]
  proc.time()-ptm
  ### calculate the distance vectors for points inside the Same process 
  # generate t sequence which represent the whole day 
  n.grid= 500
  grid=seq(0,R,length.out = n.grid)
  estResult = matrix(rep(0,n.grid*n.h),nrow = n.h) 
  for(n in 1:n.h){
    bh=h.list[n]
    # define a function for solve estimating equations
    start_point1 = sapply(grid-bh,function(x) which.max(distVS >= x))
    end_point1   = sapply(grid+bh,function(x) which.min(distVS <= x)-1)
    start_point2 = sapply(grid-bh,function(x) which.max(distVD >= x))
    end_point2   = sapply(grid+bh,function(x) which.min(distVD <= x)-1)
    #############################################################
    for (i in 1:n.grid){
      RetainPair.1=distVS[start_point1[i]:end_point1[i]]
      # retain pairs between different process
      RetainPair.2=distVD[start_point2[i]:end_point2[i]]
      ##############
      A=length(RetainPair.1)*(N-1)
      B=length(RetainPair.2)
      estResult[n,i]=A/B
    }
  }
  estResult
}
#######################
ptm<-proc.time()
sfInit(parallel=TRUE, cpus=cpu_num,type = "SOCK")
# import variables into workers from master
sfExport("lcEstimator","workpath","h.max","R",
         "TT","n.h","h.list","h.min")
# import libraries into workers from master
sfLibrary(snowfall)
#sfLibrary(flexclust)
sfLibrary(snow)
Simu<-sfLapply(1:100,lcEstimator) 
# stop multi core 
sfStop() 
proc.time()-ptm
######  MISE ######
# A,B are 500*100
n.grid= 500
grid=seq(0,R,length.out = n.grid)
gridt=R/n.grid
##### define MISE function to calculate mise ########
MISE=function(t,f.true,fn){
  se  = (fn-f.true)^2
  loc = which.max(grid[t > grid])+1
  se.loc=se[1:loc,] #
  mise=mean(colSums(se.loc*gridt))
  mise
}
##########    True value for pcf ########################
pcf=function(t){
  1+(exp(-(t^2)/((2*sigma)^2))/(rho*sqrt(4*pi)*(sigma)))
}
pcfGrid=pcf(grid)
f.true=matrix(rep(pcfGrid,100),byrow = F,nrow = n.grid)
#####################################################
Simu.array=array(unlist(Simu),c(n.h,n.grid,100)) # dimension 
save(Simu.array,file=paste0(workpath,"/local constant/local constant estimator result.RData"))
#
MISE.h=rep(0,n.h)
for (i in 1: n.h){
  Simu.h =Simu.array[i,,] 
  t=R
  MISE.h[i]=MISE(t,f.true,Simu.h)
}
estResult = Simu.array[which.min(MISE.h),,]
h.optim = h.list[which.min(MISE.h)]
h.optim
#out = list(h.optim,MISE.h,estResult)
save(estResult,file=paste0(workpath,"/local constant/local constant estimator optimal result with bandwidth h=",h.optim,".RData"))
############################################################
for (sim in 1:100){
  file=paste0(workpath,"/local constant/local constant estimator optimal result with bandwidth h=",h.optim,".RData")
  load(file)
  n.grid= 500
  grid=seq(0,R,length.out = n.grid)
  at=estResult[,sim]
  if(sim == 1){
    plot(grid,at,type="l",ylim=c(-1,16),main="Local constant estimator of pcf",xlab = "t",ylab = "pcf")
    lines(grid,rep(1,500))}
  else { lines(grid,at,type="l",ylim=c(-1,16))}
}
lines(grid,pcfGrid,col="red",type="l") 