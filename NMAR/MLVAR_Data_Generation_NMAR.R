# This script generates simulated datasets under the low-stability (i.e., low autocorrelation) condition. 
# Missing data are generated following the missing at random (MAR) mechanism.
# The data used in the model fitting with Stan and JAGS are saved in the ".Rdata" files; 
# The data used in the model fitting with Mplus are saved in the ".dat" files.

library(MASS)
library(tmvtnorm)

r1=1
N=100
O=60
high=0

set.seed(r1)

###################### Data Simulation ######################################################
N = N #number of subjects
O = O #number of time points
O = O+100
K = 2 #number of person-level covariates (including intercept)
D = 2 #number of observed variables

#A covariates matrix(including intercept)
A = matrix(rep(NA, K*N), nrow = N, ncol=K)
for (i in 1:N) {
  A[i,] <- c(1, rnorm(n = 1,mean = 0,sd = 1)) #covariates set to standard normal
}
AforInter = A[,1:K]
AforAR = A[,1:K]
AforlogSD_VAR = A[,1:K]
Aforcorr_VAR = A[,1:K]

# regression coefficients for person-specific parameters
# y1(t) = Mu1 + a1*y1(t-1) + b1*y2(t-1) + e1
# y2(t) = Mu2 + b2*y1(t-1) + a2*y2(t-1) + e2 

CoeffMu1 = c(.5,-.3)
CoeffMu2 = c(.5,.3)
if(high==0){
  Coeffa1 = c(.3,-.2)  
  Coeffa2 = c(.3,.2)
}

if(high==1){
  Coeffa1 = c(.6,-.2)  
  Coeffa2 = c(.5,.2)  
}

Coeffb1 = c(-.1,-.2) 
Coeffb2 = c(-.1,.2)
CoefflogSD_VAR_1 = c(.5,-.3) 
CoefflogSD_VAR_2 = c(.5,.3)
Coeffcorr_VAR = c(.5,-.3)

CoeffInter = matrix(rep(0,K*D),ncol=D)  
CoeffAR = array(rep(0,K*D*D),c(K,D,D))  
CoefflogSD_VAR = matrix(rep(0,K*D), ncol=D)
CoeffInter[,1] = CoeffMu1
CoeffInter[,2] = CoeffMu2
CoeffAR[,1,1] = Coeffa1
CoeffAR[,2,2] = Coeffa2
CoeffAR[,1,2] = Coeffb1
CoeffAR[,2,1] = Coeffb2
CoefflogSD_VAR[, 1] = CoefflogSD_VAR_1
CoefflogSD_VAR[, 2] = CoefflogSD_VAR_2

# Simulate covariance matrix of person-specific intercept and AR parameters 
# order of parameter: mu1, mu2, a1, a2, b1, b2
sigmaInter = 0.3
sigmaAR = 0.1
random_SD = diag(c(rep(sigmaInter,2),rep(sigmaAR,4)))
random_corr = diag(1,6)
random_corr[1,2] = -0.5 # corr(mu1,mu2)
random_corr[2,1] = random_corr[1,2] 
random_corr[1,3] = 0.5 # corr(mu1,a1)
random_corr[3,1] = random_corr[1,3]
random_corr[2,4] = 0.5 # corr(mu2,a2)
random_corr[4,2] = random_corr[2,4]
bcov = random_SD %*% random_corr %*% random_SD

# Simulate covariance matrix of VAR process
S <- array(rep(0, D*D*N), dim = c(D, D, N))
logS <- array(rep(0, D*N), dim = c(D, N))
logS_Means <- matrix(rep(0,N*D),nrow = N) 
R_Means <- c()
logS_SD <- c(.3, .3)
R_SD <- .3

#first, simulate diagonal matrices of standard deviations
for (i in 1:N) {
  for (d in 1:D) {
    logS_Means[i,d] = sum(CoefflogSD_VAR[,d]*AforlogSD_VAR[i,]);
    logS[d,i] <- rnorm(1, mean = logS_Means[i,d], sd = logS_SD[d]);
    S[d,d,i] <- exp(logS[d,i]);
  }}


#second, simulate correlation matrices
R <- array(rep(1, D*D*N), dim = c(D, D, N))
for(i in 1:N){
  R_Means[i] <- sum(Coeffcorr_VAR*Aforcorr_VAR[i,]);
  for(d in 1:(D-1)){
    for(j in (d+1):D){
      R[d,j,i] <- rnorm(n=1, mean = R_Means[i], sd = R_SD); 
      R[d,j,i] <- (exp(2*R[d,j,i])-1)/ (exp(2*R[d,j,i])+1); #Fisher Z to R
      R[j,d,i] <- R[d,j,i];
    }
  } 
}


#Third, generate covariance matrices
Sigma_VAR <- array(rep(NA, D*D*N), dim = c(D,D,N))
for(i in 1:N){
  Sigma_VAR[,,i] <- S[,,i] %*% R[,,i] %*% S[,,i]
}


# Simulate observed variables
Y = array(rep(0, N*O*D), c(N,O,D))  # Y
mus = array(rep(0, N*O*D), c(N,O,D)) # mean of Y

muInter = matrix(rep(0,N*D),nrow = N)   #mean of person-specific intercepts
intercept = matrix(rep(0,N*D),nrow = N) #person-specific intercepts

muAR = array(rep(0, N*D*D), c(N,D,D)) #mean of person-specific AR coefficients
AR = array(rep(0, N*D*D), c(N,D,D))  #person-specific AR coefficients

bmu = matrix(rep(0,N*6),nrow = N) # mean of person-specific intercepts and AR coefficients
b = matrix(rep(0,N*6),nrow = N) # person-specific intercepts and AR coefficients

for(i in 1:N){
  #values at t=1
  Y[i,1,] = rnorm(D,0,1) 
  
  for (d in 1:D){
    muInter[i,d] = sum(CoeffInter[,d]*AforInter[i,])
    for(d2 in 1:D){
      muAR[i, d, d2] = sum(CoeffAR[,d,d2]*AforAR[i,])
    }
  }
  
  bmu[i,1:2] = muInter[i,1:2]
  bmu[i,3] = muAR[i,1,1] #a1
  bmu[i,4] = muAR[i,2,2] #a2
  bmu[i,5] = muAR[i,1,2] #b1
  bmu[i,6] = muAR[i,2,1] #b2
  if(high==0){
    b[i,1:6] = mvrnorm(1,bmu[i,1:6],bcov)
  }
  if(high==1){
    b[i,1:6] = rtmvnorm(1,bmu[i,1:6],bcov,lower=c(rep(-Inf,2),rep(-0.9,4)),upper=c(rep(Inf,2),rep(0.9,4)))
  }
  intercept[i, 1:2] = b[i,1:2]
  AR[i,1,1] = b[i,3] #a1
  AR[i,2,2] = b[i,4] #a2
  AR[i,1,2] = b[i,5] #b1
  AR[i,2,1] = b[i,6] #b2

  for(t in 2:O){
    mus[i,t,] = intercept[i,] + AR[i, ,] %*% (Y[i,t-1,] - intercept[i,])
    Y[i,t,] = mvrnorm(n = 1, mu = mus[i,t,], Sigma = Sigma_VAR[,,i])
  }
}


# drop the first 100 time points
Y = Y[,101:O,]
y1 = Y[,,1]
y2 = Y[,,2]

O=O-100

# check stationarity
plot(y1[1,],type="l",ylim=c(min(y1),max(y1)))
for(i in 2:N){lines(y1[i,])}
plot(y2[1,],type="l",ylim=c(min(y2),max(y2)))
for(i in 2:N){lines(y2[i,])}


################## missing data generation ##############################################
missing_gen<-function(phi0,nmarphi11,nmarphi12,nmarphi21,nmarphi22){
  n=N
  nt=O
  logit1=matrix(rep(NA,(nt*n)),nrow=n)
  pr1=matrix(rep(NA,(nt*n)),nrow=n)
  r1=matrix(rep(NA,(nt*n)),nrow=n)
  
  logit2=matrix(rep(NA,(nt*n)),nrow=n)
  pr2=matrix(rep(NA,(nt*n)),nrow=n)
  r2=matrix(rep(NA,(nt*n)),nrow=n)
  
  for (i in 1:n){
    for (t in 2:nt){
      
      logit1[i,t]=phi0+nmarphi11*y1[i,t]+nmarphi12*y2[i,t]
      pr1[i,t]=exp(logit1[i,t])/(1+exp(logit1[i,t]))
      r1[i,t]=rbinom(1,1,pr1[i,t])
      
      logit2[i,t]=phi0+nmarphi21*y1[i,t]+nmarphi22*y2[i,t]
      pr2[i,t]=exp(logit2[i,t])/(1+exp(logit2[i,t]))
      r2[i,t]=rbinom(1,1,pr2[i,t])
    }#end of t loop
    
    r1[i,1]<-0
    r2[i,1]<-0
  }#end of i loop
  return(list(r1=r1,r2=r2))
}

nmar<-missing_gen(phi0=-1.1,nmarphi11=.6,nmarphi12=.6,nmarphi21=.6,nmarphi22=.6)


print(sum(nmar[[1]])/(N*O))
print(sum(nmar[[2]])/(N*O))

nmarmissy1<-y1
nmarmissy1[nmar[[1]]==1]<-NA
nmarmissy2<-y2
nmarmissy2[nmar[[2]]==1]<-NA

Y_obs = array(rep(NA,N*O*D), c(N,O,D))
Y_obs[,,1] = nmarmissy1
Y_obs[,,2] = nmarmissy2



############################### reformat data for Mplus ################################
cov1 = AforInter[,2]

y1 = Y_obs[,,1]
y2 = Y_obs[,,2]
y1[is.na(y1)]=99
y2[is.na(y2)]=99

mplusdata = matrix(NA,nrow = N*O, ncol = 4)
for(i in 1:N){
  mplusdata[((i-1)*O+1):(i*O),1] = rep(i,O)
  mplusdata[((i-1)*O+1):(i*O),2] = rep(cov1[i],O)
  mplusdata[((i-1)*O+1):(i*O),3] = y1[i,]
  mplusdata[((i-1)*O+1):(i*O),4] = y2[i,]
}


## save simulated data and parameters as .Rdata for Stan and JAGS
## save simulated data as .dat for Mplus
if(high==0){
  save(N, O, K, D, Y_obs, 
       AforAR, Aforcorr_VAR, AforInter, AforlogSD_VAR,
       CoeffAR, CoeffInter, Coeffcorr_VAR, CoefflogSD_VAR,
       sigmaInter, sigmaAR, R_SD, logS_SD,
       AR, intercept, R, S, Sigma_VAR, 
       muInter, muAR, logS_Means, R_Means,
       file = paste0("SimulatedData_low_N",N,"T",O,"_",r1,".Rdata"))
  write.table(mplusdata, file = paste0("SimulatedData_low_N",N,"T",O,"_",r1,".dat"), row.names=FALSE, col.names = FALSE, sep = "\t")
}

if(high==1){
  save(N, O, K, D, Y_obs, 
       AforAR, Aforcorr_VAR, AforInter, AforlogSD_VAR,
       CoeffAR, CoeffInter, Coeffcorr_VAR, CoefflogSD_VAR,
       sigmaInter, sigmaAR, R_SD, logS_SD,
       AR, intercept, R, S, Sigma_VAR, 
       muInter, muAR, logS_Means, R_Means, 
       Tmiss, Tseen, nmiss, nseen,
       file = paste0("SimulatedData_high_N",N,"T",O,"_",r1,".Rdata"))
  write.table(mplusdata, file = paste0("SimulatedData_high_N",N,"T",O,"_",r1,".dat"), row.names=FALSE, col.names = FALSE, sep = "\t")
}

