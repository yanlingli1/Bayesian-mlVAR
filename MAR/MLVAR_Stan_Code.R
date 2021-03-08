library(rstan)
library(coda)
source("postcalc.R")

#rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

r1=1
N=100
O=60
high=0


if(high==0){
  load(paste0("SimulatedData_low_N",N,"T",O,"_",r1,".Rdata"))
}

if(high==1){
  load(paste0("SimulatedData_high_N",N,"T",O,"_",r1,".Rdata"))
}

Y = Y_obs
Y_miss = which(is.na(Y), arr.ind=TRUE)
Y_miss1=Y_miss[,1]
Y_miss2=Y_miss[,2]
Y_miss3=Y_miss[,3]

y1=Y[,,1]
y2=Y[,,2]
y1[is.na(y1)] = 99
y2[is.na(y2)] = 99
Y_obs[,,1] = y1
Y_obs[,,2] = y2

# now Y stores original data; Y_obs stores data with missing flags


# the scale matrix of the IW distribution
W = diag(9)
# prior for the 1st observation
mus1 = c(0,0)
sigma_VAR_1 = diag(2)


stan_data <- list(D = dim(Y)[3],
                  Y = Y_obs,
                  AforCov = AforInter,
                  nrCoeff = dim(AforInter)[2],
                  nrCoeffInter= dim(AforInter)[2], 
                  nrCoeffAR = dim(AforAR)[2],
                  nrCoefflogSD_VAR = dim(AforlogSD_VAR)[2],
                  nrCoeffcorr_VAR = dim(Aforcorr_VAR)[2],
                  T = dim(Y)[2],
                  P = dim(Y)[1],
                  Y_miss1 = Y_miss1,
                  Y_miss2 = Y_miss2,
                  Y_miss3 = Y_miss3,
                  N_miss = dim(Y_miss)[1],
                  W = W,
                  nrW = dim(W)[1],
                  mus1 = mus1,
                  sigma_VAR_1 = sigma_VAR_1)


inits1 = list(sigma = rep(0.5,9), L = diag(9))
inits2 = list(sigma = rep(0.3,9), L = diag(9))


VARmodel<-stan(file='mlvar.stan', data=stan_data, chains=2, iter=25000, warmup=5000,
               seed=r1, init = list(inits1, inits2), cores = 2, 
               pars=c("CoeffInter","CoeffAR","CoefflogSD_VAR","Coeffcorr_VAR",
                      "sigmaInter","sigmaAR","logS_SD","R_SD","bcov",
                      "intercept", "AR", "corr_noise", "Sigma_VAR"))


codaSamples = mcmc.list(lapply(1:ncol(VARmodel) , function(x) {
  mcmc(as.array(VARmodel)[, x, ])
}))


resulttable <- zcalc(codaSamples)



warns = warnings()



## save results
if(high==0){
  save(warns, resulttable, file=paste0("Result_Stan_low_N",N,"T",O,"_",r1,".Rdata"))
}

if(high==1){
  save(warns, resulttable, file=paste0("Result_Stan_high_N",N,"T",O,"_",r1,".Rdata"))
}

