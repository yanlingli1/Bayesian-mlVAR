library(rjags)
source("postcalc.R")

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

# missing indicators
Index = matrix(NA,nrow=N,ncol=O)
for(pp in 1:N){
  for (oo in 1:O){
    Index[pp,oo] = ifelse(sum(is.na(Y_obs[pp,oo,]))==0, 0, 1)
  }
}

Tmiss = matrix(NA,nrow=N,ncol=O)
Tseen = matrix(NA,nrow=N,ncol=O)
nmiss=c()
nseen=c()
for(pp in 1:N){
  nmiss[pp] = length(which(Index[pp,] %in% 1))
  if(nmiss[pp] != 0){
    Tmiss[pp, 1:nmiss[pp]] <- which(Index[pp,] %in% 1) #location of records that have at least one missing entry
  }
  nseen[pp] = length(which(Index[pp,] %in% 0))
  Tseen[pp, 1:nseen[pp]] <- which(Index[pp,] %in% 0)  #location of fully observed data
}

# the scale matrix of the IW distribution
W = diag(c(0.5,0.5,0.1,0.1),4)
# prior for the 1st observation
mus1 = c(0,0)
prec1 = diag(2)

############################ JAGS Codes################################################

jags_data <- list(Y = Y_obs,
             AforInter = AforInter, 
             AforAR = AforAR, 
             AforlogSD_VAR = AforlogSD_VAR,
             Aforcorr_VAR = Aforcorr_VAR,
             nrCoeffcorr_VAR = dim(Aforcorr_VAR)[2],
             nrCoeffInter= dim(AforInter)[2], 
             nrCoeffAR = dim(AforAR)[2], 
             nrCoefflogSD_VAR = dim(AforlogSD_VAR)[2],
             P = dim(Y_obs)[1], 
             D = dim(Y_obs)[3], 
             nmiss = nmiss, 
             nseen = nseen, 
             Tmiss = Tmiss, 
             Tseen = Tseen,
             W = W, mus1=mus1, prec1=prec1,
             nrW = dim(W)[1])

inits1 <- list(.RNG.name="base::Wichmann-Hill", .RNG.seed=r1)
inits2 <- list(.RNG.name="base::Wichmann-Hill", .RNG.seed=r1+500)


jagsModel <- jags.model(file = "mlvar_4dim.txt", data = jags_data, inits=list(inits1,inits2), n.chains = 2, n.adapt = 4000) #n.adapt = 4000
update(jagsModel, n.iter = 1000) 


parameterlist <-c("CoeffInter","CoeffAR","CoefflogSD_VAR","Coeffcorr_VAR",
                  "sigmaInter","sigmaAR","logS_SD","R_SD","bcov",
                  "intercept", "AR", "corr_noise", "Sigma_VAR")
codaSamples <- coda.samples(jagsModel, variable.names = parameterlist, n.iter = 20000, thin = 1) 


resulttable <- zcalc(codaSamples)

warns = warnings()


## save results
if(high==0){
save(warns, resulttable, file=paste0("Result_JAGS_low_N",N,"T",O,"_",r1,".Rdata"))
}

if(high==1){
save(warns, resulttable, file=paste0("Result_JAGS_high_N",N,"T",O,"_",r1,".Rdata"))
}



