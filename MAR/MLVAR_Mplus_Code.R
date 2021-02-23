library(MplusAutomation)  # to interact with Mplus
library(coda)   
source('postcalc.R') # to summarize Bayesian estimation results
#source('mplus.r')  # to generate plots

r1=1
N=100
O=60
C="low"

data = read.table(paste0("SimulatedData_",C, "_N",N,"T",O,"_",r1,".dat"))
colnames(data) = c("subject","cov1","y1","y2")
TimePoint = c()
for(pp in 1:N){
  tmpdata = data[data$subject==pp,]
  del = 0
  oo = O
  while(tmpdata$y1[oo] == 99 & tmpdata$y2[oo] == 99){
    del = del + 1
    oo = oo - 1
  }
  TimePoint[pp] = O - del
}

input = mplusObject(TITLE = "A two-level, VAR(1) model for continuous dependent variables 
    with random intercepts and random slopes", 
                    VARIABLE = "MISSING ARE y1 y2 (99);
                      CLUSTER = subject;
                      BETWEEN = cov1; 
                      LAGGED = y1(1) y2(1);",
                    ANALYSIS = "TYPE = TWOLEVEL RANDOM;
                      ESTIMATOR = BAYES;  
                      PROCESSORS = 2;
                      FBITERATIONS = 40000;", 
                    MODEL = "%WITHIN%
                      a1 |y1 ON y1&1; 
                      a2 |y2 ON y2&1; 
                      b1 |y1 ON y2&1; 
                      b2 |y2 ON y1&1; 
                      f1 by y1@1; y1@0.2;
                      f2 by y2@1; y2@0.2;
                      s | f1 on f2;
                      v1 | f1;
                      v2 | f2;	
                      
                      %BETWEEN%
                      y1 y2 a1 a2 b1 b2 s v1 v2 (vv1-vv9);
                      y1 y2 a1 a2 b1 b2 s v1 v2 with y1 y2 a1 a2 b1 b2 s v1 v2 (c1-c36);
                      y1 y2 a1 a2 b1 b2 ON cov1 (s1-s6); 
                      [y1 y2 a1 a2 b1 b2] (s7-s12);
                      ",
                    MODELPRIORS = "s1-s12 ~ N(0,1);
                      vv1-vv2 ~ IW(0.5,12);
                      vv3-vv6 ~ IW(0.1,12);
                      vv7-vv9 ~ IW(0.5,12);
                      c1-c36 ~ IW(0,12);", 
                    OUTPUT = "TECH1 TECH8;",
                    SAVEDATA = paste0("SAVE=FSCORES(1000 1);FILE IS savedata_",C, "_N",N,"T",O,"_",r1,".csv;"),
                    PLOT = "TYPE = PLOT3;
                      FACTOR = ALL;",
                    usevariables = c("subject","cov1","y1","y2"),
                    rdata = data,
                    autov = TRUE)
res = mplusModeler(input, modelout = paste0("mlvar_",C, "_N",N,"T",O,"_",r1,".inp"), run = 1L)


reorder_subject = unique(res$results$savedata$SUBJECT)
AforlogSD_VAR = c()
for(i in 1:N){
  AforlogSD_VAR[i] = unique(data[data$subject==reorder_subject[i],2])
}
Aforcorr_VAR = AforlogSD_VAR

###################obtain posteriors of coefficients of predictors on random parameters
nchain=2
# get parameter labels
# param = mplus.list.bayesian.parameters(input) 
# param = res$results$summaries$Parameters

parameter = c("a_10","a_20","b_10","b_20", #1-4
              "mu_10","mu_20","alpha_a_1", "alpha_a_2","alpha_b_1",
              "alpha_b_2","alpha_mu_1","alpha_mu_2", #8-15
              "sd_a1", #16
              "cov_a1a2", #17
              "sd_a2", #18
              "cov_a1b1", #19
              "cov_a2b1", #20
              "sd_b1", #21
              "cov_a1b2", #22
              "cov_a2b2", #23
              "cov_b1b2", #24
              "sd_b2", #25
              "cov_mu1a1", # 44
              "cov_mu1a2", #45
              "cov_mu1b1", #46
              "cov_mu1b2", #47
              "sd_mu1",#51
              "cov_mu2a1", #52
              "cov_mu2a2",#53
              "cov_mu2b1", #54
              "cov_mu2b2", #55
              "cov_mu12", #59
              "sd_mu2") #60
# get MCMC samples for all parameters
#test = mplus.get.bayesian.parameter.data(input,1:length(param),1:nchain)
test = res$results$gh5$bayesian_data$parameters_autocorr$parameters
test[c(16,18,21,25,51,60),,]=sqrt(test[c(16,18,21,25,51,60),,]) ##transformed to sigmaAR and sigmaInter
test = test[c(1:4,8:25,44:47,51:55,59:60),,]
test=aperm(test,c(2,3,1))
iter=dim(test)[1]
# burn the first half
test=test[(iter/2+1):iter,,]
# add rownames
dimnames(test) = list(NULL, NULL, parameter)
# reformat it to a mcmc.list
codaSamples = mcmc.list(lapply(1:ncol(test) , function(x) {
  mcmc(as.array(test)[, x,])
}))
# use zcalc()
estimate <- zcalc(codaSamples)

# compare to results from Mplus
#result_Mplus=result$parameters$unstandardized


#########################obtain posteriors of random coefficients########################
# get parameter labels
# param1 = mplus.list.bayesian.plausible.labels(input)

# get MCMC samples for all parameters

tmp = res$results$gh5$bayesian_data$plausible$plausible

#ndraw = dim(tmp)[2]/100
ndraw = dim(tmp)[2]
npar = dim(tmp)[3]

fac = array(rep(NA,N*ndraw*npar), dim = c(N,ndraw,npar))

oo = 0
for(i in 1:N){
  fac[i,,] = tmp[oo+1,,]
  oo = oo + TimePoint[reorder_subject[i]]
}

fac_mu1=matrix(NA,nrow=N,ncol=ndraw)
fac_mu2=matrix(NA,nrow=N,ncol=ndraw)
fac_a1=matrix(NA,nrow=N,ncol=ndraw)
fac_a2=matrix(NA,nrow=N,ncol=ndraw)
fac_b1=matrix(NA,nrow=N,ncol=ndraw)
fac_b2=matrix(NA,nrow=N,ncol=ndraw)
fac_s=matrix(NA,nrow=N,ncol=ndraw)
fac_v1=matrix(NA,nrow=N,ncol=ndraw)
fac_v2=matrix(NA,nrow=N,ncol=ndraw)
fac_sd1=matrix(NA,nrow=N,ncol=ndraw)
fac_sd2=matrix(NA,nrow=N,ncol=ndraw)
fac_cov=matrix(NA,nrow=N,ncol=ndraw)
fac_corr=matrix(NA,nrow=N,ncol=ndraw)
fac_z=matrix(NA,nrow=N,ncol=ndraw)

for(i in 1:N){
  fac_a1[i,]=fac[i,,3]
  fac_a2[i,]=fac[i,,4]
  fac_b1[i,]=fac[i,,5]
  fac_b2[i,]=fac[i,,6]
  fac_s[i,]=fac[i,,7]
  fac_v1[i,]=fac[i,,8]
  fac_v2[i,]=fac[i,,9]
  fac_mu1[i,]=fac[i,,10]
  fac_mu2[i,]=fac[i,,11]
  fac_sd1[i,]=sqrt(fac_s[i,]^2*exp(fac_v2[i,])+exp(fac_v1[i,])+0.2)
  fac_sd2[i,]=sqrt(exp(fac_v2[i,])+0.2)
  fac_cov[i,]=fac_s[i,]*exp(fac_v2[i,])
  fac_corr[i,]=fac_cov[i,]/fac_sd1[i,]/fac_sd2[i,]
  fac_z[i,]=log((1+fac_corr[i,])/(1-fac_corr[i,]))/2  #fisher z = 1/2*ln((1+r)/(1-r))
}


parameter_add = c("logS_10","alpha_logS_1","logS_SD_1",
                  "logS_20","alpha_logS_2","logS_SD_2",
                  "corr_0","alpha_corr","R_SD",
                  "bcov[1,7]","bcov[1,8]","bcov[1,9]",
                  "bcov[2,7]","bcov[2,8]","bcov[2,9]",
                  "bcov[3,7]","bcov[3,8]","bcov[3,9]",
                  "bcov[4,7]","bcov[4,8]","bcov[4,9]",
                  "bcov[5,7]","bcov[5,8]","bcov[5,9]",
                  "bcov[6,7]","bcov[6,8]","bcov[6,9]",
                  "bcov[7,8]","bcov[7,9]","bcov[8,9]")

test_add = array(rep(NA,ndraw*1*length(parameter_add)), dim = c(ndraw,1,length(parameter_add)))

dimnames(test_add) = list(NULL, NULL, parameter_add)
for(d in 1:ndraw){
  df_random = data.frame(cbind(fac_mu1[,d],fac_mu2[,d],fac_a1[,d],fac_a2[,d],fac_b1[,d],fac_b2[,d],log(fac_sd1[,d]),log(fac_sd2[,d]),fac_z[,d],AforlogSD_VAR))
  colnames(df_random)=c("mu1","mu2","a1","a2","b1","b2","logsd1","logsd2","z","x1")
  fit_random = lm(cbind(mu1,mu2,a1,a2,b1,b2,logsd1,logsd2,z)~x1,df_random)
  bcov = cov(fit_random$residuals)
  test_add[d,,1] = fit_random$coefficients[1,7] #CoefflogSD_VAR[1,1]
  test_add[d,,2] = fit_random$coefficients[2,7]#CoefflogSD_VAR[2,1]
  test_add[d,,3] = sqrt(bcov[7,7]) #logS_SD[1]
  test_add[d,,4] = fit_random$coefficients[1,8]   #CoefflogSD_VAR[1,2]
  test_add[d,,5] = fit_random$coefficients[2,8]#CoefflogSD_VAR[2,2]
  test_add[d,,6] = sqrt(bcov[8,8]) #logS_SD[2]
  test_add[d,,7] = fit_random$coefficients[1,9] #Coeffcorr_VAR[1,1]
  test_add[d,,8] = fit_random$coefficients[2,9]#CoefflogSD_VAR[2,1]
  test_add[d,,9] = sqrt(bcov[9,9]) #R_SD
  test_add[d,,10] = bcov[1,7] #cov(fac_mu1[,d],log(fac_sd1[,d]))
  test_add[d,,11] = bcov[1,8]#cov(fac_mu1[,d],log(fac_sd2[,d]))
  test_add[d,,12] = bcov[1,9]#cov(fac_mu1[,d],fac_z[,d])
  test_add[d,,13] = bcov[2,7]#cov(fac_mu2[,d],log(fac_sd1[,d]))
  test_add[d,,14] = bcov[2,8]#cov(fac_mu2[,d],log(fac_sd2[,d]))
  test_add[d,,15] = bcov[2,9]#cov(fac_mu2[,d],fac_z[,d])
  test_add[d,,16] = bcov[3,7]#cov(fac_a1[,d],log(fac_sd1[,d]))
  test_add[d,,17] = bcov[3,8]#cov(fac_a1[,d],log(fac_sd2[,d]))
  test_add[d,,18] = bcov[3,9]#cov(fac_a1[,d],fac_z[,d])
  test_add[d,,19] = bcov[4,7]#cov(fac_a2[,d],log(fac_sd1[,d]))
  test_add[d,,20] = bcov[4,8]#cov(fac_a2[,d],log(fac_sd2[,d]))
  test_add[d,,21] = bcov[4,9]#cov(fac_a2[,d],fac_z[,d])
  test_add[d,,22] = bcov[5,7]#cov(fac_b1[,d],log(fac_sd1[,d]))
  test_add[d,,23] = bcov[5,8]#cov(fac_b1[,d],log(fac_sd2[,d]))
  test_add[d,,24] = bcov[5,9]#cov(fac_b1[,d],fac_z[,d])
  test_add[d,,25] = bcov[6,7]#cov(fac_b2[,d],log(fac_sd1[,d]))
  test_add[d,,26] = bcov[6,8]#cov(fac_b2[,d],log(fac_sd2[,d]))
  test_add[d,,27] = bcov[6,9]#cov(fac_b2[,d],fac_z[,d])
  test_add[d,,28] = bcov[7,8]#cov(log(fac_sd1[,d]),log(fac_sd2[,d]))
  test_add[d,,29] = bcov[7,9]#cov(log(fac_sd1[,d]),fac_z[,d])
  test_add[d,,30] = bcov[8,9]#cov(log(fac_sd2[,d]),fac_z[,d])
}


# reformat it to a mcmc.list
codaSamples_add = mcmc.list(lapply(1:ncol(test_add) , function(x) {
  mcmc(as.array(test_add)[, x,])
}))
# use zcalc()
estimate_add <- zcalc(codaSamples_add)
estimate_fixed = rbind(estimate,estimate_add)

rownames(estimate_fixed) = c("CoeffAR[1,1,1]", "CoeffAR[1,2,2]", "CoeffAR[1,1,2]", "CoeffAR[1,2,1]",
                             "CoeffInter[1,1]","CoeffInter[1,2]",
                             "CoeffAR[2,1,1]",
                             "CoeffAR[2,2,2]",
                             "CoeffAR[2,1,2]",
                             "CoeffAR[2,2,1]",
                             "CoeffInter[2,1]",
                             "CoeffInter[2,2]",
                             "sigmaAR[1,1]",
                             "bcov[3,4]","sigmaAR[2,2]",
                             "bcov[3,5]","bcov[4,5]","sigmaAR[1,2]",
                             "bcov[3,6]","bcov[4,6]","bcov[5,6]","sigmaAR[2,1]",
                             "bcov[1,3]","bcov[1,4]","bcov[1,5]","bcov[1,6]",
                             "sigmaInter[1]",
                             "bcov[2,3]","bcov[2,4]","bcov[2,5]","bcov[2,6]",
                             "bcov[1,2]",
                             "sigmaInter[2]",
                             "CoefflogSD_VAR[1,1]","CoefflogSD_VAR[2,1]","logS_SD[1]",
                             "CoefflogSD_VAR[1,2]","CoefflogSD_VAR[2,2]","logS_SD[2]",
                             "Coeffcorr_VAR[1,1]", "Coeffcorr_VAR[2,1]","R_SD",
                             "bcov[1,7]","bcov[1,8]","bcov[1,9]",
                             "bcov[2,7]","bcov[2,8]","bcov[2,9]",
                             "bcov[3,7]","bcov[3,8]","bcov[3,9]",
                             "bcov[4,7]","bcov[4,8]","bcov[4,9]",
                             "bcov[5,7]","bcov[5,8]","bcov[5,9]",
                             "bcov[6,7]","bcov[6,8]","bcov[6,9]",
                             "bcov[7,8]","bcov[7,9]","bcov[8,9]")


estimate_fixed = estimate_fixed[c("CoeffAR[1,1,1]", "CoeffAR[2,1,1]",
                                  "CoeffAR[1,2,1]","CoeffAR[2,2,1]",
                                  "CoeffAR[1,1,2]", "CoeffAR[2,1,2]",
                                  "CoeffAR[1,2,2]", "CoeffAR[2,2,2]",
                                  "CoeffInter[1,1]","CoeffInter[2,1]",
                                  "CoeffInter[1,2]","CoeffInter[2,2]",
                                  "CoefflogSD_VAR[1,1]","CoefflogSD_VAR[2,1]",
                                  "CoefflogSD_VAR[1,2]","CoefflogSD_VAR[2,2]",
                                  "Coeffcorr_VAR[1,1]", "Coeffcorr_VAR[2,1]",
                                  "sigmaInter[1]","sigmaInter[2]",
                                  "sigmaAR[1,1]","sigmaAR[2,1]","sigmaAR[1,2]","sigmaAR[2,2]",
                                  "logS_SD[1]","logS_SD[2]","R_SD",
                                  "bcov[1,2]","bcov[1,3]","bcov[2,4]",
                                  "bcov[1,4]","bcov[1,5]","bcov[1,6]","bcov[1,7]","bcov[1,8]","bcov[1,9]",
                                  "bcov[2,3]","bcov[2,5]","bcov[2,6]","bcov[2,7]","bcov[2,8]","bcov[2,9]",
                                  "bcov[3,4]","bcov[3,5]","bcov[3,6]","bcov[3,7]","bcov[3,8]","bcov[3,9]",
                                  "bcov[4,5]","bcov[4,6]","bcov[4,7]","bcov[4,8]","bcov[4,9]",
                                  "bcov[5,6]","bcov[5,7]","bcov[5,8]","bcov[5,9]",
                                  "bcov[6,7]","bcov[6,8]","bcov[6,9]",
                                  "bcov[7,8]","bcov[7,9]","bcov[8,9]"),]


# Check recovery of person-specific parameters
fac_add = array(rep(NA,N*ndraw*(npar+4)), dim = c(N,ndraw,(npar+4)))
fac_add[,,1:npar] = fac
fac_add[,,npar+1]=fac_sd1^2
fac_add[,,npar+2]=fac_sd2^2
fac_add[,,npar+3]=fac_cov
fac_add[,,npar+4]=fac_corr


param1 = c("F1","F2","AR[,1,1]","AR[,2,2]","AR[,1,2]","AR[,2,1]","S","V1","V2",
           "intercept[,1]","intercept[,2]","Sigma_VAR[,1,1]","Sigma_VAR[,2,2]",
           "Sigma_VAR[,1,2]","corr_noise")
test1 = array(rep(NA,ndraw*N*length(param1)), dim = c(ndraw,1,N*length(param1)))

parameter1 = c()
for(par in 1:length(param1)){
  for(i in 1:N){
    test1[,1,(par-1)*N+i] = fac_add[i,,par]
    parameter1[(par-1)*N+i]=paste0(param1[par],"[",i,"]")
  }
}


dimnames(test1) = list(NULL, NULL, parameter1)
# reformat it to a mcmc.list
codaSamples1 = mcmc.list(lapply(1:ncol(test1) , function(x) {
  mcmc(as.array(test1)[, x,])
}))
estimate1 <- zcalc(codaSamples1)

resulttable = rbind(estimate_fixed, estimate1)
save(resulttable,file = paste0("Result_Mplus_",C, "_N",N,"T",O,"_",r1,".Rdata"))



