
data {
int<lower=1> T;         // Length of time-series
int<lower=1> P;         // Number of persons
int<lower=1> N_miss;
int<lower=1> D;         // Dimensions of time-series
vector[D] Y[P,T];       // Data
int R[P,T,D]; // missing indicators
int Y_miss1[N_miss];
int Y_miss2[N_miss];
int Y_miss3[N_miss];
int<lower=1> nrW;
matrix[nrW,nrW] W;
matrix[D,D] sigma_VAR_1;
vector[D] mus1;


int<lower=1> nrCoeff;
int<lower=1> nrCoeffInter;  //Number of covariates for Inter
int<lower=1> nrCoeffAR;     //Number of covariates for AR
int<lower=1> nrCoefflogSD_VAR;   //Number of covariates for sd_noise
int<lower=1> nrCoeffcorr_VAR;    //Number of covariates for fisher_corr_noise

matrix[P,nrCoeff] AforCov; //Covariates for Inter

}


parameters {
vector[N_miss] Y_impute;

matrix[nrCoeff,9] Coeff; //coefficients of covariates
matrix[P,nrW] b;
vector<lower=0>[nrW] sigma;
cholesky_factor_corr[nrW] L;
matrix[3,D] lambda;// regression coefficients in the missing data model
}


transformed parameters {
vector[D] Y_merge[P,T];


matrix[P,nrW] bmu;

vector[D] sd_noise[P];   // prediction noise covariance matrix: standard deviations
vector<lower=-1,upper=1>[P] corr_noise; // prediction noise covariance matrix: correlations
matrix[D,D] Sigma_VAR[P];                    // prediction noise covariance matrix


Y_merge = Y;
for(u in 1:N_miss){
Y_merge[Y_miss1[u],Y_miss2[u],Y_miss3[u]] = Y_impute[u];
}

for (pp in 1:P) {
bmu[pp] = AforCov[pp,]*Coeff;


// Making prediction error covariance matrix for VAR
corr_noise[pp] = (exp(2*b[pp,9])-1)/ (exp(2*b[pp,9])+1);
sd_noise[pp,1] = exp(b[pp,7]);
sd_noise[pp,2] = exp(b[pp,8]);
Sigma_VAR[pp,1,1] = sd_noise[pp, 1] * sd_noise[pp, 1];
Sigma_VAR[pp,2,2] = sd_noise[pp, 2] * sd_noise[pp, 2];
Sigma_VAR[pp,1,2] = sd_noise[pp, 1] * corr_noise[pp] * sd_noise[pp,2]; 
Sigma_VAR[pp,2,1] = Sigma_VAR[pp,1,2];


} //close loop over persons
}



model {
vector[D] mus[P,T];
vector[D] logodds[P,T]; 
vector[D] pr[P,T]; // probability of missingness

// priors on hyperparameters: standard deviations of random parameters
L ~ lkj_corr_cholesky(2);
sigma ~ cauchy(0, 2);

// priors on hyperparameters: coefficients of covariates for random parameters
for(kk in 1:nrCoeff){
for (dd in 1:9){
Coeff[kk,dd] ~ normal(0,1);
}
}

// priors on parameters in the missing data model
for(kk in 1:3){
for(dd in 1:D){
lambda[kk,dd] ~ normal(0,1);
}}

for (pp in 1:P){
// Specify hyperpriors for multilevel structure

b[pp,1:nrW] ~ multi_normal_cholesky(bmu[pp,1:nrW],diag_pre_multiply(sigma, L)); 

//Y_merge[pp, 1,1:D] ~ multi_normal(b[pp,1:2],Sigma_VAR[pp]);
Y_merge[pp,1,1:D] ~ multi_normal(mus1, sigma_VAR_1);

for (tt in 2:T){

mus[pp,tt,1] = b[pp,1] + b[pp,3] * (Y_merge[pp,tt-1,1] - b[pp,1]) + b[pp,5] * (Y_merge[pp,tt-1,2] - b[pp,2]); 
mus[pp,tt,2] = b[pp,2] + b[pp,6] * (Y_merge[pp,tt-1,1] - b[pp,1]) + b[pp,4] * (Y_merge[pp,tt-1,2] - b[pp,2]); 
Y_merge[pp,tt] ~ multi_normal(mus[pp,tt], Sigma_VAR[pp]); 

// missing data model
for(dd in 1:D){
logodds[pp,tt,dd] = lambda[1,dd]+lambda[2,dd]*Y_merge[pp,tt,1]+
                     lambda[3,dd]*Y_merge[pp,tt,2];
pr[pp,tt,dd] = exp(logodds[pp,tt,dd])/(1+exp(logodds[pp,tt,dd]));
R[pp,tt,dd] ~ bernoulli(pr[pp,tt,dd]);
}

} //close loop over timepoints
} //close loop over persons


} 

generated quantities {

vector[D] intercept[P];         //random intercepts
matrix[D,D] AR[P];  //random auto- and cross-regression parameters

matrix[nrW,nrW] bcorr;
matrix[nrW,nrW] bcov;

real<lower=0> R_SD;
vector<lower=0>[D] logS_SD; 
matrix<lower=0>[D,D] sigmaAR;
vector<lower=0>[D] sigmaInter;

matrix[nrCoeffInter,D] CoeffInter; //coefficients of covariates for Inter
matrix[D,D] CoeffAR[nrCoeffAR];   //coefficients of covariates for AR
matrix[nrCoefflogSD_VAR,D] CoefflogSD_VAR;  //coefficients of covariates for sd_noise
matrix[nrCoeffcorr_VAR,1] Coeffcorr_VAR;

for(pp in 1:P){
intercept[pp,1] = b[pp,1];
intercept[pp,2] = b[pp,2];
AR[pp,1,1] = b[pp,3];
AR[pp,2,2] = b[pp,4];
AR[pp,1,2] = b[pp,5];
AR[pp,2,1] = b[pp,6];
}


CoeffInter = Coeff[,1:2];
CoefflogSD_VAR = Coeff[,7:8];
Coeffcorr_VAR[,1] = Coeff[,9];
for(kk in 1:nrCoeffAR){
CoeffAR[kk,1,1] = Coeff[kk,3];
CoeffAR[kk,2,2] = Coeff[kk,4];
CoeffAR[kk,1,2] = Coeff[kk,5];
CoeffAR[kk,2,1] = Coeff[kk,6];
}


bcorr = multiply_lower_tri_self_transpose(L);
bcov = quad_form_diag(bcorr, sigma);

sigmaInter[1] = sqrt(bcov[1,1]);
sigmaInter[2] = sqrt(bcov[2,2]);
sigmaAR[1,1] = sqrt(bcov[3,3]);
sigmaAR[2,2] = sqrt(bcov[4,4]);
sigmaAR[1,2] = sqrt(bcov[5,5]);
sigmaAR[2,1] = sqrt(bcov[6,6]);
logS_SD[1] = sqrt(bcov[7,7]);
logS_SD[2] = sqrt(bcov[8,8]);
R_SD = sqrt(bcov[9,9]);
}

