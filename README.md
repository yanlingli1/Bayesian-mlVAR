# README

This README file lists and describes the simulation materials and other reproducible code for fitting multilevel vector autoregressive (mlVAR) models in Stan, JAGS, and Mplus.
[citation]

The "MAR" folder contains R code used in our simulation study in which missing data geneation followed the missing at random (MAR) mechanism. The simulation materials include:
- *MLVAR_Data_Generation_MAR.R*: code for simulating data based on an mlVAR model and generating missing data following the MAR mechanism
  - *SimulatedData_low_N100T60_1.Rdata*: a simulated dataset used in model fitting with Stan and JAGS
  - *SimulatedData_low_N100T60_1.dat*: a simulated dataset used in model fitting with Mplus
- *MLVAR_Stan_Code.R*: code for fitting mlVAR with MAR missingness in Stan
  - *mlvar.stan*: Stan model script saved in a .stan file
- *MLVAR_JAGS_Code.R*: code for fitting mlVAR with MAR missingness in JAGS
  - *mlvar_4dim.txt*: JAGS model script saved in a .txt file 
- *MLVAR_Mplus_Code.R*: code for fitting mlVAR with MAR missingness in Mplus (via the `MplusAutomation` package) and implementing parameter transformations
- *postcalc.R*: code for calculating summary statistics (e.g., means, medians, modes, standard deviations, credible intervals, effective sample sizes, and Rhat statistics) of the posterior distributions

Due to space constraints, our simulation study only focused on the MAR condition. In this repo, we also provide code for fitting mlVAR with not missing at random (NMAR) missingness in Stan and JAGS. Specifically, the "NMAR" folder contains:
- *MLVAR_Data_Generation_NMAR.R*: code for simulating data based on an mlVAR model and generating missing data following the NMAR mechanism
  - *SimulatedData_low_N100T60_1.Rdata*: a simulated dataset used in model fitting with Stan and JAGS
  - *SimulatedData_low_N100T60_1.dat*: a simulated dataset used in model fitting with Mplus
- *MLVAR_Stan_Code.R*: code for fitting mlVAR with NMAR missingness in Stan
  - *mlvar.stan*: Stan model script saved in a .stan file
- *MLVAR_JAGS_Code.R*: code for fitting mlVAR with NMAR missingness in JAGS
  - *mlvar_9dim.txt*: JAGS model script saved in a .txt file (*Note*: In our simulation study, a subset of random effect covariances were fixed to zeros in the JAGS model (details and reasons can be found in the Simulation Design section in our paper). For didactic reasons, we also provide a version of JAGS scripts where all random effect covariances are freely estimated. However, it may not guarantee model convergence and the sampling efficiency can be quite low.
- *postcalc.R*: code for calculating summary statistics (e.g., means, medians, modes, standard deviations, credible intervals, effective sample sizes and Rhat statistics) of the posterior distributions


*Note*: 
1. Our code can be easily customized for fitting less complex models such as mlVAR with fixed auto- and cross-regression parameters or a fixed innovation covariance matrix. Please contact Yanling Li (yxl823@psu.edu) if you have any questions.
2. When simulating time series data, it is suggested that users plot the time series to check stationarity. 