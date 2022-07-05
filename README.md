# MVJSDMrealdata
 Multivariate Probit Joint Species Distribution Model (JSDM) for use with real covariate and observation data.
 This model allows for imperfect detection of data.
 - It expects observation data in binary form in a 3D array allowing for multiple survey sites, multiple survey replications, and multiple species
 - It assumes covariates in array format, with 2D covariate data assumed for the occupancy model and 3D covariate data assumed for the detection model.
 - It allows for testing and validation data - and estimates occupancy values for the occupancy dataset. 
 The model is implemented using MCMC analysis from JAGS and requires the jagsUI R package.  Users will need to have installed version of JAGS available.

This JSDM implements and extends the model described in : 
Pollock, L.J., Tingley, R., Morris, W.K., Golding, N., O’Hara, R.B., Parris, K.M., Vesk,
P.A. & McCarthy, M.A. (2014) Understanding co-occurrence by modelling species simultaneously with a joint species distribution model (jsdm). *Methods in Ecology and
Evolution*, **5**, 397–406.

In this project are:
Files to put covariate and observation data in a form acceptable for the JAGS model
Provides a number of Jags models which allow for variations as follows:
- Use the JSDM with or without detection
- Use a Wishart or Cholesky based prior for the covariance matrix
- Address implementation issues relating translating matrices from R to JAGS
  - this requires models to be implement for 0, 1 or multiple occupancy covariates
Runs the Jags models
Provides some basic programming of results to extract parameters, determine rhat and n.eff values for each parameter, and calculate AUC ROC.  
