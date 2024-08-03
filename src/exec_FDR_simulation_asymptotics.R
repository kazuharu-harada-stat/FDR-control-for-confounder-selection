source("./src/functions_simulation.R")

##### COMMON SETTINGs #######################
ls_settings = list(
  n_sim = 100,
  sim_label = "asymptotics", 
  seed = 1234,
  n_cores = 16,
  make_plot=FALSE,
  all_mirrors=TRUE,
  
  coef_type = "constant",
  family_Y = "gaussian",
  
  n_obs = 1000, # Number of observations
  dim_X = 100,  # Number of predictors included in model
  
  q_YA = 15,  # Number of true common predictors
  q_Y = 15,  # Number of true predictors only for Y
  q_A = 15,  # Number of true predictors only for A

  X_dist = "gaussian",
  corr_type = "Toeplitz_all",
  rho = 0.3, # correlation
  
  beta_XY = 0.08, # beta is divided by sqrt(n_obs)
  beta_XA = 0.16, # beta is divided by sqrt(n_obs)
  effsize = 0.08, # effect size is divided by sqrt(n_obs)

  tgt_FDR = 0.1
)
######################################
# list2env(ls_settings, env=environment())


# 漸近的性質の確認

for (family_Y in c("gaussian", "binomial")) {
  ls_settings$beta_XY <- ifelse(family_Y=="gaussian", 0.08, 0.16)
  ls_settings$effsize <- ifelse(family_Y=="gaussian", 0.08, 0.16)
  
  for (n_obs in c(500, 1000, 2000, 4000, 8000)) {
    ls_settings$family_Y <- family_Y
    ls_settings$n_obs <- n_obs
    
    do.call(simulate_all, c(list(algorithm = "DS", method="glm"), ls_settings))
    do.call(simulate_all, c(list(algorithm = "MDS", method="glm"), ls_settings))
    do.call(simulate_all, c(list(algorithm = "Qval", method="glm"), ls_settings))
    
    do.call(simulate_all, c(list(algorithm = "DS", method="cross"), ls_settings))
    do.call(simulate_all, c(list(algorithm = "MDS", method="cross"), ls_settings))
    do.call(simulate_all, c(list(algorithm = "Qval", method="cross"), ls_settings))
    
    do.call(simulate_all, c(list(algorithm = "DS", method="glmnet"), ls_settings))
    do.call(simulate_all, c(list(algorithm = "MDS", method="glmnet"), ls_settings))
    do.call(simulate_all, c(list(algorithm = "Qval", method="marginal"), ls_settings))
  }
}











