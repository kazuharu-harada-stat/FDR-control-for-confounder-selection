source("./src/functions_simulation.R")
source("./src/utils.R")

##### COMMON SETTINGs #######################
ls_settings = list(
  n_sim = 100,
  sim_label = "ATE_simple_low", 
  seed = 1234,
  n_cores = 16,
  make_plot=FALSE,
  all_mirrors=TRUE,
  eval_ATE = TRUE,
  
  coef_type = "constant",
  family_Y = "binomial",
  
  n_obs = 1000, # Number of observations
  dim_X = 90,  # Number of predictors included in model
  
  q_YA = 15,  # Number of true common predictors
  q_Y = 15,  # Number of true predictors only for Y
  q_A = 15,  # Number of true predictors only for A

  X_dist = "gaussian",
  corr_type = "Toeplitz_block",
  block_size = 5,
  rho = 0.3, 
  
  true_idx = "fixed",
  beta_XY = 0.10, # beta is divided by sqrt(n_obs)
  beta_XA = 0.20, # beta is divided by sqrt(n_obs)
  effsize = 0.20, # effect size is divided by sqrt(n_obs)
  
  tgt_FDR = 0.1
)
######################################
# list2env(ls_settings, env=environment())

beta_con <- c(0.05,0.10,0.15)
beta_bin <- c(0.10,0.20,0.30)

# 比較実験（低次元）

# changing_rho
for (family_Y in c("gaussian", "binomial")) {
  for (X_dist in c("gaussian", "binomial")) {
    for (eff_flag in c(1)) {
      for (rho in c(0.3)) {
        for (beta_idx in c(2)) {
          ls_settings$family_Y <- family_Y
          ls_settings$X_dist <- X_dist
          ls_settings$rho <- rho
          ls_settings$beta_XY <- ifelse(family_Y=="gaussian", beta_con[beta_idx], beta_bin[beta_idx])
          ls_settings$beta_XA <- beta_bin[beta_idx]
          ls_settings$effsize <- ifelse(family_Y=="gaussian", beta_con[beta_idx], beta_bin[beta_idx])*eff_flag*2
          
          if(family_Y=="binomial"){
            ls_settings$true_ATE <- do.call("simulate_ATE_for_bin",c(list(n_obs_sim=1000000),ls_settings))
          } else {
            ls_settings$true_ATE <- ls_settings$effsize
          }
          gc()
          
          do.call(simulate_all, c(list(algorithm = "DS", method="glm"), ls_settings))
          do.call(simulate_all, c(list(algorithm = "MDS", method="glm"), ls_settings))
          do.call(simulate_all, c(list(algorithm = "Qval", method="glm"), ls_settings))
           
          do.call(simulate_all, c(list(algorithm = "DS", method="cross"), ls_settings))
          do.call(simulate_all, c(list(algorithm = "MDS", method="cross"), ls_settings))
          do.call(simulate_all, c(list(algorithm = "Qval", method="cross"), ls_settings))
          
          do.call(simulate_all, c(list(algorithm = "DS", method="glmnet"), ls_settings))
          do.call(simulate_all, c(list(algorithm = "MDS", method="glmnet"), ls_settings))
          do.call(simulate_all, c(list(algorithm = "Qval", method="marginal"), ls_settings))
          
          do.call(simulate_all, c(list(algorithm = "ALL", method=NA), ls_settings))
        }
      }
    }
  }
}








