require(tidyverse)
require(MASS)
require(doParallel)
require(doRNG)
require(mvtnorm)


generate_data <- function(n_obs,
                          dim_X,X_dist,corr_type,block_size,rho,
                          family_Y,coef_type,true_idx,effsize,beta_XY,beta_XA,q_YA,q_Y,q_A,...){
  
  dim_X_pre <- dim_X
  dim_X <- ifelse(dim_X < 150, 150, dim_X)
  
  # Generate X
  if (corr_type == "Toeplitz_all"){
    block_size  <-  dim_X
  } else if (corr_type == "Toeplitz_block") {
    block_size  <-  block_size
  }
  
  num_blocks <- dim_X %/% block_size + 1
  
  times <- 1:block_size
  H <- abs(outer(times, times, "-"))
  V <- rho^H
  X <- foreach(j=1:num_blocks, .combine = cbind)%do%{
    gen_corr_matrix(n_obs, block_size, V, dist=X_dist)
  }
  X <- X[,1:dim_X]
  
  # Generate y
  if(true_idx=="fixed"){
    if (coef_type == "low-dim"){
      beta_Y_vec = beta_XY*sqrt(1/n_obs)
      beta_A_vec = beta_XA*sqrt(1/n_obs)
      effsize_scaled <- effsize*sqrt(1/n_obs)
    } else if (coef_type == "high-dim") {
      beta_Y_vec = beta_XY*sqrt(log(dim_X)/n_obs)
      beta_A_vec = beta_XA*sqrt(log(dim_X)/n_obs)
      effsize_scaled <- effsize*sqrt(log(dim_X)/n_obs)
    } else if (coef_type == "constant") {
      ml <- 1.25
      dg <- 0.75
      sgn_ptrn_Y <- c(ml,1,1,dg,
                      ml,1,dg,
                      1,1,1,dg,
                      -1,-1,-1,-dg)
      sgn_ptrn_A <- c(ml,1,1,dg,
                      ml,1,dg,
                      -1,-1,-1,dg,
                      1,1,1,dg)
      beta_Y_vec = beta_XY*sgn_ptrn_Y
      beta_A_vec = beta_XA*sgn_ptrn_A
      
      effsize_scaled <- effsize
    }
    
    # true_idx_OR <- 1:(q_YA + q_Y + q_A)
    if(block_size == 5){
      true_idx_rel <- c(1,2,3,4,6,11,16,21,22,26,31,36,37,41,46,
                       51,52,53,54,56,61,66,71,72,76,81,86,87,91,96,
                       101,102,103,104,106,111,116,121,122,126,131,136,137,141,146)
    } else if (block_size == 4){
      true_idx_rel <- c(1,2,3,4,5,9,13,17,18,21,25,29,30,33,37,
                       41,42,43,44,45,49,53,57,58,61,65,69,70,73,77,
                       81,82,83,84,85,89,93,97,98,101,105,109,110,113,117)
    }
    
    null_idx <- sample(setdiff(1:dim_X, true_idx_rel), dim_X_pre - 45, replace=FALSE)
    # heatmap(cor(X[,c(true_idx_rel, null_idx)]), Rowv=NA,Colv=NA)
    X <- X[,c(true_idx_rel, null_idx)]
    
    true_idx_OR = c(1:45)
    true_idx_Y = sort(true_idx_OR[1:(q_YA + q_Y)])
    true_idx_A = sort(true_idx_OR[c(1:q_YA, (q_YA + q_A + 1):length(true_idx_OR))])
    true_idx_AND <- intersect(true_idx_Y,true_idx_A)
    true_idx_OnlyY <- setdiff(true_idx_Y,true_idx_AND)
    true_idx_OnlyA <- setdiff(true_idx_A,true_idx_AND)
    
  } else if(true_idx=="random"){
    # Generate y
    if (coef_type == "low-dim"){
      beta_Y_vec = sample(c(1,-1), q_YA + q_Y, replace=T)*beta_XY*sqrt(1/n_obs)
      beta_A_vec = sample(c(1,-1), q_YA + q_A, replace=T)*beta_XA*sqrt(1/n_obs)
      effsize_scaled <- effsize*sqrt(1/n_obs)
    } else if (coef_type == "high-dim") {
      beta_Y_vec = sample(c(1,-1), q_YA + q_Y, replace=T)*beta_XY*sqrt(log(dim_X)/n_obs)
      beta_A_vec = sample(c(1,-1), q_YA + q_A, replace=T)*beta_XA*sqrt(log(dim_X)/n_obs)
      effsize_scaled <- effsize*sqrt(log(dim_X)/n_obs)
    } else if (coef_type == "constant") {
      beta_Y_vec = sample(c(1,-1), q_YA + q_Y, replace=T)*beta_XY
      beta_A_vec = sample(c(1,-1), q_YA + q_A, replace=T)*beta_XA
      effsize_scaled <- effsize
    }
    
    true_idx_OR <- sample(dim_X, q_YA + q_Y + q_A)
    true_idx_Y = sort(true_idx_OR[1:(q_YA + q_Y)])
    true_idx_A = sort(true_idx_OR[c(1:q_YA, (q_YA + q_A + 1):length(true_idx_OR))])
    true_idx_AND <- intersect(true_idx_Y,true_idx_A)
    true_idx_OnlyY <- setdiff(true_idx_Y,true_idx_AND)
    true_idx_OnlyA <- setdiff(true_idx_A,true_idx_AND)
  }
  
  beta_for_Y = numeric(dim_X_pre)
  beta_for_Y[true_idx_Y] = beta_Y_vec
  
  beta_for_A = numeric(dim_X_pre)
  beta_for_A[true_idx_A] = beta_A_vec
  
  true_PS <- sigmoid(X %*% beta_for_A)
  
  A = rbinom(n_obs, 1, true_PS)
  if(family_Y == "gaussian"){
    y = as.vector(effsize_scaled*A + X %*% beta_for_Y + rnorm(n_obs))
  } else if (family_Y == "binomial"){
    y = rbinom(n_obs, 1, sigmoid(effsize_scaled*A + X %*% beta_for_Y))
  }
  
  ls_ret <- list(
    X=X, A=A, y=y, true_PS=true_PS,
    beta_for_Y=beta_for_Y, beta_for_A=beta_for_A, effsize_scaled=effsize_scaled,
    true_idx_OR=true_idx_OR,true_idx_Y=true_idx_Y,true_idx_A=true_idx_A,true_idx_AND=true_idx_AND,
    true_idx_OnlyY=true_idx_OnlyY,true_idx_OnlyA=true_idx_OnlyA
  )
  return(ls_ret)
}

# algorithm="DS";method="glmnet"
simulate <- function(algorithm,method,tgt_FDR,all_mirrors=TRUE,eval_ATE = FALSE,
                     n_obs,dim_X,X_dist,corr_type,block_size,rho,
                     family_Y,coef_type,true_idx,effsize,beta_XY,beta_XA,q_YA,q_Y,q_A,...){
  
  ls_data <- generate_data(n_obs,dim_X,X_dist,corr_type,block_size,rho,
                           family_Y,coef_type,true_idx,effsize,beta_XY,beta_XA,q_YA,q_Y,q_A)
  list2env(ls_data, envir = environment())
  
  # ================================================================
  
  extra_args = list(...)
  
  if(algorithm == "DS"){
    if(all_mirrors == TRUE){
      ls_selected <- SelConf_DS_compare_mirror(X,A,y,family_Y=family_Y,method=method,tgt_FDR=tgt_FDR,n_grid=1000)
    }
  } else if (algorithm == "MDS"){
    if(all_mirrors == TRUE){
      ls_selected <- SelConf_MDS_compare_mirror(X,A,y,family_Y=family_Y,method=method,tgt_FDR=tgt_FDR,n_repeat=100,n_grid=1000)
    }
  } else if (algorithm == "Qval"){
    ls_selected <- SelConf_Qval(X,A,y,family_Y=family_Y,method=method,tgt_FDR=tgt_FDR)
  } else if (algorithm == "ALL"){
    ls_selected <- list(Selected = 1:ncol(X))
  }
  
  if(eval_ATE == FALSE){
    list2longdf <- function(ls_selected, true_idx, result_label){
      tmp <- lapply(ls_selected, FUN= function(item) eval_selection(dim_X, item, true_idx, result_label = result_label))
      df_eval <- do.call(rbind, lapply(names(tmp), function(s) {cbind(method_selection = s, tmp[[s]])}))
      df_ret1 <- df_eval %>%
        dplyr::select("method_selection", "result_label","FDP", "Power", "n_selected") %>% 
        pivot_longer(cols= -c("method_selection", "result_label"), names_to = "measure") %>% mutate(cval=NA)
      df_ret2 <- df_eval %>%
        dplyr::select("method_selection", "result_label","selected_idx","true_idx") %>% 
        pivot_longer(cols= -c("method_selection", "result_label"), names_to = "measure", values_to = "cval") %>% mutate(value=NA)
      df_ret <- rbind(df_ret1,df_ret2)
      return(df_ret)
    }
  } else {
    list2longdf <- function(ls_selected, true_idx, result_label){
      tmp1 <- lapply(ls_selected, FUN= function(item) eval_selection(dim_X, item, true_idx, result_label = result_label))
      tmp2 <- lapply(ls_selected, FUN= function(item) data.frame(estimate_ATE(X[,item], A, y, family=family_Y, result_label = result_label)))
      tmp <- mapply(function(v1,v2) c(v1,v2), tmp1, tmp2, SIMPLIFY = FALSE)
      df_eval <- do.call(rbind, lapply(names(tmp), function(s) {cbind(method_selection = s, data.frame(tmp[[s]]))}))
      df_ret1 <- df_eval %>%
        dplyr::select("method_selection", "result_label","FDP", "Power", "n_selected",
                      "ATE_DIF","ATE_REG","ATE_IPW","ATE_DR","ATE_trim","ATO") %>% 
        pivot_longer(cols= -c("method_selection", "result_label"), names_to = "measure") %>% mutate(cval=NA)
      df_ret2 <- df_eval %>%
        dplyr::select("method_selection", "result_label","selected_idx","true_idx") %>% 
        pivot_longer(cols= -c("method_selection", "result_label"), names_to = "measure", values_to = "cval") %>% mutate(value=NA)
      df_ret <- rbind(df_ret1,df_ret2)
      return(df_ret)
    }
  }
  
  df_eval_Y <- list2longdf(ls_selected, true_idx=true_idx_Y, result_label="Y")
  df_eval_A <- list2longdf(ls_selected, true_idx=true_idx_A, result_label="A")
  df_eval_OR <- list2longdf(ls_selected, true_idx=true_idx_OR, result_label="OR")
  df_eval_AND <- list2longdf(ls_selected, true_idx=true_idx_AND, result_label="AND")
  df_eval_OnlyY <- list2longdf(ls_selected, true_idx=true_idx_OnlyY, result_label="OnlyY")
  df_eval_OnlyA <- list2longdf(ls_selected, true_idx=true_idx_OnlyA, result_label="OnlyA")
  
  df_ret <- rbind(
    df_eval_Y, df_eval_A, df_eval_OR, df_eval_AND, df_eval_OnlyY, df_eval_OnlyA
  ) %>% 
  rename(eval_idx = result_label) %>% 
  separate(col=method_selection, sep="_", into=c("target","method1","method2","method3","method4","method5")) %>% 
  mutate(
    algorithm = algorithm,
    method0 = method,
    tgt_FDR = tgt_FDR,
    n_obs=n_obs,dim_X=dim_X,X_dist=X_dist,corr_type=corr_type,rho=rho,
    family_Y=family_Y,coef_type=coef_type,effsize=effsize,beta_XY=beta_XY,beta_XA=beta_XA,q_YA=q_YA,q_Y=q_Y,q_A=q_A
  )
  
  return(df_ret)
}

make_fname_caption <- function(timestamp,algorithm, method, n_sim, tgt_FDR, seed,
                               n_obs,dim_X,X_dist,corr_type,rho,
                               family_Y,coef_type,effsize,beta_XY,beta_XA,q_YA,q_Y,q_A,...){
  
  current_datetime <- Sys.time()
  formatted_datetime <- format(current_datetime, "%y%m%d%H%M")
  formatted_datetime2 <- format(current_datetime, "%Y-%m-%d %H:%M:%S")
  
  method_short = case_when(
    method == "glm" ~ "glm",
    method == "glmnet" ~ "las",
    method == "scalefree" ~ "scf",
    method == "cross" ~ "crs",
    method == "marginal" ~ "mar",
    method == "debiased" ~ "deb"
  )
  
  fname <- sprintf("%s_%s_Y%s_X%s_s%04d_n%04d_p%04d_q%02d_%02d_%02d_bY%04d_A%04d_e%04d_r%02d_F%02d_seed%04d",
                   toupper(substr(algorithm, 1, 2)), method_short, substr(family_Y, 1,2), substr(X_dist, 1,2),
                   n_sim, n_obs, dim_X, q_YA, q_Y, q_A,
                   round(beta_XY*100), round(beta_XA*100), round(effsize*100), rho*100, tgt_FDR*100, seed)
  
  if(timestamp == TRUE) fname <- paste0(fname, "_", formatted_datetime)
  
  caption_setting <- sprintf("Settings: Algorithm: %s(%s), Outcome Type: %s, dist. of X: %s, # or sim.=%d, n=%d, p=%d, q(common)=%d, q(Y)=%d, q(A)=%d,\n Coef. type: %s, beta(Y)=%.2f, beta(A)=%.3f, Eff. Size=%.3f, Corr. type: %s, rho=%.2f, target FDR=%.2f, random seed=%d, timestamp: %s",
                             algorithm, method, family_Y, X_dist,
                             n_sim, n_obs, dim_X, q_YA, q_Y, q_A,
                             coef_type, beta_XY, beta_XA, effsize, corr_type, rho, tgt_FDR, seed, formatted_datetime2)
  
  ret <- list(fname=fname, caption_setting=caption_setting)
  
  return(ret)
}




simulate_all <- function(n_sim, n_cores=12, seed=1234, sim_label="test", make_plot=FALSE, ...){
  source("./src/utils.R")
  
  extra_args <- list(...)
  
  folder_path <- sprintf("./output/%s", sim_label)
  
  if (!file.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
    dir.create(paste0(folder_path,"/fig"), recursive = TRUE)
  }
  
  fname_cap <- do.call("make_fname_caption", c(list(timestamp=FALSE, n_sim=n_sim, seed=seed), extra_args))
  if (file.exists(paste0(folder_path,"/", fname_cap$fname, ".csv"))){
    print(paste0("Already esists: ", fname_cap$fname))
    return(NULL)
  }
  
  objects_to_export <- c("simulate","generate_data")
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  registerDoRNG(seed)
  
  # for(k in 1:n_sim){
  df_eval_all <- foreach(k=1:n_sim, .combine = rbind, .export=objects_to_export)%dopar%{
    require(MASS)
    require(tidyverse)
    source("./src/utils.R")
    source("./src/DataSplit.R")
    source("./src/q_values_for_Causal.R")
    
    sttime <- Sys.time()
    df_eval <- do.call(simulate, extra_args)
    entime <- Sys.time()
    
    df_eval$sim_idx <- k
    df_eval$exec_time <- entime - sttime 
    
    return(df_eval)
  }
  stopCluster(cl)
  
  df_eval_all$n_sim = n_sim
  df_eval_all$sim_label = sim_label
  df_eval_all$randomseed = seed
  df_eval_all$true_ATE = extra_args$true_ATE
  
  write.csv(df_eval_all, file=paste0(folder_path,"/", fname_cap$fname, ".csv"))
  
  # VISUALIZATION
  items <- c("FDP", "Power","n_selected")
  df_hline <- data.frame(
    measure=factor(c("FDP"), levels=items),
    yint = c(extra_args$tgt_FDR)
  )
  
  df_mean <- df_eval_all %>%
    filter(measure %in% items) %>% 
    group_by(algorithm, method0, method1, method2, target, eval_idx, measure) %>%
    summarize(meanVal = mean(value)) %>% 
    mutate(
      measure = factor(measure, levels=items),
      xval = paste0(target,"_",eval_idx),
      y_facet = paste0(method1,"_",method2)
    )
  
  if (make_plot == TRUE){
    gg <- df_eval_all %>%
      filter(measure %in% items) %>% 
      mutate(
        measure = factor(measure, levels=items),
        xval = paste0(target,"_",eval_idx),
        y_facet = paste0(method1,"_",method2)
      ) %>% 
      ggplot(aes(x=xval, y=value)) +
      geom_boxplot() +
      geom_hline(data = df_hline, aes(yintercept = yint), col="tomato", linetype="dotted", linewidth=1) +
      geom_point(data = df_mean, aes(x=xval, y=meanVal), shape=4, size=2) +
      theme_bw() +
      theme(
        text = element_text(size=10),
        axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5),
        plot.caption = element_text(size = 7)
      ) +
      labs(caption = fname_cap$caption_setting, x="", y="") +
      facet_grid(measure ~ y_facet, scales="free")
    
    ggsave(paste0(folder_path, "/fig/", fname_cap$fname, ".png"), plot = gg,
           width = length(unique(df_mean$y_facet))*35, height = 180, dpi = 500, units = "mm")
  }
}


