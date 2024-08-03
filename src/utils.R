library(caret)
require(xtable)

sigmoid <- function(u){1/(1 + exp(-u))}

logit <- function(p) log(p/(1-p))

weighted_odds = function(y,w){
  return(sum(w*y)/sum(w)/sum(w*(1-y))*sum(w))
}

output_latex <- function(df, fname, file_loc="./output/"){
  latex_table <- xtable(df)
  sink(paste0(file_loc, fname))
  print.xtable(latex_table, include.rownames = FALSE)
  sink()
}


capitalizeFirst <- function(s) {
  paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}


# selected_idx=ls_selected$Y_Marginal_NA_NA_NA_NA; true_idx=true_idx_Y; p=dim_X;result_label="none"
eval_selection <- function(p, selected_idx, true_idx, result_label){
  
  true_judge <- factor(1:p %in% true_idx, levels = c(T, F))
  fit_judge <- factor(1:p %in% selected_idx, levels = c(T, F))
  eval_tab <- table(fit_judge, true_judge)
  eval_results <- confusionMatrix(eval_tab)
  
  df_eval <- data.frame(t(eval_results$byClass))
  df_eval$FDP <- eval_tab[1,2]/max(sum(eval_tab[1,]),1)
  df_eval$Power <- df_eval$Sensitivity
  df_eval$n_selected <- length(selected_idx)
  df_eval$selected_idx <- paste(selected_idx, collapse = "_")
  df_eval$true_idx <- paste(true_idx, collapse = "_")
  df_eval$result_label <- result_label
  
  return(df_eval)
}




Mirror <- function(T1, T2, mirror_type="sum", tgt_FDR=0.1, n_grid=1000){

  sign12 <- sign(T1*T2)

  if(mirror_type == "sum"){
    M_vec <- sign12 * (abs(T1) + abs(T2))
  } else if(mirror_type == "product"){
    M_vec <- sign12 * (abs(T1)*abs(T2))
  } else if(mirror_type == "2min"){
    M_vec <- sign12 * min(abs(T1),abs(T2))*2
  }
  
  ls_t <- seq(0,max(abs(M_vec)),max(abs(M_vec))/n_grid)

  tauq <- calc_threshold(M_vec, ls_t, tgt_FDR)
  S1 <- which(M_vec > tauq)
  
  ls_ret <- list(Selected = S1, Statistic = M_vec, Threshold = tauq)
  
  return(ls_ret)
}

calc_threshold <- function(M_vec, ls_t, tgt_FDR){
  FDP <- function(t){
    return(length(which(M_vec < -t))/max(length(which(M_vec > t)),1))
  }
  
  if(sum(sapply(ls_t, FDP) > tgt_FDR) == 0){
    return(0)
  } else {
    tauq_idx <- max(which(sapply(ls_t, FDP) > tgt_FDR)) + 1
    return(ls_t[tauq_idx])
  }
}


select_confounder_using_mirror <- function(MirrorY, MirrorA, tgt_FDR, n_grid=1000, c_adj=1){
  
  grid_max <- max(sqrt(MirrorY^2 + MirrorA^2))
  grid_cut <- seq(0,grid_max,grid_max/n_grid)
  
  FDP_AND = numeric(length(grid_cut))
  FDP_OR = numeric(length(grid_cut))
  
  for( i in 1:length(grid_cut) ){
    cutY = grid_cut[i]
    cutA = grid_cut[i]
    FDP_OR[i] = sum((MirrorY < -cutY) | (MirrorA < -cutA))/max(sum((MirrorY > cutY) | (MirrorA > cutA)),1)
    
    FD_Y <- sum((MirrorY < -cutY) & (MirrorA > cutA))
    FD_A <- sum((MirrorY > cutY) & (MirrorA < -cutA))
    FD_YA <- sum((MirrorY < -cutY) & (MirrorA < -cutA))
    FDP_AND[i] = (FD_Y + FD_A - c_adj*FD_YA)/max(sum((MirrorY > cutY) & (MirrorA > cutA)),1)
  }
  
  ## OR
  if(length(which(FDP_OR>tgt_FDR))==0){
    thresh_OR = 0
  }else if (all(FDP_OR>tgt_FDR)){
    thresh_OR = Inf
  } else {
    thresh_OR = grid_cut[max(which(FDP_OR > tgt_FDR)) + 1]
  }
  Selected_OR = (1:length(MirrorY))[(MirrorY > thresh_OR)|(MirrorA > thresh_OR)]
  
  ## AND
  if(length(which(FDP_AND>tgt_FDR))==0){
    thresh_AND = 0
  } else if (all(FDP_AND>tgt_FDR)){
    thresh_AND = Inf
  } else {
    thresh_AND = grid_cut[max(which(FDP_AND > tgt_FDR)) + 1]
  }
  Selected_AND = (1:length(MirrorY))[(MirrorY > thresh_AND) & (MirrorA > thresh_AND)]
  
  return(list(Selected_AND=Selected_AND, Selected_OR=Selected_OR, thresh_OR=thresh_OR, thresh_AND=thresh_AND))
}

estimate_ATE <- function(X,A,y,family="gaussian",result_label="NO LABEL"){
  
  df = data.frame(y=y, A=A, X=X)

  model_REG <- glm(y ~ ., data = df, family = family)
  
  yhat1 <- predict(model_REG, newdata = (df %>% mutate(A = 1)), type = "response")
  yhat0 <- predict(model_REG, newdata = (df %>% mutate(A = 0)), type = "response")
  yhat  <- predict(model_REG, newdata = df, type = "response")
  
  model_PS <- glm(A ~ ., data = df %>% dplyr::select(-y),family = binomial)
  PropScore <- predict(model_PS, newdata = df, type = "response")
  Weight1 <- df$A/PropScore
  Weight0 <- (1-df$A)/(1-PropScore)
  
  trim_flag <- (apply(cbind(PropScore, 1-PropScore), MARGIN = 1, min) < 0.05)*1
  Weight1_trim <- df$A/PropScore*(1-trim_flag)
  Weight0_trim <- (1-df$A)/(1-PropScore)*(1-trim_flag)
  
  Weight1_ATO <- df$A*(1-PropScore)
  Weight0_ATO <- (1-df$A)*PropScore
  
  ATE_DIF <- mean(y[A==1]) - mean(y[A==0])
  ATE_REG <- mean(yhat1 - yhat0)
  ATE_IPW <- sum(Weight1*df$y)/sum(Weight1) - sum(Weight0*df$y)/sum(Weight0)
  ATE_trim <- sum(Weight1_trim*df$y)/sum(Weight1_trim) - sum(Weight0_trim*df$y)/sum(Weight0_trim)
  ATO <- sum(Weight1_ATO*df$y)/sum(Weight1_ATO) - sum(Weight0_ATO*df$y)/sum(Weight0_ATO)
  ATE_DR <- ATE_REG + mean((Weight1 - Weight0)*(df$y - yhat))

  return(list(result_label=result_label, 
              ATE_DIF = ATE_DIF,
              ATE_REG = ATE_REG,
              ATE_IPW = ATE_IPW,
              ATE_trim = ATE_trim,
              ATE_DR = ATE_DR,
              ATO = ATO))
}

negaloglike_bin <- function(X,y,coef){
  lp <- as.vector(cbind(1,X) %*% coef)
  return( mean(- y * lp + log(1 + exp(lp))))
}


negaloglike_bin_deriv1 <- function(X,y,coef){
  lp <- as.vector(cbind(1,X) %*% coef)
  nlld1_mat <- diag(sigmoid(lp) - y) %*% X
  ls_ret <- list(
    nlld1_mean = as.vector(colMeans(nlld1_mat)),
    nlld1_cov = cov(nlld1_mat)
  )
  return(ls_ret)
}


gen_corr_matrix <- function(n_obs, mat_dim, mat_cov, dist="gaussian", names_col=NULL){
  mat_pre <- mvrnorm(n_obs, mu = rep(0,mat_dim), Sigma=mat_cov)
  if (dist == "gaussian"){
    mat <- mat_pre
  } else if (dist == "uniform"){
    mat <- (2*pnorm(mat_pre) - 1)*sqrt(3)
  } else if (dist == "binomial"){
    tmp <- pnorm(mat_pre)
    mat <- (tmp > 0.5)*2-1
  }
  
  if(!is.null(names_col)){
    colnames(mat) <- names_col
  }
  return(mat)
}

trans_corr_matrix <- function(mat_pre, dist, names_col=NULL){
  if (dist == "gaussian"){
    mat <- mat_pre
  } else if (dist == "uniform"){
    mat <- (2*pnorm(mat_pre) - 1)*sqrt(3)
  } else if (dist == "binomial"){
    tmp <- pnorm(mat_pre)
    mat <- (tmp > 0.5)*1
  }
  
  if(!is.null(names_col)){
    colnames(mat) <- names_col
  }
  return(mat)
}

eval_suffsets <- function(selected, suffsets, dim_Xm){
  selected <- sort(selected)
  sel_tmp <- unique(((selected-1) %/% dim_Xm) + 1)
  selected_X <- paste0("X",sel_tmp[sel_tmp < 9])
  
  suffset_idx_ = which(sapply(suffsets, function(s) identical(s, selected_X)))
  suffset_idx <- ifelse(length(suffset_idx_)==0, 0, suffset_idx_)
  
  selected_suffset = ifelse(suffset_idx == 0,0,1)
  ls_nsel <- setNames(lapply(1:9, function(m) sum(m == (((selected - 1) %/% dim_Xm) + 1))), paste0("nsel_X",1:9))
  
  return(c(list(suffset_idx=suffset_idx, selected_suffset=selected_suffset), ls_nsel))
}






