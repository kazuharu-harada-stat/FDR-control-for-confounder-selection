require(glmnet)
require(scalreg)


cvfit_lasso <- function(X,y,family,standardize = FALSE,lambda_thresh = "min"){
  tmp <- cv.glmnet(X, y, family = family, standardize = standardize)
  if (lambda_thresh == "min"){
    model_glmnet <- glmnet(X, y, family = family, lambda = tmp$lambda.min, standardize = standardize)
  } else if (lambda_thresh == "1se") {
    model_glmnet <- glmnet(X, y, family = family, lambda = tmp$lambda.1se, standardize = standardize)
  }
  ret <- model_glmnet
  return(ret)
}

obtain_precision_mat <- function(X, beta, family){
  
  if (family == "gaussian"){
    PrecMat <- scalreg(X)
  } else if (family == "binomial"){
    P_beta = as.vector(sigmoid(cbind(1,X) %*% beta))
    W_beta = diag(sqrt(P_beta*(1-P_beta)))
    PrecMat <- scalreg(W_beta%*%X)
  } else {
    stop("family must be gaussian or binomial")
  }
  return(PrecMat$precision)
}

debiased_lasso_for_mirror <- function(X,y,family,standardize = FALSE){
  # Step 2(a) Calc. the Lasso Estimator
  model_glmnet <- cvfit_lasso(X, y, family, standardize)
  coef_pre <- c(model_glmnet$a0, as.numeric(model_glmnet$beta))
  
  # Step 2(b) Estimate Theta following Alg. 5 (Here, we used the scalreg package in R)
  prec_mat <- obtain_precision_mat(X, beta = coef_pre, family = family)
  cov_mat <- cov(X)
  
  # Step 2(c) Calc. the debiased Lasso estimators
  if (family == "gaussian"){
    sigma_hat <- diag(prec_mat%*%cov_mat%*%t(prec_mat))
    coef_unnorm = as.vector(coef_pre[-1] + prec_mat %*% t(X) %*% (y - as.vector(cbind(1, X) %*% coef_pre))/nrow(X))
    coef_norm = coef_unnorm/sigma_hat
  } else if (family == "binomial"){
    nlld1 <- negaloglike_bin_deriv1(X,y,coef_pre)
    sigma_hat <- diag(prec_mat%*%nlld1$nlld1_cov%*%t(prec_mat))
    coef_unnorm = as.vector(coef_pre[-1] - prec_mat %*% nlld1$nlld1_mean)
    coef_norm = coef_unnorm/sigma_hat
  }
  
  ls_ref <- list(
    a0 = model_glmnet$a0,
    beta = coef_unnorm,
    beta_pre = as.numeric(model_glmnet$beta),
    beta_norm = coef_norm,
    prec_mat = prec_mat,
    sigma_hat = sigma_hat
  )
  
  return(ls_ref)
}


getStat_DS <- function(X,Z,y,family,method="glm",...){
  if(method == "glm") {
    ret <- DS_glm(X,Z,y,family)
  } else if (method == "glmnet") {
    ret <- DS_glmnet(X,Z,y,family)
  } else if (method == "cross") {
    ret <- DS_cross(X,Z,y,family)
  } else if (method == "debiased"){
    ret <- DS_debiased(X,Z,y,family)
  }
  return(ret)
}

DS_glm <- function(X,Z,y,family){
  num_row <- nrow(X)
  num_colX <- ncol(X)
  
  if (is.matrix(Z)) {
    num_colZ <- ncol(Z)
    df_tmp <- data.frame(Z=Z, X=X, y=y)
  } else if (is.vector(Z)) {
    num_colZ <- 1
    df_tmp <- data.frame(Z=Z, X=X, y=y)
  } else if (is.null(Z)) {
    num_colZ <- 0
    df_tmp <- data.frame(X=X, y=y)
  } else {stop("Z must be matrix, vector, or NULL")}
  
  # STEP 1
  split_idx <- sample(1:num_row, num_row %/% 2, replace = F)
  
  model.glm1 <- glm(y ~ ., df_tmp[split_idx, ], family = family)
  model.glm2 <- glm(y ~ ., df_tmp[-split_idx, ], family = family)
  
  coef1_se <- summary(model.glm1)$coefficients[, "Std. Error"]
  coef2_se <- summary(model.glm2)$coefficients[, "Std. Error"]
  
  coef_idx <- 1:num_colX + num_colZ
  coef1 <- model.glm1$coefficients[coef_idx + 1]/coef1_se[coef_idx + 1]
  coef2 <- model.glm2$coefficients[coef_idx + 1]/coef2_se[coef_idx + 1]
  
  ret <- list(T1 = coef1, T2 = coef2)
  
  return(ret)
}
family="binomial";
DS_glmnet <- function(X,Z,y,family){
  num_row <- nrow(X)
  num_colX <- ncol(X)
  
  if (is.matrix(Z)) {
    num_colZ <- ncol(Z)
    df_tmp <- data.frame(Z=Z, X=X, y=y)
  } else if (is.vector(Z)) {
    num_colZ <- 1
    df_tmp <- data.frame(Z=Z, X=X, y=y)
  } else if (is.null(Z)) {
    num_colZ <- 0
    df_tmp <- data.frame(X=X, y=y)
  } else {stop("Z must be matrix, vector, or NULL")}
  
  split_idx <- sample(1:num_row, num_row %/% 2, replace = F)
  ZX <- cbind(Z,X)
  ZX1 <- scale(ZX[split_idx,]); ZX2 <- scale(ZX[-split_idx,])
  
  
  
  # Step 2(a) Calc. the Lasso Estimator
  model1 <- cvfit_lasso(ZX1, y[split_idx],family)
  model2 <- cvfit_lasso(ZX2, y[-split_idx],family)
  
  if(family == "gaussian"){
    coef_se1 <- sqrt(c(var(y[split_idx] - predict(model1, newx = ZX1)))/diag(t(ZX1)%*%ZX1))
    coef_se2 <- sqrt(c(var(y[-split_idx] - predict(model2, newx = ZX2)))/diag(t(ZX2)%*%ZX2))
  } else if (family == "binomial"){
    prob1 <- predict(model1, newx = ZX1, type="response")
    prob2 <- predict(model2, newx = ZX2, type="response")
    
    coef_se1 <- sqrt(1/diag(t(ZX1) %*% diag(as.vector((y[split_idx] - prob1)^2)) %*% ZX1))
    coef_se2 <- sqrt(1/diag(t(ZX2) %*% diag(as.vector((y[-split_idx] - prob2)^2)) %*% ZX2))
  }
  
  coef_idx <- 1:num_colX + num_colZ
  coef1 <- as.numeric(model1$beta)[coef_idx]/coef_se1[coef_idx]
  coef2 <- as.numeric(model2$beta)[coef_idx]/coef_se2[coef_idx]
  
  ret <- list(T1 = coef1, T2 = coef2)
  
  return(ret)
}

DS_cross <- function(X,Z,y,family){
  num_row <- nrow(X)
  num_colX <- ncol(X)
  
  if (is.matrix(Z)) {
    num_colZ <- ncol(Z)
    df_tmp <- data.frame(Z=Z, X=X, y=y)
  } else if (is.vector(Z)) {
    num_colZ <- 1
    df_tmp <- data.frame(Z=Z, X=X, y=y)
  } else if (is.null(Z)) {
    num_colZ <- 0
    df_tmp <- data.frame(X=X, y=y)
  } else {stop("Z must be matrix, vector, or NULL")}
  
  seq_colZ <- seq_len(num_colZ)
  split_idx <- sample(1:num_row, num_row %/% 2, replace = F)
  ZX <- cbind(Z,X)
  ZX1 <- scale(ZX[split_idx,]); ZX2 <- scale(ZX[-split_idx,])
  
  
  cvlasso_type = ifelse(num_colX >= num_row/2, "1se", "min")
  
  # Step 2(a) Calc. the Lasso Estimator
  model1 <- cvfit_lasso(ZX1, y[split_idx],family,lambda_thresh = cvlasso_type)
  model2 <- cvfit_lasso(ZX2, y[-split_idx],family,lambda_thresh = cvlasso_type)
  
  sel_model1 <- sort(union(which(model1$beta != 0), seq_colZ))
  sel_model2 <- sort(union(which(model2$beta != 0), seq_colZ))
  
  length(sel_model1)
  length(sel_model2)
  
  if (length(sel_model1) == 0){
    coef1 <- rep(0,num_colZ+num_colX)
  } else {
    model12 <- glm(y[-split_idx] ~ ZX2[,sel_model1], family=family)
    coef12_se <- summary(model12)$coefficients[, "Std. Error"]
    coef12_std <- model12$coefficients/coef12_se
    coef1 <- rep(0,num_colZ+num_colX)
    coef1[sel_model1] <- coef12_std[-c(1)]
  }
  
  if (length(sel_model2) == 0){
    coef2 <- rep(0,num_colZ+num_colX)
  } else {
    model21 <- glm(y[split_idx] ~ ZX1[,sel_model2], family=family)
    coef21_se <- summary(model21)$coefficients[, "Std. Error"]
    coef21_std <- model21$coefficients/coef21_se
    coef2 <- rep(0,num_colZ+num_colX)
    coef2[sel_model2] <- coef21_std[-c(1)]
  }
  
  if(num_colZ > 0){
    coef1 <- coef1[-seq_colZ]
    coef2 <- coef2[-seq_colZ]
  }
  
  ret <- list(T1 = coef1, T2 = coef2)
  
  return(ret)
}

DS_debiased <- function(X,Z,y,family){
  num_row <- nrow(X)
  num_colX <- ncol(X)
  
  if (is.matrix(Z)) {
    num_colZ <- ncol(Z)
    df_tmp <- data.frame(Z=Z, X=X, y=y)
  } else if (is.vector(Z)) {
    num_colZ <- 1
    df_tmp <- data.frame(Z=Z, X=X, y=y)
  } else if (is.null(Z)) {
    num_colZ <- 0
    df_tmp <- data.frame(X=X, y=y)
  } else {stop("Z must be matrix, vector, or NULL")}
  
  # Step 1 Split the data set
  split_idx <- sample(1:num_row, num_row %/% 2, replace = F)
  ZX <- cbind(Z,X)
  ZX1 <- scale(ZX[split_idx,]); ZX2 <- scale(ZX[-split_idx,])
  y1 <- y[split_idx]; y2 <- y[-split_idx]
  
  # Step 2 exec. debiased Lasso
  dlasso1 <- debiased_lasso_for_mirror(ZX1,y1,family=family,standardize = FALSE)
  dlasso2 <- debiased_lasso_for_mirror(ZX2,y2,family=family,standardize = FALSE)
  
  # Step 3 Select the features
  coef_idx <- 1:num_colX + num_colZ
  coef1 <- dlasso1$beta_norm[coef_idx]
  coef2 <- dlasso2$beta_norm[coef_idx]
  
  ret <- list(T1 = coef1, T2 = coef2)
  
  return(ret)
}

getConf_Paired <- function(TY1, TY2, TA1, TA2, tgt_FDR = 0.1, mirror_type = "sum", n_grid=1000){
  
  # calc Mirror statistics
  MirrorY <- Mirror(TY1, TY2, mirror_type=mirror_type, tgt_FDR=tgt_FDR, n_grid=n_grid)$Statistic
  MirrorA <- Mirror(TA1, TA2, mirror_type=mirror_type, tgt_FDR=tgt_FDR, n_grid=n_grid)$Statistic
  
  # calc FDP
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
    FDP_AND[i] = (FD_Y + FD_A - FD_YA)/max(sum((MirrorY > cutY) & (MirrorA > cutA)),1)
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
  
  ls_ret <- list(Selected_AND=Selected_AND,
                 Selected_OR=Selected_OR,
                 Threshold_OR=thresh_OR,
                 Threshold_AND=thresh_AND)
  return(ls_ret)
}



getConf_Unified <- function(TY1, TY2, TA1, TA2, tgt_FDR = 0.1, n_grid=1000,
                            set_type = "OR", mirror_type = "sum",
                            intensity = "rad", concordance = "cos", orientation = "sin"){
  
  th1 <- atan2(TY1, TA1)
  th2 <- atan2(TY2, TA2)
  
  if(intensity == "rad"){
    r1 = apply(cbind(TA1, TY1), MARGIN = 1, function(u) sqrt(sum(u^2)))
    r2 = apply(cbind(TA2, TY2), MARGIN = 1, function(u) sqrt(sum(u^2)))
  } else if (intensity == "max"){
    r1 = apply(abs(cbind(TA1, TY1)), MARGIN = 1, max)
    r2 = apply(abs(cbind(TA2, TY2)), MARGIN = 1, max)
  } else if (intensity == "min"){
    r1 = apply(abs(cbind(TA1, TY1)), MARGIN = 1, min)
    r2 = apply(abs(cbind(TA2, TY2)), MARGIN = 1, min)
  }
  
  if(concordance == "cos"){
    gap_angle <- cos(th1 - th2)
  } else if (concordance == "sgn") {
    gap_angle <- sign(cos(th1 - th2))
  }
  
  if(mirror_type == "sum"){
    M_OR <- (r1 + r2)*gap_angle
  } else if (mirror_type == "product"){
    M_OR <- r1 * r2 * gap_angle
  }
  
  if(set_type == "OR"){
    M_vec <- M_OR
  } else if (set_type == "AND"){
    
    if(orientation == "sin"){
      abs_angle <- sin(2*th1)*sin(2*th2)
    } else if (orientation == "sgn") {
      abs_angle <- sign(sin(2*th1)*sin(2*th2))
    }
    
    M_vec <- M_OR * abs_angle
  }
  
  # thresholding
  ls_t <- seq(0,max(abs(M_vec)),max(abs(M_vec))/n_grid)
  tauq <- calc_threshold(M_vec, ls_t, tgt_FDR)
  S1 <- which(M_vec > tauq)
  
  names(M_vec) <- NULL
  names(S1) <- NULL

  ret <- list(
    Selected = S1, Statistic = M_vec, Threshold = tauq
  )
  
  return(ret)
}


SelConf_DS_compare_mirror <- function(X,A,y, family_Y="gaussian", method="glm", tgt_FDR = 0.1, n_grid=1000){
  fitY <- getStat_DS(X=X,Z=A,y=y,family=family_Y,method=method)
  fitA <- getStat_DS(X=X,Z=NULL,y=A,family="binomial",method=method)
  
  fit_Modelwise <- getConf_Modelwise(TY1=fitY$T1, TY2=fitY$T2, TA1=fitA$T1, TA2=fitA$T2,
                                     tgt_FDR = tgt_FDR, mirror_type = "sum", n_grid=n_grid)
  
  fit_Paired_sum <- getConf_Paired(TY1=fitY$T1, TY2=fitY$T2, TA1=fitA$T1, TA2=fitA$T2,
                                   tgt_FDR = tgt_FDR, mirror_type = "sum", n_grid=n_grid)
  
  fit_Paired_product <- getConf_Paired(TY1=fitY$T1, TY2=fitY$T2, TA1=fitA$T1, TA2=fitA$T2,
                                       tgt_FDR = tgt_FDR, mirror_type = "product", n_grid=n_grid)
  
  
  ls_ret <- list(
    Y_Marginal_NA_NA_NA_NA = fit_Modelwise$Selected_Y,
    A_Marginal_NA_NA_NA_NA = fit_Modelwise$Selected_A,
    
    OR_Modelwise_NA_NA_NA_NA = fit_Modelwise$Selected_OR,
    AND_Modelwise_NA_NA_NA_NA = fit_Modelwise$Selected_AND,
    
    OR_Paired_sum_NA_NA_NA = fit_Paired_sum$Selected_OR,
    AND_Paired_sum_NA_NA_NA = fit_Paired_sum$Selected_AND,
    
    OR_Paired_product_NA_NA_NA = fit_Paired_product$Selected_OR,
    AND_Paired_product_NA_NA_NA = fit_Paired_product$Selected_AND
  )
  
  # All petterns of Unified Mirror
  ## OR
  for (mirror_type in c("sum","product")) {
    for (intensity in c("rad", "max")) {
      for(concordance in c("cos", "sgn")){
        
        tmp <- getConf_Unified(TY1=fitY$T1, TY2=fitY$T2, TA1=fitA$T1, TA2=fitA$T2,
                               tgt_FDR = tgt_FDR, n_grid=n_grid, set_type = "OR",
                               mirror_type = mirror_type, intensity = intensity, concordance = concordance, orientation = NULL)
        
        item <- sprintf("OR_Unified_%s_%s_%s_NA", mirror_type, intensity, concordance)
        ls_ret[[item]] <- tmp$Selected
      }
    }
  }
  
  ## AND
  for (mirror_type in c("sum","product")) {
    for (intensity in c("rad", "min")) {
      for(concordance in c("cos", "sgn")){
        for (orientation in c("sin", "sgn")) {
          tmp <- getConf_Unified(TY1=fitY$T1, TY2=fitY$T2, TA1=fitA$T1, TA2=fitA$T2,
                                 tgt_FDR = tgt_FDR, n_grid=n_grid, set_type = "AND",
                                 mirror_type = mirror_type, intensity = intensity, concordance = concordance, orientation = orientation)
          
          item <- sprintf("AND_Unified_%s_%s_%s_%s", mirror_type, intensity, concordance, orientation)
          ls_ret[[item]] <- tmp$Selected
        }
      }
    }
  }
  
  return(ls_ret)
}



getConf_Modelwise <- function(TY1, TY2, TA1, TA2, tgt_FDR = 0.1, mirror_type = "sum", n_grid=1000){
  
  # calc Mirror statistics
  MirrorY <- Mirror(TY1, TY2, mirror_type=mirror_type, tgt_FDR=tgt_FDR, n_grid=n_grid)$Selected
  MirrorA <- Mirror(TA1, TA2, mirror_type=mirror_type, tgt_FDR=tgt_FDR, n_grid=n_grid)$Selected
  
  names(MirrorY) <- NULL; names(MirrorA) <- NULL
  
  ls_ret <- list(
    Selected_AND=intersect(MirrorY, MirrorA),
    Selected_OR=union(MirrorY, MirrorA),
    Selected_Y=MirrorY,
    Selected_A=MirrorA
  )
  return(ls_ret)
}



SelConf_MDS_compare_mirror <- function(X,A,y,family_Y="gaussian", method="glm", tgt_FDR = 0.1, n_repeat=100, n_grid=1000){
  num_colX <- ncol(X)

  df_DS_all <- foreach(itr = 1:n_repeat, .combine=rbind)%do%{
    ls_DS <- SelConf_DS_compare_mirror(X,A,y,family_Y,method, tgt_FDR,n_grid)
    
    df_DS_inc <- foreach(m = names(ls_DS), .combine=rbind)%do%{

      df_tmp <- data.frame(
        method = m,
        idx_X = 1:num_colX,
        included = (1:num_colX %in% ls_DS[[m]])*1,
        value_I = (1:num_colX %in% ls_DS[[m]])*1 / max(length(ls_DS[[m]]),1)
      )
      return(df_tmp)
    }
    df_DS_inc$idx_MDS <- itr
    return(df_DS_inc)
  }

  df_MDS <- df_DS_all %>%
    group_by(method, idx_X) %>%
    summarise(
      cnt_row = n(),
      mean_I = mean(value_I)
    ) %>%
    group_by(method) %>%
    arrange(method, mean_I) %>%
    mutate(
      cumsum_I = cumsum(mean_I),
      flag_FDP = (cumsum_I <= tgt_FDR)
    ) %>% filter(flag_FDP == FALSE)
  
  # write.csv(df_MDS, file="test_03.csv")
  # write.csv(df_DS_all, file="test_DS_03.csv")

  ls_ret <- split(df_MDS$idx_X, df_MDS$method)

  items <- unique(df_DS_all$method)

  for (item in items) {
    if(!(item %in% names(ls_ret))) {
      ls_ret[[item]] <- numeric(0)
    }
  }

  return(ls_ret)
}




