require(mvtnorm)
# require(hdi)



getStat_Pval <- function(X,Z,y,family,method="glm"){
  if(method == "glm") {
    ret <- est_glm(X,Z,y,family)
  } else if (method == "cross"){
    ret <- est_cross(X,Z,y,family)
  } else if (method == "debiased"){
    ret <- est_debiased(X,Z,y,family)
  } else if (method == "marginal"){
    ret <- est_marginal(X,Z,y,family)
  }
  return(ret)
}

est_glm <- function(X,Z,y,family){
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
  
  model_glm <- glm(y ~ ., df_tmp, family = family)
  
  coef_idx = 1:num_colX + num_colZ + 1
  
  ls_ret <- list(
    coef = summary(model_glm)$coefficients[coef_idx,1],
    SE = summary(model_glm)$coefficients[coef_idx,2],
    z_values = summary(model_glm)$coefficients[coef_idx,3],
    p_values = summary(model_glm)$coefficients[coef_idx,4]
  )
  
  return(ls_ret)
}


est_debiased <- function(X,Z,y,family){
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
  ZX <- cbind(Z,X)
  ZX_std <- scale(ZX)

  # Step 2 exec. debiased Lasso
  dlasso <- debiased_lasso_for_mirror(ZX_std,y,family=family,standardize = FALSE)
  
  # Step 3 Select the features
  coef_idx <- 1:num_colX + num_colZ
  
  ret <- list(
    coef = dlasso$beta[coef_idx],
    SE = dlasso$sigma_hat[coef_idx],
    z_values = dlasso$beta_norm[coef_idx],
    p_values = pnorm(-abs(dlasso$beta_norm[coef_idx]))*2
  )
  
  return(ret)
}


est_cross <- function(X,Z,y,family){
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
  model_lasso <- cvfit_lasso(ZX1, y[split_idx],family,lambda_thresh = cvlasso_type)
  selected <- sort(union(which(model_lasso$beta != 0), seq_colZ))
  
  # Step 2(b) Estimate parameters
  coef <- rep(0,num_colX + num_colZ)
  SE <- rep(999,num_colX + num_colZ)
  z_values <- rep(0,num_colX + num_colZ)
  p_values <- rep(1,num_colX + num_colZ)
  
  if (length(selected) > 0){
    model_glm <- glm(y[-split_idx] ~ ZX2[,selected], family=family)
    
    coef[selected] <- model_glm$coefficients[-c(1)]
    SE[selected] <- summary(model_glm)$coefficients[-c(1), 2]
    z_values[selected] <- summary(model_glm)$coefficients[-c(1), 3]
    p_values[selected] <- summary(model_glm)$coefficients[-c(1), 4]
  }
  
  if(num_colZ > 0){
    coef <- coef[-seq_colZ]
    SE <- SE[-seq_colZ]
    z_values <- z_values[-seq_colZ]
    p_values <- p_values[-seq_colZ]
  }
  
  ls_ret <- list(
    coef = coef,
    SE = SE,
    z_values = z_values,
    p_values = p_values
  )
  
  return(ls_ret)
}


est_marginal <- function(X,Z,y,family){
  num_row <- nrow(X)
  num_colX <- ncol(X)
  
  if (is.matrix(Z)) {
    num_colZ <- ncol(Z)
  } else if (is.vector(Z)) {
    num_colZ <- 1
  } else if (is.null(Z)) {
    num_colZ <- 0
  } else {stop("Z must be matrix, vector, or NULL")}
  
  coef <- rep(0,num_colX)
  SE <- rep(999,num_colX)
  z_values <- rep(0,num_colX)
  p_values <- rep(1,num_colX)

  
  for (idx in 1:num_colX) {
    if( num_colZ == 0 ) {
      model_glm <- glm(y ~ X[,idx], family=family)
    } else {
      model_glm <- glm(y ~ X[,idx] + Z, family=family)
    }
    
    coef[idx] <- model_glm$coefficients[2]
    SE[idx] <- summary(model_glm)$coefficients[2,2]
    z_values[idx] <- summary(model_glm)$coefficients[2,3]
    p_values[idx] <- summary(model_glm)$coefficients[2,4]
  }
  
  ls_ret <- list(
    coef = coef,
    SE = SE,
    z_values = z_values,
    p_values = p_values
  )
  
  return(ls_ret)
}


SelConf_Qval <- function(X,A,y,family_Y="gaussian", method="glm", tgt_FDR = 0.1){
  
  fitY <- getStat_Pval(X,A,y,family=family_Y,method=method)
  fitA <- getStat_Pval(X,NULL,A,family=family_Y,method=method)
  
  Selected_Y_pval <- which(fitY$p_values < 0.05)
  Selected_A_pval <- which(fitA$p_values < 0.05)

  Selected_Y_BHq <- which(p.adjust(fitY$p_values, method="BH") < tgt_FDR)
  Selected_A_BHq <- which(p.adjust(fitA$p_values, method="BH") < tgt_FDR)
  
  Selected_Y_BYq <- which(p.adjust(fitY$p_values, method="BY") < tgt_FDR)
  Selected_A_BYq <- which(p.adjust(fitA$p_values, method="BY") < tgt_FDR)
  
  z_values_OR <- abs(cbind(fitY$z_values, fitA$z_values))
  p_values_OR <- apply(z_values_OR, 1, function(u) pmvnorm(lower = u, mean = c(0,0), sigma = diag(2)))*4
  Selected_OR_pval <- which(p_values_OR < 0.05)
  Selected_OR_BHq <- which(p.adjust(p_values_OR, method="BH") < tgt_FDR)
  Selected_OR_BYq <- which(p.adjust(p_values_OR, method="BY") < tgt_FDR)
  
  p_values_OR_min <- pmin(fitY$p_values, fitA$p_values)
  Selected_OR_pval_min <- which(p_values_OR_min < 0.05)
  Selected_OR_BHq_min <- which(p.adjust(p_values_OR_min, method="BH") < tgt_FDR)
  Selected_OR_BYq_min <- which(p.adjust(p_values_OR_min, method="BY") < tgt_FDR)

  p_values_AND <- pmax(fitY$p_values, fitA$p_values)
  Selected_AND_pval <- which(p_values_AND < 0.05)
  Selected_AND_BHq <- which(p.adjust(p_values_AND, method="BH") < tgt_FDR)
  Selected_AND_BYq <- which(p.adjust(p_values_AND, method="BY") < tgt_FDR)
  
  ls_ret <- list(
    Y_Marginal_pval = Selected_Y_pval,
    A_Marginal_pval = Selected_A_pval,
    
    OR_Modelwise_pval = union(Selected_Y_pval, Selected_A_pval),
    AND_Modelwise_pval = intersect(Selected_Y_pval, Selected_A_pval),
    
    OR_Joint_pval = Selected_OR_pval,
    AND_Joint_pval = Selected_AND_pval,
    
    OR_min_pval = Selected_OR_pval_min,

    Y_Marginal_BH = Selected_Y_BHq,
    A_Marginal_BH = Selected_A_BHq,
    
    OR_Modelwise_BH = union(Selected_Y_BHq, Selected_A_BHq),
    AND_Modelwise_BH = intersect(Selected_Y_BHq, Selected_A_BHq),
    
    OR_Joint_BH = Selected_OR_BHq,
    AND_Joint_BH = Selected_AND_BHq,
    
    OR_min_BH = Selected_OR_BHq_min,
    
    Y_Marginal_BY = Selected_Y_BYq,
    A_Marginal_BY = Selected_A_BYq,
    
    OR_Modelwise_BY = union(Selected_Y_BYq, Selected_A_BYq),
    AND_Modelwise_BY = intersect(Selected_Y_BYq, Selected_A_BYq),
    
    OR_Joint_BY = Selected_OR_BYq,
    AND_Joint_BY = Selected_AND_BYq,
    
    OR_min_BY = Selected_OR_BYq_min
  )
  
  return(ls_ret)
}




