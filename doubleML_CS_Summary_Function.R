

#------------------------------------------------------------------------------
#                   Double ML Cross-Sectional DiD 
#-----------------------------------------------------------------------------

doubleML_did_cs <- function(data, y_var, d_var, t_var, x_vars, methods, cf = 5, 
                            feature_eng = FALSE, bin_cut = 0.02, corr_cut = 0.9, cluster_var = NULL) {
  
  Y  <- data[[y_var]]
  D  <- as.numeric(data[[d_var]])
  Tt <- as.numeric(data[[t_var]])
  X_0  <- as.matrix(data[, x_vars])
  n  <- length(Y)
  if(!is.null(cluster_var)){
    C <- data[[cluster_var]]
  }
  
  if (feature_eng) {
    X_expanded <- design_matrix(X_0, int = "all", int_d = 2, 
                                poly = "all", poly_d = 2)
    X <- data_screen(X_expanded, treat = D, 
                     bin_cut = bin_cut, corr_cut = corr_cut, 
                     quiet = TRUE)
  } else {
    X <- X_0
  }
  
  # Kombiniertes Treatment: (D,T) 
  DT_combined <- 2*D + Tt  
  
  # ---- 1) m(X)
  np_m <- nuisance_parameters(
    NuPa = "D.hat", X = X, D = D,
    methods = methods, cf = cf,
    stacking = "short", stratify = TRUE, ensemble_type = "nnls",
    quiet = TRUE
  )
  m_hat <- as.numeric(np_m$nuisance_parameters$D.hat)
  
  # Kombiniertes Treatment: (D,T) 
  DT_combined <- 2*D + Tt  
  # ---- 2) g(d,t,X) 
  np_g <- nuisance_parameters(
    NuPa = "Y.hat.d",      
    X = X,
    Y = Y,
    D = DT_combined,       
    methods = methods,
    cf = cf,
    stacking = "short",
    stratify = TRUE,
    ensemble_type = "nnls",
    quiet = TRUE
  )
  
  #
  g00_hat <- np_g$nuisance_parameters$Y.hat.d[, 1]  # D=0,T=0
  g01_hat <- np_g$nuisance_parameters$Y.hat.d[, 2]  # D=0,T=1
  g10_hat <- np_g$nuisance_parameters$Y.hat.d[, 3]  # D=1,T=0
  g11_hat <- np_g$nuisance_parameters$Y.hat.d[, 4]  # D=1,T=1
  
  # ---- 3) Score ----
  
  psi_a <- -D / mean(D)
  
  psi_b <- D / mean(D) * (g11_hat - g10_hat - (g01_hat - g00_hat)) + 
    (Tt * D / mean(D * Tt)) * (Y - g11_hat) -
    (D * (1 - Tt) / mean(D * (1 - Tt))) * (Y - g10_hat) - 
    (m_hat * (1 - D) * Tt / (1 - m_hat)) / mean(m_hat * (1 - D) * Tt / (1 - m_hat)) * (Y - g01_hat) + 
    (m_hat * (1 - D) * (1 - Tt) / (1 - m_hat)) / mean(m_hat * (1 - D) * (1 - Tt) / (1 - m_hat)) * (Y - g00_hat)
  
  theta_hat <- -mean(psi_b) / mean(psi_a)
  
  psi <- psi_a * theta_hat + psi_b
  influence_fun <- -1/mean(psi_a) * psi
  if(is.null(cluster_var)){
    se_hat <- sqrt(mean(psi^2) / (mean(psi_a)^2 * n))
    clustered = FALSE
  } else{
    clusters <- unique(C)
    cluster_sums_squared <- rep(NA, length(clusters))
    index = 1
    for(i in clusters){
      cluster_sums_squared[index] <- sum(influence_fun[C == i])^2
      index = index + 1
    }
    var_hat <- 1/n * sum(cluster_sums_squared) 
    se_hat <- sqrt(var_hat/n) 
    clustered = TRUE
  }
  # ---- 4) Return ----
  out <- list(
    coefficients = c(ATT = theta_hat),
    se = se_hat,
    clustered = clustered,
    t = theta_hat / se_hat,
    p = 2 * (1 - pnorm(abs(theta_hat / se_hat))),
    ci_lower = theta_hat - 1.96 * se_hat,
    ci_upper = theta_hat + 1.96 * se_hat,
    N = n,
    N_treated = sum(D),
    N_control = sum(1 - D),
    N_post = sum(Tt),
    N_pre = sum(1 - Tt),
    y_var = y_var,
    d_var = d_var,
    t_var = t_var,
    x_vars = x_vars,
    n_covariates = ncol(X),
    cf = cf,
    methods = names(methods),
    feature_eng = feature_eng
  )
  
  class(out) <- "doubleml_did_cs"
  return(out)
}



#-----------------------------------------------------------
#                    summary function                     
#----------------------------------------------------------


summary.doubleml_did_cs <- function(object, digits = 4, ...) {
  
  results <- data.frame(
    Estimate = round(object$coefficients["ATT"], digits),
    `Std. Error` = round(object$se, digits),
    `t value` = round(object$t, digits),
    `Pr(>|t|)` = format.pval(object$p, digits = digits),
    `95% CI Lower` = round(object$ci_lower, digits),
    `95% CI Upper` = round(object$ci_upper, digits),
    check.names = FALSE
  )
  
  rownames(results) <- "ATT"
  
  print(results, quote = FALSE)
  
  invisible(object)
}

print.doubleml_did_cs <- function(x, ...) {
  cat("\nDoubleML DID-CS Estimate:\n")
  cat("  ATT:", round(x$coefficients["ATT"], 4), 
      "(SE:", round(x$se, 4), ")\n")
  cat("  p-value:", format.pval(x$p, digits = 4), "\n")
  cat("  N:", x$N, "\n\n")
  invisible(x)
}
