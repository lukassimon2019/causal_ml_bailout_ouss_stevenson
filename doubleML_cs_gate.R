doubleML_did_cs_gate <- function(data, y_var, d_var, t_var, x_vars, methods, cf = 5,
                                 feature_eng = FALSE, bin_cut = 0.02, corr_cut = 0.9,
                                 cluster_var = NULL,
                                 heterogeneity = FALSE, gate_var = NULL) {
  
  Y  <- data[[y_var]]
  D  <- as.numeric(data[[d_var]])
  Tt <- as.numeric(data[[t_var]])
  X_0 <- as.matrix(data[, x_vars])
  n  <- length(Y)
  
  if (!is.null(cluster_var)) {
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
  DT_combined <- 2 * D + Tt
  
  # ---- 1) m(X)
  np_m <- nuisance_parameters(
    NuPa = "D.hat", X = X, D = D,
    methods = methods, cf = cf,
    stacking = "short", stratify = TRUE, ensemble_type = "nnls",
    quiet = TRUE
  )
  m_hat <- as.numeric(np_m$nuisance_parameters$D.hat)
  
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
  
  g00_hat <- np_g$nuisance_parameters$Y.hat.d[, 1]  # D=0,T=0
  g01_hat <- np_g$nuisance_parameters$Y.hat.d[, 2]  # D=0,T=1
  g10_hat <- np_g$nuisance_parameters$Y.hat.d[, 3]  # D=1,T=0
  g11_hat <- np_g$nuisance_parameters$Y.hat.d[, 4]  # D=1,T=1
  
  # ---- 3) Score (ATT)
  psi_a <- -D / mean(D)
  
  psi_b <- D / mean(D) * (g11_hat - g10_hat - (g01_hat - g00_hat)) +
    (Tt * D / mean(D * Tt)) * (Y - g11_hat) -
    (D * (1 - Tt) / mean(D * (1 - Tt))) * (Y - g10_hat) -
    (m_hat * (1 - D) * Tt / (1 - m_hat)) /
    mean(m_hat * (1 - D) * Tt / (1 - m_hat)) * (Y - g01_hat) +
    (m_hat * (1 - D) * (1 - Tt) / (1 - m_hat)) /
    mean(m_hat * (1 - D) * (1 - Tt) / (1 - m_hat)) * (Y - g00_hat)
  
  theta_hat <- -mean(psi_b) / mean(psi_a)
  
  psi <- psi_a * theta_hat + psi_b
  influence_fun <- -1 / mean(psi_a) * psi
  
  if (is.null(cluster_var)) {
    se_hat <- sqrt(mean(psi^2) / (mean(psi_a)^2 * n))
    clustered <- FALSE
  } else {
    clusters <- unique(C)
    cluster_sums_squared <- rep(NA_real_, length(clusters))
    index <- 1
    for (cl in clusters) {
      cluster_sums_squared[index] <- sum(influence_fun[C == cl])^2
      index <- index + 1
    }
    var_hat <- 1 / n * sum(cluster_sums_squared)
    se_hat <- sqrt(var_hat / n)
    clustered <- TRUE
  }
  
  # ---- 4) Return object (base)
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
    feature_eng = feature_eng,
    cluster_var = cluster_var
  )
  
  # ------ 5) GATE  -------
  if (heterogeneity) {
    
    if (is.null(gate_var)) stop("heterogeneity = TRUE requires gate_var")
    
    G <- as.factor(data[[gate_var]])
    
    # ATT -> restrict to treated units
    idx <- (D == 1)
    n1 <- sum(idx)
    
    if (n1 == 0) stop("No treated units found (D==1). Cannot compute GATE.")
    
    psi_a_g <- psi_a[idx]
    psi_b_g <- psi_b[idx]
    G_g <- droplevels(G[idx])
    
    # pseudo-outcome 
    y_tilde <- -psi_b_g / psi_a_g
    
    # regression for group means (no intercept -> one coef per group)
    gate_fit <- lm(y_tilde ~ 0 + G_g)
    tau_hat <- coef(gate_fit)
    
    # ensure names align with group labels
    groups <- levels(G_g)
    tau_hat <- tau_hat[paste0("G_g", groups)]
    names(tau_hat) <- groups
    
    # Influence function per group (defined on treated subsample)
    IF_mat <- matrix(0, nrow = n1, ncol = length(groups))
    colnames(IF_mat) <- groups
    
    for (g in groups) {
      I_g <- as.numeric(G_g == g)
      p_g_hat <- mean(I_g)
      IF_mat[, g] <- I_g / p_g_hat * (y_tilde - tau_hat[g])
    }
    
    # Clustered SE using Logic from above 
    if (is.null(cluster_var)) {
      gate_se <- apply(IF_mat, 2, function(IF_g) sqrt(mean(IF_g^2) / n1))
      gate_clustered <- FALSE
    } else {
      C_g <- C[idx]
      clusters_g <- unique(C_g)
      gate_se <- rep(NA_real_, length(groups))
      names(gate_se) <- groups
      
      for (g in groups) {
        IF_g <- IF_mat[, g]
        cluster_sums_squared <- rep(NA_real_, length(clusters_g))
        j <- 1
        for (cl in clusters_g) {
          cluster_sums_squared[j] <- sum(IF_g[C_g == cl])^2
          j <- j + 1
        }
        var_hat_g <- 1 / n1 * sum(cluster_sums_squared)
        gate_se[g] <- sqrt(var_hat_g / n1)
      }
      gate_clustered <- TRUE
    }
    
    gate_t <- tau_hat / gate_se
    gate_p <- 2 * (1 - pnorm(abs(gate_t)))
    gate_ci_lower <- tau_hat - 1.96 * gate_se
    gate_ci_upper <- tau_hat + 1.96 * gate_se
    
    out$gate <- list(
      tau = tau_hat,
      se = gate_se,
      clustered = gate_clustered,
      t = gate_t,
      p = gate_p,
      ci_lower = gate_ci_lower,
      ci_upper = gate_ci_upper,
      groups = groups,
      N_treated = n1,
      gate_var = gate_var,
      influence_function = IF_mat  # (treated-subsample IFs; columns = groups)
    )
  }
  
  class(out) <- "doubleml_did_cs_gate"
  return(out)
}