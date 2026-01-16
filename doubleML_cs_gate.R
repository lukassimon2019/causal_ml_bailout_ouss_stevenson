
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
    groups <- levels(G)
    
    # containers
    tau_hat <- rep(NA_real_, length(groups))
    names(tau_hat) <- groups
    
    IF_mat <- matrix(NA_real_, nrow = n, ncol = length(groups))
    colnames(IF_mat) <- groups
    
    for (j in seq_along(groups)) {
      g <- groups[j]
      I_g <- as.numeric(G == g)
      
      denom <- mean(I_g * psi_a)
      if (abs(denom) < 1e-12) {
        stop(paste0("Group ", g, " has near-zero denom E[I_g*psi_a]."))
      }
      
      num <- mean(I_g * psi_b)
      tau_hat[g] <- - num / denom
      
      # Influence function for tau_g
      psi_g <- I_g * (psi_a * tau_hat[g] + psi_b)
      IF_mat[, j] <- -1 / denom * psi_g
    }
    
    # SEs using the same logic as ATT, applied per group
    gate_se <- rep(NA_real_, length(groups))
    names(gate_se) <- groups
    
    if (is.null(cluster_var)) {
      for (j in seq_along(groups)) {
        gate_se[j] <- sqrt(mean(IF_mat[, j]^2) / n)
      }
      gate_clustered <- FALSE
    } else {
      clusters <- unique(C)
      for (j in seq_along(groups)) {
        infl <- IF_mat[, j]
        cluster_sums_squared <- rep(NA_real_, length(clusters))
        k <- 1
        for (cl in clusters) {
          cluster_sums_squared[k] <- sum(infl[C == cl])^2
          k <- k + 1
        }
        var_hat_g <- 1 / n * sum(cluster_sums_squared)
        gate_se[j] <- sqrt(var_hat_g / n)
      }
      gate_clustered <- TRUE
    }
    
    gate_t <- tau_hat / gate_se
    gate_p <- 2 * (1 - pnorm(abs(gate_t)))
    
    out$gate <- list(
      tau = tau_hat,
      se = gate_se,
      clustered = gate_clustered,
      t = gate_t,
      p = gate_p,
      ci_lower = tau_hat - 1.96 * gate_se,
      ci_upper = tau_hat + 1.96 * gate_se,
      groups = groups,
      N = n,
      gate_var = gate_var,
      influence_function = IF_mat
    )
  }
  class(out) <- "doubleML_did_cs_gate"
  return(out)
}





#-----------------------------------------------------------
#   summary + print methods for DID-CS with optional GATE
#   class: "doubleML_did_cs_gate"
#-----------------------------------------------------------

summary.doubleML_did_cs_gate <- function(object, digits = 4, ...) {
  
  # --- ATT table
  att_tab <- data.frame(
    Estimate       = round(object$coefficients["ATT"], digits),
    `Std. Error`   = round(object$se, digits),
    `t value`      = round(object$t, digits),
    `Pr(>|t|)`     = format.pval(object$p, digits = digits),
    `95% CI Lower` = round(object$ci_lower, digits),
    `95% CI Upper` = round(object$ci_upper, digits),
    check.names = FALSE
  )
  rownames(att_tab) <- "ATT"
  
  cat("\n===============================================\n")
  cat("  DoubleML DID Callaway-Sant'Anna Estimation\n")
  cat("===============================================\n")
  
  cat("\nSample Information:\n")
  cat(sprintf("  Total N: %d\n", object$N))
  cat(sprintf("  Treated: %d (%.1f%%)\n", object$N_treated, 
              100 * object$N_treated / object$N))
  cat(sprintf("  Control: %d (%.1f%%)\n", object$N_control, 
              100 * object$N_control / object$N))
  cat(sprintf("  Post-period: %d (%.1f%%)\n", object$N_post, 
              100 * object$N_post / object$N))
  
  cat("\nModel Specification:\n")
  cat(sprintf("  Outcome: %s\n", object$y_var))
  cat(sprintf("  Treatment: %s\n", object$d_var))
  cat(sprintf("  Time: %s\n", object$t_var))
  cat(sprintf("  Covariates: %d\n", object$n_covariates))
  cat(sprintf("  Cross-fitting folds: %d\n", object$cf))
  cat(sprintf("  ML methods: %s\n", paste(object$methods, collapse = ", ")))
  if (!is.null(object$cluster_var)) {
    cat(sprintf("  Clustered SEs: %s\n", object$cluster_var))
  }
  
  cat("\n--- Average Treatment Effect on the Treated (ATT) ---\n")
  print(att_tab, quote = FALSE)
  
  # --- Optional GATE table
  if (!is.null(object$gate) && is.list(object$gate)) {
    
    gate <- object$gate
    groups <- gate$groups
    if (is.null(groups)) groups <- names(gate$tau)
    
    gate_tab <- data.frame(
      Group          = groups,
      Estimate       = round(as.numeric(gate$tau[groups]), digits),
      `Std. Error`   = round(as.numeric(gate$se[groups]), digits),
      `t value`      = round(as.numeric(gate$t[groups]), digits),
      `Pr(>|t|)`     = format.pval(as.numeric(gate$p[groups]), digits = digits),
      `95% CI Lower` = round(as.numeric(gate$ci_lower[groups]), digits),
      `95% CI Upper` = round(as.numeric(gate$ci_upper[groups]), digits),
      check.names = FALSE
    )
    
    cat("\n--- Group Average Treatment Effects (GATE) ---\n")
    cat(sprintf("Grouping variable: %s\n\n", gate$gate_var))
    print(gate_tab, quote = FALSE, row.names = FALSE)
    
    # Add significance stars
    cat("\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  }
  
  cat("\n===============================================\n")
  invisible(object)
}

print.doubleML_did_cs_gate <- function(x, digits = 4, ...) {
  cat("\n===============================================\n")
  cat("  DoubleML DID-CS Estimate\n")
  cat("===============================================\n")
  
  cat("\nOutcome:", x$y_var, "\n")
  cat("Sample: N =", x$N, 
      sprintf("(Treated = %d, Control = %d)\n", x$N_treated, x$N_control))
  
  cat("\nAverage Treatment Effect on the Treated (ATT):\n")
  cat(sprintf("  Estimate: %.*f\n", digits, x$coefficients["ATT"]))
  cat(sprintf("  Std. Error: %.*f\n", digits, x$se))
  cat(sprintf("  t-value: %.*f\n", digits, x$t))
  cat(sprintf("  p-value: %s\n", format.pval(x$p, digits = digits)))
  cat(sprintf("  95%% CI: [%.*f, %.*f]\n", 
              digits, x$ci_lower, digits, x$ci_upper))
  
  # Add significance indicator
  sig <- ifelse(x$p < 0.001, "***",
                ifelse(x$p < 0.01, "**",
                       ifelse(x$p < 0.05, "*",
                              ifelse(x$p < 0.1, ".", ""))))
  if (sig != "") {
    cat(sprintf("  Significance: %s\n", sig))
  }
  
  if (!is.null(x$gate) && is.list(x$gate)) {
    cat("\n--- Group Average Treatment Effects (GATE) ---\n")
    cat(sprintf("Grouping variable: %s\n\n", x$gate$gate_var))
    
    gate <- x$gate
    groups <- gate$groups
    if (is.null(groups)) groups <- names(gate$tau)
    
    for (g in groups) {
      est <- as.numeric(gate$tau[g])
      se <- as.numeric(gate$se[g])
      p <- as.numeric(gate$p[g])
      
      sig_g <- ifelse(p < 0.001, "***",
                      ifelse(p < 0.01, "**",
                             ifelse(p < 0.05, "*",
                                    ifelse(p < 0.1, ".", ""))))
      
      cat(sprintf("  %s: %.*f (SE: %.*f) %s\n", 
                  g, digits, est, digits, se, sig_g))
    }
  }
  
  cat("\n===============================================\n")
  invisible(x)
}








