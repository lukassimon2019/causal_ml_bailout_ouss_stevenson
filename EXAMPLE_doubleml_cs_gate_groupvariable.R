library(dplyr)
library(ranger)
library(glmnet)
library(NuisanceParameters)
library(haven)

CashBail <- read_stata("CashBail.dta") 

#------------------------------------------------------------------------------
#                   Double ML Cross-Sectional DiD - Our function 
#-----------------------------------------------------------------------------


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










#-------------------------------------------------------------------------------
#                     Manipulate the Dataframe 
#-------------------------------------------------------------------------------

CashBail <- CashBail %>%
  mutate(
    yq_numeric = as.numeric(as.factor(yq)),
    # Creating past_fel from prior_fel 
    past_fel = as.numeric(prior_fel > 0),
    
    # Creating Denied 
    Denied = as.numeric(initialbailtype == "Denied"),
    
    # Addressing NAs
    income_missing = as.numeric(is.na(defendantzipmedianhouseholdincom)),
    poverty_missing = as.numeric(is.na(PctBelPov)),
    
    defendantzipmedianhouseholdincom = ifelse(is.na(defendantzipmedianhouseholdincom), 
                                              mean(defendantzipmedianhouseholdincom, na.rm = TRUE), 
                                              defendantzipmedianhouseholdincom),
    
    PctBelPov = ifelse(is.na(PctBelPov), 
                       mean(PctBelPov, na.rm = TRUE), 
                       PctBelPov),
    
    PostElig = EligibleOffense * Post,
    
    off = as.factor(off),
    yq = as.factor(yq),
    DOW = as.factor(DOW), 
    CommissionerName_input = as.factor(CommissionerName_input) 
  )



# --------------------------------------------------------------------
# NOW THE IMPORTANT PART 
# ---------------------------------------------------------------------


## -----------------------------
## 1) Create ONE mutually exclusive race group variable (like the paper)
## -----------------------------
CashBail$racegrp <- NA_character_

CashBail$racegrp[CashBail$defendantisblack == 1 & CashBail$Hisp == 0] <- "Black_nonHisp"
CashBail$racegrp[CashBail$White == 1          & CashBail$Hisp == 0]   <- "White_nonHisp"
CashBail$racegrp[CashBail$Hisp == 1]                                 <- "Hisp"

# keep only those 3 groups
CashBail_race <- subset(CashBail, !is.na(racegrp))
CashBail_race$racegrp <- factor(CashBail_race$racegrp,
                                levels = c("Black_nonHisp","White_nonHisp","Hisp"))

print(table(CashBail_race$racegrp))


## -----------------------------
## 2) Covariates (IMPORTANT: exclude race indicators now)
## -----------------------------
covariates <- c(
  "defendantageatarrest", "male",
  "defendantzipmedianhouseholdincom", "PctBelPov",
  "PD", "felony", "prior_FTA", "prior_9y", "past_fel",
  "income_missing", "poverty_missing"
)

## (optional but safe) drop any covariate that is constant in the full estimation sample
drop_constant_cols <- function(df, vars) {
  vars[sapply(vars, function(v) {
    x <- df[[v]]
    x <- x[!is.na(x)]
    length(unique(x)) > 1
  })]
}
cov_use <- drop_constant_cols(CashBail_race, covariates)


## -----------------------------
## 3) Methods in the format your function expects
## -----------------------------
my_methods <- list(
  D.hat = list(
    rf = list(method = "ranger", tuning = "fold")
  ),
  Y.hat.d = list(
    rf = list(method = "ranger", tuning = "fold")
  )
)


## -----------------------------
## 4) Run one pooled GATT with gate_var = "racegrp"
## -----------------------------
res_gatt <- doubleML_did_cs_gate(
  data = CashBail_race,
  y_var = "ROR",
  d_var = "EligibleOffense",
  t_var = "Post",
  x_vars = cov_use,
  methods = my_methods,
  cf = 5,
  cluster_var = "off",
  heterogeneity = TRUE,
  gate_var = "racegrp"
)

## -----------------------------
## 5) Print GATT results by group
## -----------------------------
print(res_gatt$gate$tau)
print(res_gatt$gate$se)
print(res_gatt$gate$ci_lower)
print(res_gatt$gate$ci_upper)