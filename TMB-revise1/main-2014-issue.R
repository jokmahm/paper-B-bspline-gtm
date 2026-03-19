# ============================================================
# COMPLETE R SCRIPT aligned with gtm_year.cpp
# FIXES INCLUDED:
# 1) p1, p2 fixed per year
# 2) m fixed per year
# 3) scale s_y profiled in C++
# 4) estimated: sigma + spline coeffs
# 5) spline weights normalized in C++
# 6) inactive spline coeffs mapped out in R
# 7) one active spline coefficient fixed as a reference in R
# 8) robust fallback SEs for q_full, s_y, sigma using:
#    pseudo-inverse Hessian + numerical Jacobian (delta method)
# 9) robust naming of q_full[1],...,q_full[K] in summary tables
# ============================================================

library(TMB)
library(splines2)
library(numDeriv)

# ----------------------------
# Load data
# ----------------------------
load("data_N32.RData")  # data_N32: 32 x 4 x 48
load("data_C32.RData")  # data_C32: 32 x 4 x 48
load("data_B.RData")    # object should be data_B

# ----------------------------
# Length grid
# ----------------------------
l_min <- 5
l_max <- 21
step  <- 0.5

x_edges <- seq(l_min, l_max, by = step)
x_mid   <- seq(l_min + step / 2, l_max - step / 2, by = step)
stopifnot(length(x_mid) == 32)

# ----------------------------
# B-spline setup
# ----------------------------
degree <- 2
K <- 32

n_int <- K - degree - 1
int_knots <- seq(min(x_edges), max(x_edges), length.out = n_int + 2)[-c(1, n_int + 2)]

Bfun <- function(x) {
  bSpline(
    x,
    knots = int_knots,
    degree = degree,
    intercept = TRUE,
    Boundary.knots = range(x_edges)
  )
}

grid  <- seq(min(x_edges), max(x_edges), length.out = 4001)
Bgrid <- Bfun(grid)

trapz <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

# I2[k,i] = ∫_{l_i^-}^{l_i^+} B_k(u) du
I2 <- matrix(0, nrow = K, ncol = length(x_mid))
for (i in seq_along(x_mid)) {
  a <- x_mid[i] - step / 2
  b <- x_mid[i] + step / 2
  idx <- which(grid >= a & grid <= b)
  for (k in 1:K) I2[k, i] <- trapz(grid[idx], Bgrid[idx, k])
}

# I1[k,j] = ∫_{l_j^-}^{l_max^+} B_k(u) du
I1 <- matrix(0, nrow = K, ncol = length(x_mid))
max_edge <- max(x_edges)
for (j in seq_along(x_mid)) {
  a <- x_mid[j] - step / 2
  b <- max_edge
  idx <- which(grid >= a & grid <= b)
  for (k in 1:K) I1[k, j] <- trapz(grid[idx], Bgrid[idx, k])
}

stopifnot(all(dim(I1) == c(K, 32)), all(dim(I2) == c(K, 32)))

# ----------------------------
# Compile + load C++
# ----------------------------
dll_name <- dynlib("gtm_year")
try(dyn.unload(dll_name), silent = TRUE)
compile("gtm_year.cpp")
dyn.load(dll_name)

# ----------------------------
# Helper: extract age/year slice
# ----------------------------
get_age_year <- function(arr, age_idx, year) {
  arr[, age_idx, year - 1971]
}

# ----------------------------
# Fixed per-year parameter vectors
# ----------------------------
p1_vec <- c(
  1.85e-01, 1.55e-01, 1.58e-01, 2.09e-01, 7.87e-01, 1.07e+00, 2.65e-01,
  5.88e-01, 7.84e-02, 9.17e-02, 1.01e-01, 1.10e-01, 1.23e-01, 2.32e-01,
  1.96e-01, 2.30e-01, 1.20e+00, 8.83e-02, 8.96e-02, 1.30e-01, 1.70e-01,
  3.77e-01, 7.20e-02, 2.22e-01, 4.29e-01, 1.37e-01, 6.81e-02, 7.35e-01,
  6.18e-02, 2.17e-01, 1.17e-01
)

p2_vec <- c(
  1.48e+01, 1.28e+01, 1.47e+01, 1.28e+01, 1.24e+01, 1.21e+01, 1.38e+01,
  1.30e+01, 1.31e+01, 1.50e+01, 1.39e+01, 1.60e+01, 1.28e+01, 1.36e+01,
  1.24e+01, 1.30e+01, 1.40e+01, 1.38e+01, 1.27e+01, 1.25e+01, 1.22e+01,
  1.43e+01, 1.45e+01, 1.40e+01, 1.46e+01, 1.27e+01, 1.39e+01, 1.22e+01,
  1.58e+01, 1.30e+01, 1.23e+01
)

m_vec <- c(
  1.43e-01, 8.80e-02, 7.52e-02, 8.87e-02, 1.72e-02, 2.20e-02, 9.52e-02,
  8.35e-03, 1.40e-01, 8.92e-03, 9.73e-02, 2.16e-01, 1.09e-01, 1.05e-01,
  9.72e-02, 3.94e-02, 1.86e-02, 4.22e-02, 1.78e-01, 1.48e-01, 1.19e-01,
  1.16e-01, 8.48e-02, 1.23e-01, 9.35e-02, 1.71e-01, 1.02e-01, 1.74e-02,
  1.01e-01, 7.36e-02, 9.22e-02
)

# ----------------------------
# Years
# ----------------------------
year_start  <- 1988
year_finish <- 2017
years <- year_start:(year_finish - 1)

stopifnot(
  length(p1_vec) >= length(years),
  length(p2_vec) >= length(years),
  length(m_vec)  >= length(years)
)

p1_by_year <- setNames(p1_vec[seq_along(years)], as.character(years))
p2_by_year <- setNames(p2_vec[seq_along(years)], as.character(years))
m_by_year  <- setNames(m_vec[seq_along(years)],  as.character(years))

# ============================================================
# Robust helper functions for fallback SE computation
# ============================================================

# ----------------------------
# Stable symmetric pseudo-inverse
# ----------------------------
pinv_sym <- function(H, tol = NULL) {
  H <- 0.5 * (H + t(H))
  ee <- eigen(H, symmetric = TRUE)
  
  if (is.null(tol)) {
    tol <- max(1e-10, sqrt(.Machine$double.eps) * max(abs(ee$values), 1))
  }
  
  keep <- abs(ee$values) > tol
  if (!any(keep)) {
    stop("All Hessian eigenvalues are below tolerance; cannot form pseudo-inverse.")
  }
  
  V <- ee$vectors[, keep, drop = FALSE]
  d <- ee$values[keep]
  
  V %*% diag(1 / d, nrow = length(d)) %*% t(V)
}

# ----------------------------
# Evaluate derived quantities at arbitrary free-parameter vector
# Requires REPORT(q_full), REPORT(s_y), REPORT(sigma) in C++
# ----------------------------
eval_derived <- function(obj, par) {
  obj$fn(par)   # update internal state in obj
  rr <- obj$report()
  
  out <- c(rr$q_full, rr$s_y, rr$sigma)
  names(out) <- c(
    paste0("q_full[", seq_along(rr$q_full), "]"),
    "s_y",
    "sigma"
  )
  out
}

# ----------------------------
# Delta-method fallback SE
# ----------------------------
delta_se_fallback <- function(fit,
                              hess_tol = NULL,
                              zero_clamp = 1e-14) {
  obj <- fit$obj
  par_hat <- fit$opt$par
  
  H <- try(obj$he(par_hat), silent = TRUE)
  if (inherits(H, "try-error") || any(!is.finite(H))) {
    stop("Failed to compute finite Hessian for fallback SE.")
  }
  
  Vtheta <- pinv_sym(H, tol = hess_tol)
  
  gfun <- function(par) eval_derived(obj, par)
  J <- numDeriv::jacobian(gfun, par_hat, method = "Richardson")
  
  Vg <- J %*% Vtheta %*% t(J)
  Vg <- 0.5 * (Vg + t(Vg))
  
  vars <- diag(Vg)
  vars[vars < 0 & abs(vars) < zero_clamp] <- 0
  
  se <- rep(NaN, length(vars))
  ok <- is.finite(vars) & (vars >= 0)
  se[ok] <- sqrt(vars[ok])
  
  est <- gfun(par_hat)
  
  data.frame(
    name = names(est),
    estimate = as.numeric(est),
    se_fallback = se,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# ----------------------------
# Robust extraction:
# - rebuild q_full names if summary() only shows "q_full"
# - patch NaN SE by fallback
# ----------------------------
extract_report_table_robust <- function(fit) {
  s <- summary(fit$rep)
  
  df <- data.frame(
    raw_name = rownames(s),
    estimate = s[, "Estimate"],
    se = s[, "Std. Error"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  df$name <- df$raw_name
  
  q_idx_plain <- which(df$raw_name == "q_full")
  if (length(q_idx_plain) > 0) {
    df$name[q_idx_plain] <- paste0("q_full[", seq_along(q_idx_plain), "]")
  }
  
  keep <- grepl("^q_full\\[", df$name) | df$name %in% c("s_y", "sigma")
  df <- df[keep, c("name", "estimate", "se"), drop = FALSE]
  
  bad <- !is.finite(df$se) | is.nan(df$se)
  
  if (any(bad)) {
    fb <- delta_se_fallback(fit)
    m <- match(df$name, fb$name)
    
    replace_idx <- which(bad & !is.na(m))
    if (length(replace_idx) > 0) {
      df$se[replace_idx] <- fb$se_fallback[m[replace_idx]]
    }
    
    bad_est <- !is.finite(df$estimate)
    replace_est_idx <- which(bad_est & !is.na(m))
    if (length(replace_est_idx) > 0) {
      df$estimate[replace_est_idx] <- fb$estimate[m[replace_est_idx]]
    }
  }
  
  df$z <- df$estimate / df$se
  df$p_value <- 2 * (1 - pnorm(abs(df$z)))
  
  # For coefficients with exactly zero or non-finite SE, do not force z/p
  bad2 <- !is.finite(df$se) | (df$se <= 0)
  df$z[bad2] <- NA_real_
  df$p_value[bad2] <- NA_real_
  
  df
}

# ----------------------------
# Diagnostic: count NaN/invalid SE in raw sdreport output
# ----------------------------
count_bad_se_raw <- function(fit) {
  s0 <- summary(fit$rep)
  df0 <- data.frame(
    raw_name = rownames(s0),
    se = s0[, "Std. Error"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  df0$name <- df0$raw_name
  q_idx_plain <- which(df0$raw_name == "q_full")
  if (length(q_idx_plain) > 0) {
    df0$name[q_idx_plain] <- paste0("q_full[", seq_along(q_idx_plain), "]")
  }
  
  keep <- grepl("^q_full\\[", df0$name) | df0$name %in% c("s_y", "sigma")
  sum(!is.finite(df0$se[keep]) | is.nan(df0$se[keep]))
}

# ============================================================
# Fit one year
# ============================================================
fit_one_year <- function(N_prev_raw, N_obs_raw, C_prev_raw,
                         x_mid, I1, I2,
                         p1_fix, p2_fix, m_fix,
                         eps_pos = 1e-12,
                         eps_denom = 1e-12,
                         obs_threshold = 0,
                         n_start = 20,
                         seed = 1,
                         control = list(eval.max = 5e4, iter.max = 5e4)) {
  
  stopifnot(
    length(N_prev_raw) == 32,
    length(N_obs_raw)  == 32,
    length(C_prev_raw) == 32
  )
  
  K <- nrow(I1)
  stopifnot(K == nrow(I2), ncol(I1) == 32, ncol(I2) == 32)
  
  use_idx <- which(N_obs_raw > obs_threshold)
  if (length(use_idx) == 0) stop("No active bins (all N_obs_raw <= obs_threshold).")
  
  w <- rep(0, 32)
  w[use_idx] <- 1
  
  q_mask <- rep(0, K)
  q_mask[use_idx] <- 1
  
  max_prev <- max(N_prev_raw)
  if (max_prev == 0) max_prev <- 1
  
  max_obs <- max(N_obs_raw)
  if (max_obs == 0) max_obs <- 1
  
  N_prev <- N_prev_raw / max_prev
  N_obs  <- N_obs_raw  / max_obs
  C_prev <- C_prev_raw / max_prev
  
  obs_scale <- sd(N_obs[w > 0])
  if (!is.finite(obs_scale) || obs_scale <= 0) obs_scale <- 0.1
  
  data <- list(
    N_prev = as.numeric(N_prev),
    N_obs  = as.numeric(N_obs),
    C_prev = as.numeric(C_prev),
    x_mid  = as.numeric(x_mid),
    I1     = I1,
    I2     = I2,
    eps_pos   = as.numeric(eps_pos),
    eps_denom = as.numeric(eps_denom),
    w = as.numeric(w),
    obs_scale = as.numeric(obs_scale),
    q_mask = as.numeric(q_mask),
    p1_fix = as.numeric(p1_fix),
    p2_fix = as.numeric(p2_fix),
    m_fix  = as.numeric(m_fix)
  )
  
  parameters0 <- list(
    log_sigma  = 0,
    q_raw_full = rep(log(0.1), K)
  )
  
  build_obj <- function(par_list) {
    active <- which(q_mask == 1)
    if (length(active) == 0) stop("No active spline coefficients.")
    
    ref_idx  <- active[1]
    free_idx <- setdiff(active, ref_idx)
    
    q_map <- rep(NA_integer_, length(par_list$q_raw_full))
    if (length(free_idx) > 0) {
      q_map[free_idx] <- seq_along(free_idx)
    }
    q_map <- factor(q_map)
    
    par_list$q_raw_full[ref_idx] <- 0
    
    MakeADFun(
      data = data,
      parameters = par_list,
      DLL = "gtm_year",
      silent = TRUE,
      map = list(q_raw_full = q_map)
    )
  }
  
  set.seed(seed)
  best <- NULL
  best_objval <- Inf
  
  for (s in seq_len(n_start)) {
    par_s <- parameters0
    par_s$log_sigma  <- rnorm(1, 0, 1)
    par_s$q_raw_full <- rnorm(K, mean = log(0.1), sd = 1)
    
    active   <- which(q_mask == 1)
    ref_idx  <- active[1]
    inactive <- which(q_mask == 0)
    
    if (length(inactive) > 0) par_s$q_raw_full[inactive] <- log(1e-8)
    par_s$q_raw_full[ref_idx] <- 0
    
    obj_s <- build_obj(par_s)
    
    opt_s <- try(
      nlminb(
        start     = obj_s$par,
        objective = obj_s$fn,
        gradient  = obj_s$gr,
        control   = control
      ),
      silent = TRUE
    )
    
    if (inherits(opt_s, "try-error")) next
    
    if (is.finite(opt_s$objective) && opt_s$objective < best_objval) {
      best <- list(obj = obj_s, opt = opt_s)
      best_objval <- opt_s$objective
    }
  }
  
  if (is.null(best)) stop("All multi-start attempts failed for this year.")
  
  obj <- best$obj
  opt <- best$opt
  
  H <- try(obj$he(opt$par), silent = TRUE)
  
  rep <- if (inherits(H, "try-error") || any(!is.finite(H))) {
    warning("AD Hessian failed; trying sdreport without supplied Hessian.")
    sdreport(obj, par.fixed = opt$par)
  } else {
    sdreport(obj, par.fixed = opt$par, hessian.fixed = H)
  }
  
  list(
    opt = opt,
    rep = rep,
    obj = obj,
    obs_scale = obs_scale,
    use_idx = use_idx,
    q_mask = q_mask,
    p1_fix = p1_fix,
    p2_fix = p2_fix,
    m_fix  = m_fix
  )
}

# ============================================================
# Run across years
# ============================================================
fits <- vector("list", length(years))
names(fits) <- as.character(years)

for (ii in seq_along(years)) {
  y <- years[ii]
  
  N_prev_raw <- get_age_year(data_N32, 2, y)
  N_obs_raw  <- get_age_year(data_N32, 3, y + 1)
  C_prev_raw <- get_age_year(data_C32, 2, y + 1)
  
  fits[[ii]] <- fit_one_year(
    N_prev_raw = N_prev_raw,
    N_obs_raw  = N_obs_raw,
    C_prev_raw = C_prev_raw,
    x_mid = x_mid,
    I1 = I1,
    I2 = I2,
    p1_fix = p1_by_year[as.character(y)],
    p2_fix = p2_by_year[as.character(y)],
    m_fix  = m_by_year[as.character(y)],
    obs_threshold = 1e-12,
    n_start = 20,
    seed = 1000 + ii,
    control = list(eval.max = 5e4, iter.max = 5e4)
  )
  
  max_grad  <- max(abs(fits[[ii]]$obj$gr(fits[[ii]]$opt$par)))
  conv_code <- fits[[ii]]$opt$convergence
  msg       <- fits[[ii]]$opt$message
  ok        <- (conv_code == 0) || (max_grad < 1e-8)
  
  cat(
    "Year", y,
    "ok:", ok,
    "conv_code:", conv_code,
    "max|grad|:", signif(max_grad, 6),
    "message:", msg, "\n"
  )
}

# ============================================================
# Raw NaN-SE counts before fallback
# ============================================================
for (yy in names(fits)) {
  bad0 <- count_bad_se_raw(fits[[yy]])
  cat("Year", yy, "- NaN/invalid SE count before fallback:", bad0, "\n")
}

# ============================================================
# Build result tables for all years safely
# ============================================================
tabs <- list()

for (yy in names(fits)) {
  cat("Processing year", yy, "...\n")
  
  tabs[[yy]] <- tryCatch(
    extract_report_table_robust(fits[[yy]]),
    error = function(e) {
      cat("Year", yy, "failed:", conditionMessage(e), "\n")
      NULL
    }
  )
}

# Show failed years, if any
failed_years <- names(tabs)[vapply(tabs, is.null, logical(1))]
cat("Failed years:", if (length(failed_years) == 0) "None" else paste(failed_years, collapse = ", "), "\n")

# ============================================================
# Remaining bad SE counts after fallback
# ============================================================
for (yy in names(tabs)) {
  if (is.null(tabs[[yy]])) next
  n_bad <- sum(!is.finite(tabs[[yy]]$se) | is.nan(tabs[[yy]]$se))
  cat("Year", yy, "- remaining bad SE after fallback:", n_bad, "\n")
}

# ============================================================
# Example extraction for one year
# ============================================================
tab_1988 <- tabs[["1988"]]
print(head(tab_1988, 40))

sigma_1988 <- tab_1988$estimate[tab_1988$name == "sigma"]
s_y_1988   <- tab_1988$estimate[tab_1988$name == "s_y"]

cat("1988: s_y =", s_y_1988, " sigma =", sigma_1988, "\n")

rep_list_1988 <- fits[["1988"]]$obj$report()
G_1988       <- rep_list_1988$G
N_pred0_1988 <- rep_list_1988$N_pred0
N_pred_1988  <- rep_list_1988$N_pred
inside_1988  <- rep_list_1988$inside

# ============================================================
# Combine all years in one long table
# ============================================================
all_results <- do.call(
  rbind,
  lapply(names(tabs), function(yy) {
    if (is.null(tabs[[yy]])) return(NULL)
    out <- tabs[[yy]]
    out$year <- as.integer(yy)
    out
  })
)

all_results <- all_results[, c("year", "name", "estimate", "se", "z", "p_value")]
print(head(all_results, 50))

# ============================================================
# Optional helpers
# ============================================================

# Get one year's robust table
get_year_table <- function(year) {
  tabs[[as.character(year)]]
}

# Save all results if desired
# write.csv(all_results, "gtm_all_results.csv", row.names = FALSE)


# ============================================================
# Build N3 and GT2 from fitted TMB objects
# ============================================================

build_result_arrays <- function(fits, years_fit) {
  n_years <- length(years_fit)
  n_class <- 32
  
  N3  <- matrix(NA_real_, nrow = n_class, ncol = n_years)
  GT2 <- array(NA_real_, dim = c(n_class, n_class, n_years))
  
  colnames(N3) <- as.character(years_fit + 1)   # observed year
  dimnames(GT2) <- list(
    paste0("to_", 1:n_class),
    paste0("from_", 1:n_class),
    as.character(years_fit)
  )
  
  for (j in seq_along(years_fit)) {
    fit_year <- years_fit[j]
    rr <- fits[[as.character(fit_year)]]$obj$report()
    
    N3[, j]    <- as.numeric(rr$N_pred)
    GT2[, , j] <- rr$G
  }
  
  list(N3 = N3, GT2 = GT2)
}

# Example:
res_arrays <- build_result_arrays(fits, years)
N3  <- res_arrays$N3
GT2 <- res_arrays$GT2

# ============================================================
# Functon to campair years
# ============================================================

testGPT_R2 <- function(start_obs, year_obs, x, N3, data_N32) {
  col_idx <- year_obs - start_obs + 1
  arr_idx <- year_obs - 1971
  
  if (col_idx < 1 || col_idx > ncol(N3)) {
    stop("Requested year_obs is outside the columns of N3.")
  }
  if (arr_idx < 1 || arr_idx > dim(data_N32)[3]) {
    stop("Requested year_obs is outside data_N32 year range.")
  }
  
  n3_true <- data_N32[, 3, arr_idx]
  n3_true <- n3_true / max(n3_true)
  
  n3_est <- N3[, col_idx]
  dx <- if (length(x) > 1) min(diff(x)) else 0.5
  
  plot(
    x, n3_true,
    type = "h",
    lwd = 10,
    lend = "butt",
    col = "grey75",
    xlab = "Length",
    ylab = "Normalized fish abundance",
    ylim = c(0, 1.05),
    main = paste("Year:", year_obs)
  )
  
  lines(x, n3_est, lwd = 2, col = "blue")
  
  legend(
    "topright",
    legend = c("Empirical data", "Estimated by Model"),
    col = c("grey75", "blue"),
    lwd = c(8, 2),
    bty = "n"
  )
  
  box(bty = "l")
}


obs_years <- years + 1
start_obs <- min(obs_years)   # 1989

testGPT_R2(
  start_obs = start_obs,
  year_obs = 2015,
  x = x_mid,
  N3 = N3,
  data_N32 = data_N32
)





library(ggplot2)
library(reshape2)



# ============================================================
# k_mean_GTM2_R
# - clusters yearly GTMs using row sums of each GTM
# - makes elbow plot
# - makes cluster profile plots
# - overlays cluster labels on biomass time series
# ============================================================
k_mean_GTM2_R_better <- function(GT2,
                                 years,
                                 data_B,
                                 optimal_k = 2,
                                 seed = 100,
                                 nstart = 50,
                                 save_plots = FALSE,
                                 out_dir = ".") {
  
  stopifnot(length(dim(GT2)) == 3)
  stopifnot(dim(GT2)[1] == 32, dim(GT2)[2] == 32)
  stopifnot(dim(GT2)[3] == length(years))
  
  n3 <- dim(GT2)[3]
  
  GT <- matrix(0, nrow = n3, ncol = 32)
  for (i in seq_len(n3)) {
    GT[i, ] <- rowSums(GT2[, , i])
  }
  
  # ----------------------------
  # Elbow method
  # ----------------------------
  k_values <- 1:min(10, n3)
  wss <- numeric(length(k_values))
  
  set.seed(seed)
  for (ii in seq_along(k_values)) {
    k <- k_values[ii]
    km <- kmeans(GT, centers = k, nstart = nstart)
    wss[ii] <- sum(km$withinss)
  }
  
  df_elbow <- data.frame(k = k_values, wss = wss)
  
  p_elbow <- ggplot(df_elbow, aes(x = k, y = wss)) +
    geom_line(linewidth = 0.9, color = col_model) +
    geom_point(size = 2.8, color = col_model) +
    scale_x_continuous(breaks = k_values) +
    labs(
      x = "Number of clusters (k)",
      y = "Within-cluster sum of squares"
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )
  
  print(p_elbow)
  
  if (save_plots) {
    ggsave(
      filename = file.path(out_dir, "elbow_plot_gtm.png"),
      plot = p_elbow, width = 7, height = 5, dpi = 300
    )
  }
  
  # ----------------------------
  # Final kmeans
  # ----------------------------
  set.seed(seed)
  km_final <- kmeans(GT, centers = optimal_k, nstart = nstart)
  idx <- km_final$cluster
  centers <- km_final$centers
  
  # ----------------------------
  # Cluster profile plots
  # ----------------------------
  for (cl in seq_len(optimal_k)) {
    cluster_indices <- which(idx == cl)
    cluster_mat <- GT[cluster_indices, , drop = FALSE]
    
    df_cluster <- as.data.frame(cluster_mat)
    colnames(df_cluster) <- paste0("L", 1:32)
    df_cluster$year <- years[cluster_indices]
    
    df_long <- reshape2::melt(
      df_cluster,
      id.vars = "year",
      variable.name = "length_class",
      value.name = "row_sum"
    )
    df_long$length_class_num <- as.numeric(gsub("L", "", df_long$length_class))
    
    mean_profile <- colMeans(cluster_mat)
    df_mean <- data.frame(
      length_class_num = 1:32,
      mean_row_sum = mean_profile
    )
    
    max_indices <- apply(cluster_mat, 1, which.max)
    avg_index <- mean(max_indices)
    print(avg_index)
    
    cluster_color <- if (cl == 1) col_cluster1 else col_cluster2
    
    p_cluster <- ggplot() +
      geom_line(
        data = df_long,
        aes(x = length_class_num, y = row_sum, group = factor(year)),
        linewidth = 0.8, alpha = 0.7, color = col_grayline
      ) +
      geom_line(
        data = df_mean,
        aes(x = length_class_num, y = mean_row_sum),
        linewidth = 1.3, color = cluster_color
      ) +
      geom_vline(
        xintercept = avg_index,
        linetype = "dashed",
        linewidth = 0.8,
        color = col_vline
      ) +
      scale_x_continuous(breaks = seq(0, 32, by = 4), limits = c(1, 32)) +
      labs(
        title = paste("Cluster", cl),
        x = "Length class",
        y = "CGTR"
      ) +
      theme_classic(base_size = 13) +
      theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        legend.position = "none"
      )
    
    print(p_cluster)
    
    if (save_plots) {
      ggsave(
        filename = file.path(out_dir, paste0("cluster_profile_", cl, ".png")),
        plot = p_cluster, width = 8, height = 5, dpi = 300
      )
    }
  }
  
  # ----------------------------
  # Biomass with cluster labels
  # ----------------------------
  b_idx <- years - 1971 - 1
  
  df_bio <- data.frame(
    year = years,
    age3_biomass = as.numeric(data_B[3, b_idx]),
    total_biomass = as.numeric(data_B[5, b_idx]),
    cluster = factor(idx)
  )
  
  p_age3 <- ggplot(df_bio, aes(x = year, y = age3_biomass)) +
    geom_line(linewidth = 0.9, color = col_biomass) +
    geom_point(size = 2.6, color = col_biomass) +
    geom_text(aes(label = cluster), vjust = -0.8, size = 4.5, fontface = "bold", color = col_vline) +
    labs(
      x = "Year",
      y = "Age-3 biomass (t)"
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )
  
  print(p_age3)
  
  if (save_plots) {
    ggsave(
      filename = file.path(out_dir, "age3_biomass_clusters.png"),
      plot = p_age3, width = 9, height = 5, dpi = 300
    )
  }
  
  p_total <- ggplot(df_bio, aes(x = year, y = total_biomass)) +
    geom_line(linewidth = 0.9, color = col_biomass) +
    geom_point(size = 2.6, color = col_biomass) +
    geom_text(aes(label = cluster), vjust = -0.8, size = 4.5, fontface = "bold", color = col_vline) +
    labs(
      x = "Year",
      y = "Total biomass (t)"
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11)
    )
  
  print(p_total)
  
  if (save_plots) {
    ggsave(
      filename = file.path(out_dir, "total_biomass_clusters.png"),
      plot = p_total, width = 9, height = 5, dpi = 300
    )
  }
  
  invisible(list(
    idx = idx,
    centers = centers,
    GT = GT,
    k_values = k_values,
    wss = wss
  ))
}

dir.create("results_plots", showWarnings = FALSE)
# ----------------------------
# Color palette
# ----------------------------
col_empirical <- "#B8C4D6"   # soft blue-gray
col_model     <- "#2C7FB8"   # medium blue

col_cluster1  <- "#CC79A7"   # soft magenta
col_cluster2  <- "#D55E00"   # warm orange-red

col_biomass   <- "#1F3B73"   # dark blue
col_vline     <- "#C44E52"   # muted red
col_grayline  <- "#B0B0B0"   # neutral gray

##1) K-means clustering plots
km_out <- k_mean_GTM2_R_better(
  GT2 = GT2,
  years = years,
  data_B = data_B,
  optimal_k = 2,
  seed = 100,
  save_plots = TRUE,
  out_dir = "results_plots"
)






# ============================================================
# Maturity function in R
# ============================================================
maturity_vec_R <- function(x_mid, p1, p2) {
  1 / (1 + exp(4 * p1 * (p2 - x_mid)))
}

# ============================================================
# Project one year using sampled fit parameters and evaluation-year inputs
# This mirrors your current TMB logic:
#   Surv = N_prev * (1-rmat)
#   inside = exp(-12m)*Surv - exp(-6m)*C_prev
#   N_pred0 = G %*% inside
#   N_pred  = s_y * N_pred0
# ============================================================
project_with_sampled_fit <- function(sample_fit,
                                     N_prev_raw,
                                     C_prev_raw,
                                     x_mid,
                                     eps_pos = 1e-12) {
  
  rr <- sample_fit$obj$report()
  G <- rr$G
  s_y <- as.numeric(rr$s_y)
  
  p1 <- as.numeric(sample_fit$p1_fix)
  p2 <- as.numeric(sample_fit$p2_fix)
  m  <- as.numeric(sample_fit$m_fix)
  
  max_prev <- max(N_prev_raw)
  if (!is.finite(max_prev) || max_prev <= 0) max_prev <- 1
  
  N_prev <- N_prev_raw / max_prev
  C_prev <- C_prev_raw / max_prev
  
  rmat <- maturity_vec_R(x_mid, p1, p2)
  Surv <- N_prev * (1 - rmat)
  inside <- exp(-m) * Surv - exp(-0.5 * m) * C_prev
  
  N_pred0 <- as.numeric(G %*% inside)
  N_pred0 <- pmax(N_pred0, eps_pos)
  N_pred  <- s_y * N_pred0
  
  N_pred
}



# ============================================================
# monte_carlo_compare_R
# - cluster_idx comes from k_mean_GTM2_R(...)$idx
# - uses sampled FULL fitted year objects, not only G
# - compares normalized model prediction with normalized A3
# ============================================================
monte_carlo_compare_R_better <- function(cluster_idx,
                                         fits,
                                         years_fit,
                                         data_N32,
                                         data_C32,
                                         data_B,
                                         x_mid,
                                         num_samples = 100,
                                         seed = 12345,
                                         save_plots = FALSE,
                                         out_dir = ".",
                                         file_prefix = "mc") {
  
  if (is.list(cluster_idx) && !is.null(cluster_idx$idx)) {
    cluster_idx <- cluster_idx$idx
  }
  cluster_idx <- as.integer(cluster_idx)
  
  stopifnot(length(cluster_idx) == length(years_fit))
  
  idx_cluster1 <- which(cluster_idx == 1)
  idx_cluster2 <- which(cluster_idx == 2)
  
  if (length(idx_cluster1) == 0 || length(idx_cluster2) == 0) {
    stop("Both cluster 1 and cluster 2 must contain at least one year.")
  }
  
  # ----------------------------
  # Colors
  # ----------------------------
  col_biomass <- "blue"
  col_c1 <- "magenta"
  col_c2 <- "red"
  fill_c1 <- adjustcolor(col_c1, alpha.f = 0.18)
  fill_c2 <- adjustcolor(col_c2, alpha.f = 0.18)
  
  y_values <- years_fit
  num_years <- length(y_values)
  
  N_l <- array(NA_real_, dim = c(32, num_years, num_samples))
  N_h <- array(NA_real_, dim = c(32, num_years, num_samples))
  A3  <- matrix(NA_real_, nrow = 32, ncol = num_years)
  
  norm_diff_Nl_samples <- matrix(NA_real_, nrow = num_samples, ncol = num_years)
  norm_diff_Nh_samples <- matrix(NA_real_, nrow = num_samples, ncol = num_years)
  
  set.seed(seed)
  
  for (i in seq_along(y_values)) {
    y <- y_values[i]
    
    # same indexing as your fitted model
    N_prev_raw <- data_N32[, 2, y - 1971]
    N_obs_raw  <- data_N32[, 3, (y + 1) - 1971]
    C_prev_raw <- data_C32[, 2, (y + 1) - 1971]
    
    max_obs <- max(N_obs_raw)
    if (!is.finite(max_obs) || max_obs <= 0) max_obs <- 1
    A3[, i] <- N_obs_raw / max_obs
    
    for (j in seq_len(num_samples)) {
      rand_idx_l <- sample(idx_cluster1, size = 1)
      rand_idx_h <- sample(idx_cluster2, size = 1)
      
      fit_l <- fits[[as.character(years_fit[rand_idx_l])]]
      fit_h <- fits[[as.character(years_fit[rand_idx_h])]]
      
      N_l[, i, j] <- project_with_sampled_fit(
        sample_fit = fit_l,
        N_prev_raw = N_prev_raw,
        C_prev_raw = C_prev_raw,
        x_mid = x_mid
      )
      
      N_h[, i, j] <- project_with_sampled_fit(
        sample_fit = fit_h,
        N_prev_raw = N_prev_raw,
        C_prev_raw = C_prev_raw,
        x_mid = x_mid
      )
      
      norm_diff_Nl_samples[j, i] <- sqrt(sum((A3[, i] - N_l[, i, j])^2))
      norm_diff_Nh_samples[j, i] <- sqrt(sum((A3[, i] - N_h[, i, j])^2))
    }
  }
  
  # ----------------------------
  # Biomass aligned with fit years
  # ----------------------------
  b_idx <- y_values - 1971
  total_biomass <- as.numeric(data_B[5, b_idx])
  age3_biomass  <- as.numeric(data_B[3, b_idx])
  
  norm_total_biomass <- total_biomass / max(total_biomass)
  norm_age3_biomass  <- age3_biomass  / max(age3_biomass)
  
  # ----------------------------
  # Normalize errors globally
  # ----------------------------
  max_norm_diffs <- max(c(norm_diff_Nl_samples, norm_diff_Nh_samples), na.rm = TRUE)
  norm_diff_Nl_samples <- norm_diff_Nl_samples / max_norm_diffs
  norm_diff_Nh_samples <- norm_diff_Nh_samples / max_norm_diffs
  
  mean_Nl <- colMeans(norm_diff_Nl_samples, na.rm = TRUE)
  mean_Nh <- colMeans(norm_diff_Nh_samples, na.rm = TRUE)
  std_Nl  <- apply(norm_diff_Nl_samples, 2, sd, na.rm = TRUE)
  std_Nh  <- apply(norm_diff_Nh_samples, 2, sd, na.rm = TRUE)
  
  # ---------------------------------------------------------
  # IMPORTANT:
  # MATLAB code used prctile(mean +/- std, ...) which on a
  # vector gives a single scalar. That is not ideal for a band.
  # Here we use pointwise empirical quantiles, which gives
  # a proper year-by-year shaded region.
  # ---------------------------------------------------------
  lower_Nl <- apply(norm_diff_Nl_samples, 2, quantile, probs = 0.025, na.rm = TRUE)
  upper_Nl <- apply(norm_diff_Nl_samples, 2, quantile, probs = 0.975, na.rm = TRUE)
  lower_Nh <- apply(norm_diff_Nh_samples, 2, quantile, probs = 0.025, na.rm = TRUE)
  upper_Nh <- apply(norm_diff_Nh_samples, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  # ----------------------------
  # Internal plotting helper
  # ----------------------------
  draw_dual_axis_plot <- function(x_years, left_y, left_label, main_title,
                                  mean1, mean2, low1, up1, low2, up2,
                                  file = NULL) {
    
    if (!is.null(file)) {
      png(file, width = 1100, height = 650, res = 130)
    }
    
    old_par <- par(no.readonly = TRUE)
    on.exit({
      par(old_par)
      if (!is.null(file)) dev.off()
    }, add = TRUE)
    
    par(
      mar = c(4.5, 5.2, 3.0, 5.2),
      mgp = c(2.6, 0.8, 0),
      tcl = -0.3,
      xpd = FALSE
    )
    
    # Left axis: biomass
    plot(
      x_years, left_y,
      type = "o",
      pch = 1,
      lwd = 1.6,
      col = col_biomass,
      xlab = "Year",
      ylab = "",
      ylim = c(0, 1),
      main = main_title,
      cex.main = 1.15,
      cex.lab = 1.0,
      cex.axis = 0.9,
      axes = FALSE
    )
    
    axis(1)
    axis(2, col.axis = col_biomass)
    
    mtext(
      left_label,
      side = 2,
      line = 3.2,
      col = col_biomass,
      cex = 1
    )
    box()
    grid(nx = NA, ny = NULL, col = "grey85", lty = 1)
    
    # Right axis: errors
    par(new = TRUE)
    y_right_lim <- range(c(low1, up1, low2, up2, mean1, mean2), na.rm = TRUE)
    
    plot(
      x_years, mean1,
      type = "n",
      axes = FALSE,
      xlab = "",
      ylab = "",
      ylim = y_right_lim
    )
    
    polygon(
      c(x_years, rev(x_years)),
      c(up1, rev(low1)),
      col = fill_c1,
      border = NA
    )
    polygon(
      c(x_years, rev(x_years)),
      c(up2, rev(low2)),
      col = fill_c2,
      border = NA
    )
    
    lines(x_years, mean1, col = col_c1, lty = 2, lwd = 1.2)
    points(x_years, mean1, col = col_c1, pch = 15, cex = 0.8)
    
    lines(x_years, mean2, col = col_c2, lty = 2, lwd = 1.2)
    points(x_years, mean2, col = col_c2, pch = 17, cex = 0.8)
    
    axis(4, col.axis = col_c2)
    mtext("Normalized Error in Norm", side = 4, line = 3.2, col = col_c2)
    
    legend(
      "topright",
      inset = c(0.01, 0.01),
      legend = c(left_label, "Err(cluster 1)", "Err(cluster 2)"),
      col = c(col_biomass, col_c1, col_c2),
      lty = c(1, 2, 2),
      lwd = c(1.6, 1.2, 1.2),
      pch = c(1, 15, 17),
      pt.cex = 0.8,
      cex = 0.82,
      bty = "o",
      bg = "white"
    )
  }
  
  # ----------------------------
  # Plot 1: Total biomass
  # ----------------------------
  file1 <- if (save_plots) file.path(out_dir, paste0(file_prefix, "_total_biomass_vs_error.png")) else NULL
  draw_dual_axis_plot(
    x_years = y_values,
    left_y = norm_total_biomass,
    left_label = "Normalized Total Biomass",
    main_title = "Total Biomass vs Cluster-based Error",
    mean1 = mean_Nl,
    mean2 = mean_Nh,
    low1 = lower_Nl,
    up1 = upper_Nl,
    low2 = lower_Nh,
    up2 = upper_Nh,
    file = file1
  )
  
  # ----------------------------
  # Plot 2: Age-3 biomass
  # ----------------------------
  file2 <- if (save_plots) file.path(out_dir, paste0(file_prefix, "_age3_biomass_vs_error.png")) else NULL
  draw_dual_axis_plot(
    x_years = y_values,
    left_y = norm_age3_biomass,
    left_label = "Normalized Age-3 Biomass",
    main_title = "Age-3 Biomass vs Cluster-based Error",
    mean1 = mean_Nl,
    mean2 = mean_Nh,
    low1 = lower_Nl,
    up1 = upper_Nl,
    low2 = lower_Nh,
    up2 = upper_Nh,
    file = file2
  )
  
  invisible(list(
    N_l = N_l,
    N_h = N_h,
    A3 = A3,
    norm_diff_Nl_samples = norm_diff_Nl_samples,
    norm_diff_Nh_samples = norm_diff_Nh_samples,
    mean_Nl = mean_Nl,
    mean_Nh = mean_Nh,
    std_Nl = std_Nl,
    std_Nh = std_Nh,
    lower_Nl = lower_Nl,
    upper_Nl = upper_Nl,
    lower_Nh = lower_Nh,
    upper_Nh = upper_Nh
  ))
}

idx <- km_out$idx
# 2) Monte Carlo comparison
mc_out <- monte_carlo_compare_R_better(
  cluster_idx = idx,
  fits = fits,
  years_fit = years,
  data_N32 = data_N32,
  data_C32 = data_C32,
  data_B = data_B,
  x_mid = x_mid,
  num_samples = 100,
  seed = 12345,
  save_plots = TRUE,
  out_dir = "results_plots",
  file_prefix = "mc"
)

