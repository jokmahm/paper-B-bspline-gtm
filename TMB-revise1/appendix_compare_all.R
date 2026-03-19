library(TMB)
library(splines2)
library(numDeriv)

# ============================================================
# Load data
# ============================================================
load("data_N32.RData")
load("data_C32.RData")
load("data_B.RData")

# ============================================================
# Grid
# ============================================================
l_min <- 5
l_max <- 21
step  <- 0.5

x_edges <- seq(l_min, l_max, by = step)
x_mid   <- seq(l_min + step/2, l_max - step/2, by = step)

stopifnot(length(x_mid) == 32, length(x_edges) == 33)

# ============================================================
# B-spline matrices
# ============================================================
degree <- 2
K_bs <- 32

n_int <- K_bs - degree - 1
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

I2 <- matrix(0, nrow = K_bs, ncol = length(x_mid))
for (i in seq_along(x_mid)) {
  a <- x_mid[i] - step/2
  b <- x_mid[i] + step/2
  idx <- which(grid >= a & grid <= b)
  for (k in 1:K_bs) I2[k, i] <- trapz(grid[idx], Bgrid[idx, k])
}

I1 <- matrix(0, nrow = K_bs, ncol = length(x_mid))
max_edge <- max(x_edges)
for (j in seq_along(x_mid)) {
  a <- x_mid[j] - step/2
  b <- max_edge
  idx <- which(grid >= a & grid <= b)
  for (k in 1:K_bs) I1[k, j] <- trapz(grid[idx], Bgrid[idx, k])
}

# ============================================================
# Compile and load DLLs
# ============================================================
compile_and_load <- function(base, ext = ".cpp") {
  dll <- dynlib(base)
  try(dyn.unload(dll), silent = TRUE)
  compile(paste0(base, ext))
  dyn.load(dll)
}

# B-spline uses the fixed file for 2014 -> 2015
compile_and_load("gtm_bspline_year", ext = ".cpp")

# Comparator methods unchanged
compile_and_load("gtm_normal_year", ext = ".cpp")
compile_and_load("gtm_gamma_year", ext = ".cpp")
compile_and_load("gtm_lognormal_year", ext = ".cpp")

# ============================================================
# Helper
# ============================================================
get_age_year <- function(arr, age_idx, year) {
  arr[, age_idx, year - 1971]
}

# ============================================================
# Fixed per-year parameters
# ============================================================
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

year_start  <- 1988
year_finish <- 2017
years <- year_start:(year_finish - 1)

p1_by_year <- setNames(p1_vec[seq_along(years)], as.character(years))
p2_by_year <- setNames(p2_vec[seq_along(years)], as.character(years))
m_by_year  <- setNames(m_vec[seq_along(years)],  as.character(years))

# ============================================================
# B-spline fit function using gtm_bspline_year
# Includes 2014 -> 2015 fix only
# ============================================================
fit_one_year_bspline <- function(N_prev_raw, N_obs_raw, C_prev_raw,
                                 x_mid, I1, I2,
                                 p1_fix, p2_fix, m_fix,
                                 year_y,
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
  
  # Apply the fix only for y = 2014 (projection to 2015)
  use_support_mask <- (year_y == 2014)
  lambda_smooth    <- if (year_y == 2014) 0.01 else 0
  
  if (!use_support_mask) {
    q_mask <- rep(0, K)
    q_mask[use_idx] <- 1
  } else {
    support_tol <- 1e-12
    q_mask <- rep(0, K)
    for (k in 1:K) {
      contributes <- any(I2[k, use_idx] > support_tol) || any(I1[k, use_idx] > support_tol)
      q_mask[k] <- as.numeric(contributes)
    }
  }
  
  max_prev <- max(N_prev_raw)
  if (max_prev <= 0) max_prev <- 1
  
  max_obs <- max(N_obs_raw)
  if (max_obs <= 0) max_obs <- 1
  
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
    m_fix  = as.numeric(m_fix),
    lambda_smooth = as.numeric(lambda_smooth)
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
      DLL = "gtm_bspline_year",
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
  
  if (is.null(best)) stop("All B-spline multi-start attempts failed for this year.")
  
  H <- try(best$obj$he(best$opt$par), silent = TRUE)
  
  rep <- if (inherits(H, "try-error") || any(!is.finite(H))) {
    sdreport(best$obj, par.fixed = best$opt$par, getJointPrecision = FALSE)
  } else {
    sdreport(best$obj,
             par.fixed = best$opt$par,
             hessian.fixed = H,
             getJointPrecision = FALSE)
  }
  
  list(
    obj = best$obj,
    opt = best$opt,
    rep = rep,
    method = "bspline",
    lambda_smooth = lambda_smooth,
    use_support_mask = use_support_mask
  )
}

# ============================================================
# Generic comparator builder
# ============================================================
build_comparator_obj <- function(method,
                                 N_prev_raw, N_obs_raw, C_prev_raw,
                                 x_mid, x_edges,
                                 p1_fix, p2_fix, m_fix,
                                 obs_threshold = 1e-12,
                                 eps_pos = 1e-12,
                                 eps_shape = 1e-10,
                                 start_par = NULL) {
  
  use_idx <- which(N_obs_raw > obs_threshold)
  if (length(use_idx) == 0) stop("No active bins.")
  
  w <- rep(0, length(N_obs_raw))
  w[use_idx] <- 1
  
  max_prev <- max(N_prev_raw)
  if (max_prev <= 0) max_prev <- 1
  
  max_obs <- max(N_obs_raw)
  if (max_obs <= 0) max_obs <- 1
  
  N_prev <- N_prev_raw / max_prev
  N_obs  <- N_obs_raw  / max_obs
  C_prev <- C_prev_raw / max_prev
  
  obs_scale <- sd(N_obs[w > 0])
  if (!is.finite(obs_scale) || obs_scale <= 0) obs_scale <- 0.1
  
  if (method == "normal") {
    DLL <- "gtm_normal_year"
    if (is.null(start_par)) {
      start_par <- list(
        log_sigma = 0,
        u1 = log(1.05 * max(x_mid) - max(x_mid)),
        u2 = log(0.3),
        u3 = log(1.0),
        u4 = log(0.1),
        u5 = 0
      )
    }
    data <- list(
      N_prev = as.numeric(N_prev),
      N_obs  = as.numeric(N_obs),
      C_prev = as.numeric(C_prev),
      x_mid  = as.numeric(x_mid),
      x_edges = as.numeric(x_edges),
      w = as.numeric(w),
      obs_scale = as.numeric(obs_scale),
      eps_pos = as.numeric(eps_pos),
      eps_var = as.numeric(eps_shape),
      p1_fix = as.numeric(p1_fix),
      p2_fix = as.numeric(p2_fix),
      m_fix  = as.numeric(m_fix)
    )
  }
  
  if (method == "gamma") {
    DLL <- "gtm_gamma_year"
    if (is.null(start_par)) {
      start_par <- list(
        log_sigma = 0,
        z_Linf = 0,
        z_K = 0,
        z_CV = 0
      )
    }
    data <- list(
      N_prev = as.numeric(N_prev),
      N_obs  = as.numeric(N_obs),
      C_prev = as.numeric(C_prev),
      x_mid  = as.numeric(x_mid),
      x_edges = as.numeric(x_edges),
      w = as.numeric(w),
      obs_scale = as.numeric(obs_scale),
      eps_pos = as.numeric(eps_pos),
      eps_mu = as.numeric(eps_shape),
      p1_fix = as.numeric(p1_fix),
      p2_fix = as.numeric(p2_fix),
      m_fix  = as.numeric(m_fix)
    )
  }
  
  if (method == "lognormal") {
    DLL <- "gtm_lognormal_year"
    if (is.null(start_par)) {
      start_par <- list(
        log_sigma = 0,
        z_Linf = 0,
        z_K = 0,
        z_CV = 0
      )
    }
    data <- list(
      N_prev = as.numeric(N_prev),
      N_obs  = as.numeric(N_obs),
      C_prev = as.numeric(C_prev),
      x_mid  = as.numeric(x_mid),
      x_edges = as.numeric(x_edges),
      w = as.numeric(w),
      obs_scale = as.numeric(obs_scale),
      eps_pos = as.numeric(eps_pos),
      eps_mu = as.numeric(eps_shape),
      p1_fix = as.numeric(p1_fix),
      p2_fix = as.numeric(p2_fix),
      m_fix  = as.numeric(m_fix)
    )
  }
  
  obj <- MakeADFun(
    data = data,
    parameters = start_par,
    DLL = DLL,
    silent = TRUE
  )
  
  list(
    obj = obj,
    obs_scale = obs_scale,
    use_idx = use_idx,
    N_prev = N_prev,
    N_obs = N_obs,
    C_prev = C_prev
  )
}

# ============================================================
# Multi-start fit for comparator
# ============================================================
fit_one_year_comparator <- function(method,
                                    N_prev_raw, N_obs_raw, C_prev_raw,
                                    x_mid, x_edges,
                                    p1_fix, p2_fix, m_fix,
                                    n_start = 20,
                                    seed = 1,
                                    control = list(eval.max = 5e4, iter.max = 5e4)) {
  
  set.seed(seed)
  
  best <- NULL
  best_obj <- Inf
  fail_log <- character(0)
  
  for (s in seq_len(n_start)) {
    
    if (method == "normal") {
      start_par <- list(
        log_sigma = rnorm(1, 0, 0.3),
        u1 = log(1.05 * max(x_mid) - max(x_mid)) + rnorm(1, 0, 0.2),
        u2 = log(0.3) + rnorm(1, 0, 0.3),
        u3 = log(1.0) + rnorm(1, 0, 0.3),
        u4 = log(0.1) + rnorm(1, 0, 0.3),
        u5 = rnorm(1, 0, 0.3)
      )
    }
    
    if (method == "gamma") {
      center <- list(log_sigma = 0, z_Linf = 0, z_K = -0.5, z_CV = -0.5)
      start_par <- list(
        log_sigma = center$log_sigma + rnorm(1, 0, 0.1),
        z_Linf    = center$z_Linf    + rnorm(1, 0, 0.3),
        z_K       = center$z_K       + rnorm(1, 0, 0.3),
        z_CV      = center$z_CV      + rnorm(1, 0, 0.3)
      )
    }
    
    if (method == "lognormal") {
      center <- list(log_sigma = 0, z_Linf = 0, z_K = -0.5, z_CV = -0.5)
      start_par <- list(
        log_sigma = center$log_sigma + rnorm(1, 0, 0.1),
        z_Linf    = center$z_Linf    + rnorm(1, 0, 0.3),
        z_K       = center$z_K       + rnorm(1, 0, 0.3),
        z_CV      = center$z_CV      + rnorm(1, 0, 0.3)
      )
    }
    
    tmp <- try(
      build_comparator_obj(
        method = method,
        N_prev_raw = N_prev_raw,
        N_obs_raw = N_obs_raw,
        C_prev_raw = C_prev_raw,
        x_mid = x_mid,
        x_edges = x_edges,
        p1_fix = p1_fix,
        p2_fix = p2_fix,
        m_fix  = m_fix,
        start_par = start_par
      ),
      silent = TRUE
    )
    
    if (inherits(tmp, "try-error")) {
      fail_log <- c(fail_log, paste("start", s, "- MakeADFun failed"))
      next
    }
    
    f0 <- try(tmp$obj$fn(tmp$obj$par), silent = TRUE)
    if (inherits(f0, "try-error") || !is.finite(f0)) {
      fail_log <- c(fail_log, paste("start", s, "- initial fn not finite"))
      next
    }
    
    opt <- try(
      if (method == "normal") {
        g0 <- try(tmp$obj$gr(tmp$obj$par), silent = TRUE)
        if (inherits(g0, "try-error") || any(!is.finite(g0))) {
          fail_log <- c(fail_log, paste("start", s, "- initial gr not finite"))
          next
        }
        nlminb(
          start = tmp$obj$par,
          objective = tmp$obj$fn,
          gradient = tmp$obj$gr,
          control = control
        )
      } else {
        nlminb(
          start = tmp$obj$par,
          objective = tmp$obj$fn,
          control = control
        )
      },
      silent = TRUE
    )
    
    if (inherits(opt, "try-error")) {
      fail_log <- c(fail_log, paste("start", s, "- nlminb failed"))
      next
    }
    
    if (!is.finite(opt$objective)) {
      fail_log <- c(fail_log, paste("start", s, "- objective not finite"))
      next
    }
    
    if (opt$objective < best_obj) {
      best <- list(obj = tmp$obj, opt = opt)
      best_obj <- opt$objective
    }
  }
  
  if (is.null(best)) {
    cat("Failure log for method =", method, "\n")
    print(unique(fail_log))
    stop(paste("All", method, "multi-starts failed."))
  }
  
  H <- try(best$obj$he(best$opt$par), silent = TRUE)
  
  rep <- if (inherits(H, "try-error") || any(!is.finite(H))) {
    sdreport(best$obj, par.fixed = best$opt$par, getJointPrecision = FALSE)
  } else {
    sdreport(best$obj,
             par.fixed = best$opt$par,
             hessian.fixed = H,
             getJointPrecision = FALSE)
  }
  
  list(obj = best$obj, opt = best$opt, rep = rep, method = method)
}

# ============================================================
# Extract method-specific reports
# ============================================================
extract_method_report <- function(fit) {
  rr <- fit$obj$report(fit$opt$par)
  
  out <- list(
    G = rr$G,
    N_pred = as.numeric(rr$N_pred),
    N_pred0 = as.numeric(rr$N_pred0),
    s_y = as.numeric(rr$s_y),
    sigma = as.numeric(rr$sigma)
  )
  
  if (fit$method == "normal") {
    out$theta <- c(
      Linf = as.numeric(rr$Linf),
      K = as.numeric(rr$K),
      sLinf = as.numeric(rr$sLinf),
      sK = as.numeric(rr$sK),
      rho = as.numeric(rr$rho)
    )
  }
  
  if (fit$method == "gamma") {
    out$theta <- c(
      Linf = as.numeric(rr$Linf),
      K = as.numeric(rr$K),
      CV = as.numeric(rr$CV)
    )
  }
  
  if (fit$method == "lognormal") {
    out$theta <- c(
      Linf = as.numeric(rr$Linf),
      K = as.numeric(rr$K),
      CV = as.numeric(rr$CV)
    )
  }
  
  out
}

# ============================================================
# Run all years, all methods
# ============================================================
methods <- c("bspline", "normal", "gamma", "lognormal")

results <- list(
  bspline   = vector("list", length(years)),
  normal    = vector("list", length(years)),
  gamma     = vector("list", length(years)),
  lognormal = vector("list", length(years))
)

names(results$bspline)   <- as.character(years)
names(results$normal)    <- as.character(years)
names(results$gamma)     <- as.character(years)
names(results$lognormal) <- as.character(years)

Err <- list(
  bspline   = rep(NA_real_, length(years)),
  normal    = rep(NA_real_, length(years)),
  gamma     = rep(NA_real_, length(years)),
  lognormal = rep(NA_real_, length(years))
)

RunTime <- list(
  bspline   = rep(NA_real_, length(years)),
  normal    = rep(NA_real_, length(years)),
  gamma     = rep(NA_real_, length(years)),
  lognormal = rep(NA_real_, length(years))
)

Nhat <- list(
  bspline   = matrix(NA_real_, nrow = length(x_mid), ncol = length(years) + 1),
  normal    = matrix(NA_real_, nrow = length(x_mid), ncol = length(years) + 1),
  gamma     = matrix(NA_real_, nrow = length(x_mid), ncol = length(years) + 1),
  lognormal = matrix(NA_real_, nrow = length(x_mid), ncol = length(years) + 1)
)

initial_obs <- get_age_year(data_N32, 3, year_start)
sc0 <- max(initial_obs)
if (sc0 <= 0) sc0 <- 1
for (m in methods) Nhat[[m]][, 1] <- initial_obs / sc0

# ============================================================
# Optional test for gamma
# ============================================================
ii <- 1
y <- years[ii]

N_prev_raw <- get_age_year(data_N32, 2, y)
N_obs_raw  <- get_age_year(data_N32, 3, y + 1)
C_prev_raw <- get_age_year(data_C32, 2, y + 1)

fit_gamma_test <- fit_one_year_comparator(
  method = "gamma",
  N_prev_raw = N_prev_raw,
  N_obs_raw = N_obs_raw,
  C_prev_raw = C_prev_raw,
  x_mid = x_mid,
  x_edges = x_edges,
  p1_fix = p1_by_year[as.character(y)],
  p2_fix = p2_by_year[as.character(y)],
  m_fix  = m_by_year[as.character(y)],
  n_start = 5,
  seed = 4001
)

# ============================================================
# Main loop
# ============================================================
for (ii in seq_along(years)) {
  y <- years[ii]
  
  N_prev_raw <- get_age_year(data_N32, 2, y)
  N_obs_raw  <- get_age_year(data_N32, 3, y + 1)
  C_prev_raw <- get_age_year(data_C32, 2, y + 1)
  
  p1_fix <- p1_by_year[as.character(y)]
  p2_fix <- p2_by_year[as.character(y)]
  m_fix  <- m_by_year[as.character(y)]
  
  obs_sc <- max(N_obs_raw)
  if (obs_sc <= 0) obs_sc <- 1
  N_obs_norm <- N_obs_raw / obs_sc
  
  # ---------------- B-spline ----------------
  t0 <- proc.time()[3]
  results$bspline[[ii]] <- fit_one_year_bspline(
    N_prev_raw = N_prev_raw,
    N_obs_raw  = N_obs_raw,
    C_prev_raw = C_prev_raw,
    x_mid = x_mid,
    I1 = I1,
    I2 = I2,
    p1_fix = p1_fix,
    p2_fix = p2_fix,
    m_fix  = m_fix,
    year_y = y,
    obs_threshold = 1e-12,
    n_start = 20,
    seed = 2000 + ii,
    control = list(eval.max = 5e4, iter.max = 5e4)
  )
  RunTime$bspline[ii] <- proc.time()[3] - t0
  
  rr_bs <- results$bspline[[ii]]$obj$report(results$bspline[[ii]]$opt$par)
  Nhat$bspline[, ii + 1] <- pmax(as.numeric(rr_bs$N_pred), 0)
  Err$bspline[ii] <- norm(Nhat$bspline[, ii + 1] - N_obs_norm, type = "2")
  
  if (y == 2014) {
    cat("B-spline special fix applied for 2014 -> 2015:",
        "lambda_smooth =", results$bspline[[ii]]$lambda_smooth,
        ", support_mask =", results$bspline[[ii]]$use_support_mask, "\n")
  }
  
  # ---------------- Normal ----------------
  t0 <- proc.time()[3]
  results$normal[[ii]] <- fit_one_year_comparator(
    method = "normal",
    N_prev_raw = N_prev_raw,
    N_obs_raw = N_obs_raw,
    C_prev_raw = C_prev_raw,
    x_mid = x_mid,
    x_edges = x_edges,
    p1_fix = p1_fix,
    p2_fix = p2_fix,
    m_fix  = m_fix,
    n_start = 20,
    seed = 3000 + ii
  )
  RunTime$normal[ii] <- proc.time()[3] - t0
  
  rr_n <- extract_method_report(results$normal[[ii]])
  Nhat$normal[, ii + 1] <- pmax(rr_n$N_pred, 0)
  Err$normal[ii] <- norm(Nhat$normal[, ii + 1] - N_obs_norm, type = "2")
  
  # ---------------- Gamma ----------------
  t0 <- proc.time()[3]
  results$gamma[[ii]] <- fit_one_year_comparator(
    method = "gamma",
    N_prev_raw = N_prev_raw,
    N_obs_raw = N_obs_raw,
    C_prev_raw = C_prev_raw,
    x_mid = x_mid,
    x_edges = x_edges,
    p1_fix = p1_fix,
    p2_fix = p2_fix,
    m_fix  = m_fix,
    n_start = 20,
    seed = 4000 + ii
  )
  RunTime$gamma[ii] <- proc.time()[3] - t0
  
  rr_g <- extract_method_report(results$gamma[[ii]])
  Nhat$gamma[, ii + 1] <- pmax(rr_g$N_pred, 0)
  Err$gamma[ii] <- norm(Nhat$gamma[, ii + 1] - N_obs_norm, type = "2")
  
  # ---------------- Log-normal ----------------
  t0 <- proc.time()[3]
  results$lognormal[[ii]] <- fit_one_year_comparator(
    method = "lognormal",
    N_prev_raw = N_prev_raw,
    N_obs_raw = N_obs_raw,
    C_prev_raw = C_prev_raw,
    x_mid = x_mid,
    x_edges = x_edges,
    p1_fix = p1_fix,
    p2_fix = p2_fix,
    m_fix  = m_fix,
    n_start = 20,
    seed = 5000 + ii
  )
  RunTime$lognormal[ii] <- proc.time()[3] - t0
  
  rr_ln <- extract_method_report(results$lognormal[[ii]])
  Nhat$lognormal[, ii + 1] <- pmax(rr_ln$N_pred, 0)
  Err$lognormal[ii] <- norm(Nhat$lognormal[, ii + 1] - N_obs_norm, type = "2")
  
  cat("Finished year", y, "\n")
}

# ============================================================
# Summary
# ============================================================
summary_table <- data.frame(
  year = years + 1,
  bspline = Err$bspline,
  normal = Err$normal,
  gamma = Err$gamma,
  lognormal = Err$lognormal
)

print(summary_table)

cat("\nMedian projection error:\n")
cat("B-spline   :", median(Err$bspline, na.rm = TRUE), "\n")
cat("Normal     :", median(Err$normal, na.rm = TRUE), "\n")
cat("Gamma      :", median(Err$gamma, na.rm = TRUE), "\n")
cat("Log-normal :", median(Err$lognormal, na.rm = TRUE), "\n")

cat("\nMean projection error:\n")
cat("B-spline   :", mean(Err$bspline, na.rm = TRUE), "\n")
cat("Normal     :", mean(Err$normal, na.rm = TRUE), "\n")
cat("Gamma      :", mean(Err$gamma, na.rm = TRUE), "\n")
cat("Log-normal :", mean(Err$lognormal, na.rm = TRUE), "\n")

cat("\nAverage runtime (sec):\n")
cat("B-spline   :", mean(RunTime$bspline, na.rm = TRUE), "\n")
cat("Normal     :", mean(RunTime$normal, na.rm = TRUE), "\n")
cat("Gamma      :", mean(RunTime$gamma, na.rm = TRUE), "\n")
cat("Log-normal :", mean(RunTime$lognormal, na.rm = TRUE), "\n")

# ============================================================
# Observed matrix for plotting
# ============================================================
N_obs_mat <- matrix(NA_real_, nrow = length(x_mid), ncol = length(years) + 1)
N_obs_mat[, 1] <- initial_obs / sc0

for (ii in seq_along(years)) {
  y <- years[ii]
  obs_raw <- get_age_year(data_N32, 3, y + 1)
  sc <- max(obs_raw)
  if (sc <= 0) sc <- 1
  N_obs_mat[, ii + 1] <- obs_raw / sc
}


##---------------------------------
#create the table
##---------------------------------

# Create output folder
dir.create("tables_dir", showWarnings = FALSE)

# Copy and format the table
tab <- summary_table

# Round to 4 decimals and keep trailing zeros
tab$bspline   <- format(round(tab$bspline,   4), nsmall = 4)
tab$normal    <- format(round(tab$normal,    4), nsmall = 4)
tab$gamma     <- format(round(tab$gamma,     4), nsmall = 4)
tab$lognormal <- format(round(tab$lognormal, 4), nsmall = 4)

# Output file
out_file <- file.path("tables_dir", "gtm_errors_table.tex")

# Build LaTeX code
latex_lines <- c(
  "\\begin{table}",
  "\\centering",
  "\\footnotesize",
  "\\caption{Per-year projection error ($\\bigl\\| \\, \\widehat{\\mathbf N}_{y,a}(\\theta) - \\mathbf N_{y,a} \\, \\bigr\\|_2$) for each GTM. $\\mathbf N_{y,a}$ considered in normalized scale.}",
  "\\label{tab:gtm-errors}",
  "\\begin{tabular}{r | r r r r}",
  "\\toprule",
  "Year & B-spline & Normal & Gamma & Log-Normal \\\\",
  "\\midrule"
)

# Add rows automatically
for (i in seq_len(nrow(tab))) {
  latex_lines <- c(
    latex_lines,
    paste0(
      tab$year[i], " & ",
      tab$bspline[i], " & ",
      tab$normal[i], " & ",
      tab$gamma[i], " & ",
      tab$lognormal[i], " \\\\"
    )
  )
}

# Close table
latex_lines <- c(
  latex_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

# Write to file
writeLines(latex_lines, out_file)


#

cat("LaTeX table saved to:\n", out_file, "\n")


##-------------------------------------------
#create the projecttion error comparison plot
##-------------------------------------------

out_file <- file.path("results_plots", "projection_error_comparison.png")

png(out_file, width = 1250, height = 650, res = 140)

par(mar = c(5, 5, 4, 9) + 0.1)

plot(summary_table$year, summary_table$bspline,
     type = "n",
     xlab = "Year",
     ylab = "Projection error",
     main = "Projection error by year",
     xlim = c(1985, 2020),
     ylim = c(0, 0.9),
     xaxt = "n",
     yaxt = "n")

axis(1, at = seq(1985, 2020, by = 5))
axis(2, las = 1)

abline(v = seq(1985, 2020, by = 5), col = "gray85", lty = 3, lwd = 0.7)
abline(h = seq(0, 0.8, by = 0.2), col = "gray85", lty = 3, lwd = 0.7)

box()

lines(summary_table$year, summary_table$bspline,
      type = "o", pch = 1, lty = 1, lwd = 1.7, cex = 0.75, col = "deepskyblue4")

lines(summary_table$year, summary_table$normal,
      type = "o", pch = 0, lty = 2, lwd = 1.7, cex = 0.75, col = "tomato2")

lines(summary_table$year, summary_table$gamma,
      type = "o", pch = 2, lty = 3, lwd = 1.7, cex = 0.75, col = "goldenrod2")

lines(summary_table$year, summary_table$lognormal,
      type = "o", pch = 5, lty = 4, lwd = 1.7, cex = 0.75, col = "purple3")

legend(x = 2015.8, y = 0.87,
       legend = c("B-spline", "Normal", "Gamma", "Log-Normal"),
       col = c("deepskyblue4", "tomato2", "goldenrod2", "purple3"),
       pch = c(1, 0, 2, 5),
       lty = c(1, 2, 3, 4),
       lwd = 1.7,
       pt.cex = 0.6,
       cex = 0.65,
       bty = "n",
       xjust = 0)

dev.off()

##---------------------------------
#create the year comparison plots
##---------------------------------

# Example existing plotting hooks
# compare_and_plot_grid(N_obs_mat, Nhat$bspline, Nhat$normal, Nhat$gamma, Nhat$lognormal)
# plot_error_timeseries(Err$bspline, Err$normal, Err$gamma, Err$lognormal)
# Path for saving


library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ============================================================
# Helper: build long data for plotting
# Assumes:
#   N_obs_mat           : n x (nYears+1)
#   Nhat$bspline        : n x (nYears+1)
#   Nhat$normal         : n x (nYears+1)
#   Nhat$gamma          : n x (nYears+1)
#   Nhat$lognormal      : n x (nYears+1)
#   x_mid               : length n
# First column is the initial year and is skipped.
# ============================================================
build_comparison_plot_data <- function(N_obs_mat,
                                       Nhat_bs,
                                       Nhat_normal,
                                       Nhat_gamma,
                                       Nhat_lognormal,
                                       x_mid,
                                       year_start = 1988,
                                       clip_negative = TRUE) {
  stopifnot(
    is.matrix(N_obs_mat),
    is.matrix(Nhat_bs),
    is.matrix(Nhat_normal),
    is.matrix(Nhat_gamma),
    is.matrix(Nhat_lognormal),
    nrow(N_obs_mat) == length(x_mid),
    all(dim(N_obs_mat) == dim(Nhat_bs)),
    all(dim(N_obs_mat) == dim(Nhat_normal)),
    all(dim(N_obs_mat) == dim(Nhat_gamma)),
    all(dim(N_obs_mat) == dim(Nhat_lognormal))
  )
  
  clip0 <- function(z) if (clip_negative) pmax(z, 0) else z
  
  yrs <- (year_start + 1):(year_start + ncol(N_obs_mat) - 1)
  cols_to_use <- 2:ncol(N_obs_mat)
  
  df_obs <- do.call(
    rbind,
    lapply(seq_along(cols_to_use), function(k) {
      cc <- cols_to_use[k]
      data.frame(
        year = yrs[k],
        x = x_mid,
        value = clip0(N_obs_mat[, cc]),
        method = "Observed"
      )
    })
  )
  
  df_bs <- do.call(
    rbind,
    lapply(seq_along(cols_to_use), function(k) {
      cc <- cols_to_use[k]
      data.frame(
        year = yrs[k],
        x = x_mid,
        value = clip0(Nhat_bs[, cc]),
        method = "B-spline"
      )
    })
  )
  
  df_n <- do.call(
    rbind,
    lapply(seq_along(cols_to_use), function(k) {
      cc <- cols_to_use[k]
      data.frame(
        year = yrs[k],
        x = x_mid,
        value = clip0(Nhat_normal[, cc]),
        method = "Normal"
      )
    })
  )
  
  df_g <- do.call(
    rbind,
    lapply(seq_along(cols_to_use), function(k) {
      cc <- cols_to_use[k]
      data.frame(
        year = yrs[k],
        x = x_mid,
        value = clip0(Nhat_gamma[, cc]),
        method = "Gamma"
      )
    })
  )
  
  df_ln <- do.call(
    rbind,
    lapply(seq_along(cols_to_use), function(k) {
      cc <- cols_to_use[k]
      data.frame(
        year = yrs[k],
        x = x_mid,
        value = clip0(Nhat_lognormal[, cc]),
        method = "Log-Normal"
      )
    })
  )
  
  bind_rows(df_obs, df_bs, df_n, df_g, df_ln) %>%
    mutate(
      method = factor(
        method,
        levels = c("Observed", "B-spline", "Normal", "Gamma", "Log-Normal")
      ),
      year = factor(year, levels = yrs)
    )
}

# ============================================================
# Helper: one faceted page
# ============================================================
make_comparison_page <- function(df_page,
                                 step = 0.5,
                                 ncol = 3,
                                 base_size = 11,
                                 show_legend = TRUE) {
  df_obs  <- df_page %>% filter(method == "Observed")
  df_line <- df_page %>% filter(method != "Observed")
  
  ggplot() +
    geom_col(
      data = df_obs,
      aes(x = x, y = value, fill = method),
      width = step * 0.85,
      alpha = 0.55,
      color = NA
    ) +
    geom_line(
      data = df_line,
      aes(x = x, y = value, color = method, linetype = method),
      linewidth = 0.6
    ) +
    facet_wrap(~ year, ncol = ncol, scales = "fixed") +
    scale_fill_manual(
      values = c("Observed" = "#9BB7FF"),
      drop = FALSE
    ) +
    scale_color_manual(
      values = c(
        "B-spline"   = "black",
        "Normal"     = "#D55E00",
        "Gamma"      = "#0072B2",
        "Log-Normal" = "#009E73"
      ),
      drop = FALSE
    ) +
    scale_linetype_manual(
      values = c(
        "B-spline"   = "solid",
        "Normal"     = "dashed",
        "Gamma"      = "dotted",
        "Log-Normal" = "dotdash"
      ),
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = pretty(df_page$x, n = 5)
    ) +
    labs(
      x = "Length class midpoint",
      y = "Normalized abundance"
    ) +
    theme_bw(base_size = base_size) +
    theme(
      legend.position = if (show_legend) "top" else "none",
      legend.title = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.8, "lines"),
      axis.title.x = element_text(size = base_size),
      axis.title.y = element_text(size = base_size),
      strip.text = element_text(face = "bold", size = base_size - 1)
    ) +
    guides(
      fill = guide_legend(order = 1),
      color = guide_legend(order = 2),
      linetype = guide_legend(order = 2)
    )
}

# ============================================================
# Main function:
# Creates 2 or 3 figures automatically
# ============================================================
plot_comparison_pages <- function(N_obs_mat,
                                  Nhat_bs,
                                  Nhat_normal,
                                  Nhat_gamma,
                                  Nhat_lognormal,
                                  x_mid,
                                  year_start = 1988,
                                  step = 0.5,
                                  n_pages = 3,
                                  ncol = 3,
                                  base_size = 11,
                                  save_prefix = NULL,
                                  width = 11,
                                  height = 8.5,
                                  dpi = 300) {
  df_all <- build_comparison_plot_data(
    N_obs_mat = N_obs_mat,
    Nhat_bs = Nhat_bs,
    Nhat_normal = Nhat_normal,
    Nhat_gamma = Nhat_gamma,
    Nhat_lognormal = Nhat_lognormal,
    x_mid = x_mid,
    year_start = year_start,
    clip_negative = TRUE
  )
  
  years_num <- as.integer(as.character(unique(df_all$year)))
  years_num <- sort(years_num)
  
  # split years into 2 or 3 groups
  split_idx <- split(years_num, cut(seq_along(years_num), breaks = n_pages, labels = FALSE))
  
  plot_list <- vector("list", length(split_idx))
  
  for (i in seq_along(split_idx)) {
    yrs_i <- split_idx[[i]]
    df_i <- df_all %>% filter(as.integer(as.character(year)) %in% yrs_i)
    
    p_i <- make_comparison_page(
      df_page = df_i,
      step = step,
      ncol = ncol,
      base_size = base_size,
      show_legend = TRUE
    ) +
      plot_annotation(
        title = paste0("Observed vs projected length distributions (Years ",
                       min(yrs_i), "–", max(yrs_i), ")")
      )
    
    plot_list[[i]] <- p_i
    
    if (!is.null(save_prefix)) {
      ggsave(
        filename = paste0(save_prefix, "_part", i, ".png"),
        plot = p_i,
        width = width,
        height = height,
        dpi = dpi
      )
    }
  }
  
  plot_list
}


plots_cmp <- plot_comparison_pages(
  N_obs_mat = N_obs_mat,
  Nhat_bs = Nhat$bspline,
  Nhat_normal = Nhat$normal,
  Nhat_gamma = Nhat$gamma,
  Nhat_lognormal = Nhat$lognormal,
  x_mid = x_mid,
  year_start = 1988,
  step = 0.5,
  n_pages = 3,
  save_prefix = file.path("results_plots", "appendix_projection_comparison")
)

plots_cmp[[1]]
plots_cmp[[2]]
plots_cmp[[3]]

pdf(file.path("results_plots","appendix_projection_comparison.pdf"),
    width = 11, height = 8.5)

for(p in plots_cmp) print(p)

dev.off()


###---------------------------------
## Same comparison plots in two pages
###---------------------------------



# ============================================================
# Build long data
# ============================================================
build_comparison_plot_data <- function(N_obs_mat,
                                       Nhat_bs,
                                       Nhat_normal,
                                       Nhat_gamma,
                                       Nhat_lognormal,
                                       x_mid,
                                       year_start = 1988,
                                       clip_negative = TRUE) {
  stopifnot(
    is.matrix(N_obs_mat),
    is.matrix(Nhat_bs),
    is.matrix(Nhat_normal),
    is.matrix(Nhat_gamma),
    is.matrix(Nhat_lognormal),
    nrow(N_obs_mat) == length(x_mid),
    all(dim(N_obs_mat) == dim(Nhat_bs)),
    all(dim(N_obs_mat) == dim(Nhat_normal)),
    all(dim(N_obs_mat) == dim(Nhat_gamma)),
    all(dim(N_obs_mat) == dim(Nhat_lognormal))
  )
  
  clip0 <- function(z) if (clip_negative) pmax(z, 0) else z
  
  yrs <- (year_start + 1):(year_start + ncol(N_obs_mat) - 1)
  cols_to_use <- 2:ncol(N_obs_mat)
  
  make_df <- function(mat, method_name) {
    do.call(
      rbind,
      lapply(seq_along(cols_to_use), function(k) {
        cc <- cols_to_use[k]
        data.frame(
          year = yrs[k],
          x = x_mid,
          value = clip0(mat[, cc]),
          method = method_name
        )
      })
    )
  }
  
  bind_rows(
    make_df(N_obs_mat, "Observed"),
    make_df(Nhat_bs, "B-spline"),
    make_df(Nhat_normal, "Normal"),
    make_df(Nhat_gamma, "Gamma"),
    make_df(Nhat_lognormal, "Log-Normal")
  ) %>%
    mutate(
      method = factor(
        method,
        levels = c("Observed", "B-spline", "Normal", "Gamma", "Log-Normal")
      )
    )
}

# ============================================================
# One yearly panel
# ============================================================
make_year_plot <- function(df_year,
                           step = 0.5,
                           base_size = 8.5,
                           show_legend = FALSE) {
  df_obs  <- df_year %>% filter(method == "Observed")
  df_line <- df_year %>% filter(method != "Observed")
  yy <- unique(df_year$year)
  
  ggplot() +
    geom_col(
      data = df_obs,
      aes(x = x, y = value, fill = method),
      width = step * 0.85,
      alpha = 0.50,
      color = NA
    ) +
    geom_line(
      data = df_line,
      aes(x = x, y = value, color = method, linetype = method),
      linewidth = 0.45
    ) +
    scale_fill_manual(
      values = c("Observed" = "#9BB7FF"),
      drop = FALSE
    ) +
    scale_color_manual(
      values = c(
        "B-spline"   = "black",
        "Normal"     = "#D55E00",
        "Gamma"      = "#0072B2",
        "Log-Normal" = "#009E73"
      ),
      drop = FALSE
    ) +
    scale_linetype_manual(
      values = c(
        "B-spline"   = "solid",
        "Normal"     = "dashed",
        "Gamma"      = "dotted",
        "Log-Normal" = "dotdash"
      ),
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = c(6, 10, 14, 18),
      limits = range(df_year$x)
    ) +
    labs(
      title = as.character(yy),
      x = "Length midpoint",
      y = "Normalized abundance"
    ) +
    theme_bw(base_size = base_size) +
    theme(
      legend.position = if (show_legend) "top" else "none",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1),
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 0.2),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = base_size - 0.2),
      axis.title.y = element_text(size = base_size - 0.2),
      axis.text.x = element_text(size = base_size - 1),
      axis.text.y = element_text(size = base_size - 1),
      plot.margin = margin(3, 3, 3, 3)
    ) +
    guides(
      fill = guide_legend(order = 1),
      color = guide_legend(order = 2),
      linetype = guide_legend(order = 2)
    )
}

# ============================================================
# Blank spacer with no panel border
# ============================================================
make_empty_plot <- function() {
  ggplot() +
    theme_void()
}

# ============================================================
# Fixed 2-page paper layout:
# page 1 = 18 panels
# page 2 = 11 panels + 7 empty spaces
# so panel size remains identical
# ============================================================
plot_comparison_pages_2paper_fixed <- function(N_obs_mat,
                                               Nhat_bs,
                                               Nhat_normal,
                                               Nhat_gamma,
                                               Nhat_lognormal,
                                               x_mid,
                                               year_start = 1988,
                                               step = 0.5,
                                               base_size = 8.5,
                                               ncol = 3,
                                               save_prefix = NULL,
                                               width = 8.27,
                                               height = 11.69,
                                               dpi = 300) {
  df_all <- build_comparison_plot_data(
    N_obs_mat = N_obs_mat,
    Nhat_bs = Nhat_bs,
    Nhat_normal = Nhat_normal,
    Nhat_gamma = Nhat_gamma,
    Nhat_lognormal = Nhat_lognormal,
    x_mid = x_mid,
    year_start = year_start,
    clip_negative = TRUE
  )
  
  years_num <- sort(unique(df_all$year))
  stopifnot(length(years_num) == 29)
  
  yrs_page1 <- years_num[1:18]   # 1989-2006
  yrs_page2 <- years_num[19:29]  # 2007-2017
  
  # build one plot per year
  year_plot_list <- lapply(years_num, function(yy) {
    make_year_plot(
      df_year = df_all %>% filter(year == yy),
      step = step,
      base_size = base_size,
      show_legend = FALSE
    )
  })
  names(year_plot_list) <- as.character(years_num)
  
  # one legend donor plot
  legend_plot <- make_year_plot(
    df_year = df_all %>% filter(year == years_num[1]),
    step = step,
    base_size = base_size,
    show_legend = TRUE
  )
  
  # page 1: 18 real plots
  plots_page1 <- year_plot_list[as.character(yrs_page1)]
  
  page1_grid <- wrap_plots(plots_page1, ncol = ncol, guides = "collect") &
    theme(legend.position = "top")
  
  page1 <- legend_plot + plot_layout(guides = "collect") &
    theme(legend.position = "top")
  
  page1 <- (page1_grid) +
    plot_annotation(
      title = paste0("Observed vs projected length distributions (", min(yrs_page1), "–", max(yrs_page1), ")")
    ) &
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      legend.position = "top"
    )
  
  # page 2: 11 real plots + 7 blanks to preserve same panel size
  plots_page2 <- year_plot_list[as.character(yrs_page2)]
  n_blank <- 18 - length(plots_page2)
  blank_list <- replicate(n_blank, make_empty_plot(), simplify = FALSE)
  
  page2_grid <- wrap_plots(c(plots_page2, blank_list), ncol = ncol, guides = "collect") &
    theme(legend.position = "top")
  
  page2 <- (page2_grid) +
    plot_annotation(
      title = paste0("Observed vs projected length distributions (", min(yrs_page2), "–", max(yrs_page2), ")")
    ) &
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      legend.position = "top"
    )
  
  if (!is.null(save_prefix)) {
    ggsave(
      filename = paste0(save_prefix, "_page1.png"),
      plot = page1,
      width = width,
      height = height,
      dpi = dpi,
      units = "in"
    )
    
    ggsave(
      filename = paste0(save_prefix, "_page2.png"),
      plot = page2,
      width = width,
      height = height,
      dpi = dpi,
      units = "in"
    )
    
    pdf(
      file = paste0(save_prefix, ".pdf"),
      width = width,
      height = height,
      onefile = TRUE
    )
    print(page1)
    print(page2)
    dev.off()
  }
  
  list(page1 = page1, page2 = page2)
}

plots_cmp <- plot_comparison_pages_2paper_fixed(
  N_obs_mat = N_obs_mat,
  Nhat_bs = Nhat$bspline,
  Nhat_normal = Nhat$normal,
  Nhat_gamma = Nhat$gamma,
  Nhat_lognormal = Nhat$lognormal,
  x_mid = x_mid,
  year_start = 1988,
  step = 0.5,
  base_size = 8.5,
  ncol = 3,
  save_prefix = file.path("results_plots", "appendix_projection_comparison"),
  width = 8.27,
  height = 11.69,
  dpi = 300
)

plots_cmp$page1
plots_cmp$page2




##########################################################
##########################################################
## Out-of-sample comparison 
##########################################################
##########################################################

# ============================================================
# Fit all methods year-by-year and store GTMs
# ============================================================

fits_bspline   <- vector("list", length(years))
fits_normal    <- vector("list", length(years))
fits_gamma     <- vector("list", length(years))
fits_lognormal <- vector("list", length(years))

G_bspline   <- vector("list", length(years))
G_normal    <- vector("list", length(years))
G_gamma     <- vector("list", length(years))
G_lognormal <- vector("list", length(years))

names(fits_bspline)   <- names(fits_normal) <- names(fits_gamma) <- names(fits_lognormal) <- as.character(years)
names(G_bspline)      <- names(G_normal)    <- names(G_gamma)    <- names(G_lognormal)    <- as.character(years)

for (ii in seq_along(years)) {
  y <- years[ii]
  cat("Fitting year:", y, "\n")
  
  N_prev_raw <- get_age_year(data_N32, 2, y)
  N_obs_raw  <- get_age_year(data_N32, 3, y + 1)
  C_prev_raw <- get_age_year(data_C32, 2, y + 1)
  
  p1_fix <- p1_by_year[as.character(y)]
  p2_fix <- p2_by_year[as.character(y)]
  m_fix  <- m_by_year[as.character(y)]
  
  # ---------------- B-spline ----------------
  fits_bspline[[ii]] <- fit_one_year_bspline(
    N_prev_raw = N_prev_raw,
    N_obs_raw  = N_obs_raw,
    C_prev_raw = C_prev_raw,
    x_mid = x_mid,
    I1 = I1,
    I2 = I2,
    p1_fix = p1_fix,
    p2_fix = p2_fix,
    m_fix  = m_fix,
    year_y = y,
    n_start = 20,
    seed = 1000 + ii
  )
  G_bspline[[ii]] <- extract_method_report(fits_bspline[[ii]])$G
  
  # ---------------- Normal ----------------
  fits_normal[[ii]] <- fit_one_year_comparator(
    method = "normal",
    N_prev_raw = N_prev_raw,
    N_obs_raw  = N_obs_raw,
    C_prev_raw = C_prev_raw,
    x_mid = x_mid,
    x_edges = x_edges,
    p1_fix = p1_fix,
    p2_fix = p2_fix,
    m_fix  = m_fix,
    n_start = 20,
    seed = 2000 + ii
  )
  G_normal[[ii]] <- extract_method_report(fits_normal[[ii]])$G
  
  # ---------------- Gamma ----------------
  fits_gamma[[ii]] <- fit_one_year_comparator(
    method = "gamma",
    N_prev_raw = N_prev_raw,
    N_obs_raw  = N_obs_raw,
    C_prev_raw = C_prev_raw,
    x_mid = x_mid,
    x_edges = x_edges,
    p1_fix = p1_fix,
    p2_fix = p2_fix,
    m_fix  = m_fix,
    n_start = 20,
    seed = 3000 + ii
  )
  G_gamma[[ii]] <- extract_method_report(fits_gamma[[ii]])$G
  
  # ---------------- Log-normal ----------------
  fits_lognormal[[ii]] <- fit_one_year_comparator(
    method = "lognormal",
    N_prev_raw = N_prev_raw,
    N_obs_raw  = N_obs_raw,
    C_prev_raw = C_prev_raw,
    x_mid = x_mid,
    x_edges = x_edges,
    p1_fix = p1_fix,
    p2_fix = p2_fix,
    m_fix  = m_fix,
    n_start = 20,
    seed = 4000 + ii
  )
  G_lognormal[[ii]] <- extract_method_report(fits_lognormal[[ii]])$G
}

# ============================================================
# Helpers for leave-one-year-out cluster-mean validation
# ============================================================

maturity_vec_R <- function(x_mid, p1, p2) {
  1 / (1 + exp(4 * p1 * (p2 - x_mid)))
}

profile_scale <- function(N_obs, N_pred0, w = NULL) {
  if (is.null(w)) w <- rep(1, length(N_obs))
  idx <- which(w > 0)
  num <- sum(N_obs[idx] * N_pred0[idx])
  den <- sum(N_pred0[idx]^2)
  s_y <- if (den > 1e-12) num / den else 1
  max(s_y, 1e-12)
}

compute_inside_state_norm <- function(year_y) {
  N_prev_raw <- get_age_year(data_N32, 2, year_y)
  N_obs_raw  <- get_age_year(data_N32, 3, year_y + 1)
  C_prev_raw <- get_age_year(data_C32, 2, year_y + 1)
  
  p1 <- unname(p1_by_year[as.character(year_y)])
  p2 <- unname(p2_by_year[as.character(year_y)])
  m  <- unname(m_by_year[as.character(year_y)])
  
  max_prev <- max(N_prev_raw)
  if (!is.finite(max_prev) || max_prev <= 0) max_prev <- 1
  
  max_obs <- max(N_obs_raw)
  if (!is.finite(max_obs) || max_obs <= 0) max_obs <- 1
  
  # normalize exactly as in the estimation code
  N_prev <- N_prev_raw / max_prev
  N_obs  <- N_obs_raw  / max_obs
  C_prev <- C_prev_raw / max_prev
  
  rmat <- maturity_vec_R(x_mid, p1, p2)
  N_imm <- N_prev * (1 - rmat)
  inside <- exp(-m) * N_imm - exp(-0.5 * m) * C_prev
  
  list(
    N_obs = N_obs,
    inside = inside,
    max_prev = max_prev,
    max_obs = max_obs
  )
}

mean_G_cluster_loo <- function(G_list, clusters, holdout_idx) {
  cl <- clusters[holdout_idx]
  train_idx <- which(clusters == cl & seq_along(clusters) != holdout_idx)
  
  if (length(train_idx) == 0) stop("No training years left in this cluster.")
  
  G_arr <- array(NA_real_, dim = c(32, 32, length(train_idx)))
  for (jj in seq_along(train_idx)) {
    G_arr[, , jj] <- G_list[[train_idx[jj]]]
  }
  apply(G_arr, c(1, 2), mean)
}

# ============================================================
# LOO validation using cluster-specific mean GTMs
# ============================================================

validate_method_loo <- function(G_list, method_name, clusters, obs_threshold = 1e-12) {
  out <- data.frame(
    method = method_name,
    holdout_year = years,
    cluster = clusters,
    s_y = NA_real_,
    error_L2 = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (ii in seq_along(years)) {
    y <- years[ii]
    
    G_bar <- mean_G_cluster_loo(G_list, clusters, ii)
    st <- compute_inside_state_norm(y)
    
    N_obs <- st$N_obs
    inside <- st$inside
    
    N_pred0 <- as.numeric(G_bar %*% inside)
    N_pred0[N_pred0 < 0] <- 0
    
    w <- as.numeric(N_obs > obs_threshold)
    if (sum(w) == 0) w <- rep(1, length(N_obs))
    
    s_y <- profile_scale(N_obs, N_pred0, w)
    N_pred <- s_y * N_pred0
    
    err <- sqrt(sum(((N_obs - N_pred) * w)^2))
    
    out$s_y[ii] <- s_y
    out$error_L2[ii] <- err
  }
  
  out
}

# ============================================================
# Use the real cluster labels from k-means
# ============================================================

clusters <- km_out$idx
print(clusters)   # should be only 1 and 2

cv_bspline   <- validate_method_loo(G_bspline,   "B-spline",   clusters)
cv_normal    <- validate_method_loo(G_normal,    "Normal",     clusters)
cv_gamma     <- validate_method_loo(G_gamma,     "Gamma",      clusters)
cv_lognormal <- validate_method_loo(G_lognormal, "Log-normal", clusters)

cv_all <- rbind(cv_bspline, cv_normal, cv_gamma, cv_lognormal)

cat("\n--- Mean out-of-sample error by method ---\n")
print(aggregate(error_L2 ~ method, data = cv_all, mean))

cat("\n--- Median out-of-sample error by method ---\n")
print(aggregate(error_L2 ~ method, data = cv_all, median))

cat("\n--- SD out-of-sample error by method ---\n")
print(aggregate(error_L2 ~ method, data = cv_all, sd))

#######
# plot
library(ggplot2)

p_cv <- ggplot(cv_all, aes(x = holdout_year, y = error_L2, color = method, linetype = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  labs(
    x = "Held-out year",
    y = "Out-of-sample L2 error (normalized scale)"
  ) +
  theme_classic(base_size = 13)

print(p_cv)

ggsave(
  filename = file.path("results_plots", "cv_cluster_mean_gtm_errors.png"),
  plot = p_cv,
  width = 8,
  height = 5,
  dpi = 300
)



################################################
###############################################
plot_bspline_prediction <- function(year_y,
                                    G_bspline,
                                    clusters,
                                    years,
                                    x_mid) {
  
  ii <- which(years == year_y)
  
  # mean GTM of same cluster excluding this year
  G_bar <- mean_G_cluster_loo(G_bspline, clusters, ii)
  
  # build normalized state
  st <- compute_inside_state_norm(year_y)
  
  N_obs <- st$N_obs
  inside <- st$inside
  
  # prediction
  N_pred0 <- as.numeric(G_bar %*% inside)
  N_pred0[N_pred0 < 0] <- 0
  
  s_y <- profile_scale(N_obs, N_pred0)
  N_pred <- s_y * N_pred0
  
  df <- data.frame(
    length = x_mid,
    Observed = N_obs,
    Predicted = N_pred
  )
  
  library(tidyr)
  df_long <- pivot_longer(df, -length,
                          names_to = "type",
                          values_to = "value")
  
  library(ggplot2)
  
  ggplot(df_long, aes(x = length, y = value, color = type)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    theme_classic(base_size = 14) +
    labs(
      title = paste("B-spline cluster-mean prediction vs observed:", year_y + 1),
      x = "Length class midpoint",
      y = "Normalized abundance"
    ) +
    scale_color_manual(values = c("Observed"="black","Predicted"="red"))
}

plot_bspline_prediction(1992, G_bspline, clusters, years, x_mid)

plot_bspline_prediction(1996, G_bspline, clusters, years, x_mid)
