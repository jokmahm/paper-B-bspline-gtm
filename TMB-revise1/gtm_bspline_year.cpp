#include <TMB.hpp>

template<class Type>
Type posfun_local(const Type &x, const Type &eps, Type &pen) {
  if (x >= eps) return x;
  Type y = eps / (Type(2) - x / eps);
  pen += pow(x - eps, 2);
  return y;
}

template<class Type>
vector<Type> maturity_vec(const vector<Type>& x_mid, const Type& p1, const Type& p2) {
  vector<Type> r(x_mid.size());
  for (int i = 0; i < x_mid.size(); i++) {
    r(i) = Type(1) / (Type(1) + exp(Type(4) * p1 * (p2 - x_mid(i))));
  }
  return r;
}

template<class Type>
matrix<Type> build_G(const vector<Type>& q_full,
                     const matrix<Type>& I1,
                     const matrix<Type>& I2,
                     const Type& eps_denom,
                     Type& pen) {
  int K = q_full.size();
  int n = I1.cols();
  
  matrix<Type> G(n, n);
  G.setZero();
  
  // Denominator for each source class j
  vector<Type> S(n);
  for (int j = 0; j < n; j++) {
    Type Sj = Type(0);
    for (int k = 0; k < K; k++) Sj += q_full(k) * I1(k, j);
    S(j) = posfun_local(Sj, eps_denom, pen);
  }
  
  // Numerator for each target class i
  vector<Type> num_i(n);
  for (int i = 0; i < n; i++) {
    Type tmp = Type(0);
    for (int k = 0; k < K; k++) tmp += q_full(k) * I2(k, i);
    num_i(i) = tmp;
  }
  
  // Lower-triangular GTM
  for (int j = 0; j < n; j++) {
    for (int i = j; i < n; i++) {
      G(i, j) = num_i(i) / S(j);
    }
  }
  
  return G;
}

template<class Type>
Type objective_function<Type>::operator() () {
  // --------------------
  // DATA
  // --------------------
  DATA_VECTOR(N_prev);     // length n = 32
  DATA_VECTOR(N_obs);      // length n = 32
  DATA_VECTOR(C_prev);     // length n = 32
  DATA_VECTOR(x_mid);      // length n = 32
  DATA_MATRIX(I1);         // K x n
  DATA_MATRIX(I2);         // K x n
  DATA_SCALAR(eps_pos);
  DATA_SCALAR(eps_denom);
  DATA_VECTOR(w);          // length n, active bins
  DATA_SCALAR(obs_scale);
  
  DATA_SCALAR(p1_fix);
  DATA_SCALAR(p2_fix);
  DATA_SCALAR(m_fix);
  
  DATA_VECTOR(q_mask);         // length K, 0/1
  DATA_SCALAR(lambda_smooth);  // smoothness penalty weight
  
  // --------------------
  // PARAMETERS
  // --------------------
  PARAMETER(log_sigma);
  PARAMETER_VECTOR(q_raw_full);
  
  // --------------------
  // TRANSFORMS
  // --------------------
  Type sigma = obs_scale * exp(log_sigma);
  
  int K = q_raw_full.size();
  vector<Type> q_full(K);
  
  for (int k = 0; k < K; k++) {
    Type q_raw_k = q_raw_full(k);
    
    // Clip exponent to reduce overflow/underflow during Hessian evaluation
    if (q_raw_k > Type(30))  q_raw_k = Type(30);
    if (q_raw_k < Type(-30)) q_raw_k = Type(-30);
    
    Type qpos = exp(q_raw_k);
    q_full(k) = q_mask(k) * qpos;
  }
  
  // Normalize to remove scale non-identifiability with s_y
  Type q_sum = Type(0);
  for (int k = 0; k < K; k++) q_sum += q_full(k);
  
  if (q_sum > Type(0)) {
    for (int k = 0; k < K; k++) q_full(k) /= q_sum;
  }
  
  // --------------------
  // MODEL
  // --------------------
  vector<Type> rmat = maturity_vec(x_mid, p1_fix, p2_fix);
  
  vector<Type> Surv(N_prev.size());
  for (int i = 0; i < N_prev.size(); i++) {
    Surv(i) = N_prev(i) * (Type(1) - rmat(i));
  }
  
  vector<Type> M(N_prev.size());
  for (int i = 0; i < N_prev.size(); i++) {
    M(i) = m_fix;
  }
  
  Type pen = Type(0);
  
  matrix<Type> G = build_G(q_full, I1, I2, eps_denom, pen);
  
  vector<Type> inside(N_prev.size());
  for (int i = 0; i < N_prev.size(); i++) {
    Type e6  = exp(Type(-0.5) * M(i));
    Type e12 = exp(Type(-1.0) * M(i));
    inside(i) = e12 * Surv(i) - e6 * C_prev(i);
  }
  
  vector<Type> N_pred0 = G * inside;
  
  for (int i = 0; i < N_pred0.size(); i++) {
    N_pred0(i) = posfun_local(N_pred0(i), eps_pos, pen);
  }
  
  // --------------------
  // PROFILED SCALE s_y
  // --------------------
  Type num = Type(0);
  Type den = Type(0);
  
  for (int i = 0; i < N_obs.size(); i++) {
    if (w(i) > Type(0)) {
      num += N_pred0(i) * N_obs(i);
      den += N_pred0(i) * N_pred0(i);
    }
  }
  
  Type s_y = Type(1.0);
  if (den > Type(1e-12)) s_y = num / den;
  if (s_y < Type(1e-12)) s_y = Type(1e-12);
  
  vector<Type> N_pred = s_y * N_pred0;
  
  // --------------------
  // LIKELIHOOD
  // --------------------
  Type nll = Type(0);
  for (int i = 0; i < N_obs.size(); i++) {
    if (w(i) > Type(0)) {
      nll -= dnorm(N_obs(i), N_pred(i), sigma, true);
    }
  }
  
  // positivity penalty from posfun
  nll += Type(1e-4) * pen;
  
  // optional smoothness penalty on q_raw_full
  Type pen_smooth = Type(0);
  if (lambda_smooth > Type(0) && K >= 3) {
    for (int k = 1; k < K - 1; k++) {
      Type d2 = q_raw_full(k + 1) - Type(2.0) * q_raw_full(k) + q_raw_full(k - 1);
      pen_smooth += d2 * d2;
    }
    nll += lambda_smooth * pen_smooth;
  }
  
  // --------------------
  // REPORT
  // --------------------
  REPORT(q_full);
  REPORT(s_y);
  REPORT(sigma);
  
  ADREPORT(q_full);
  ADREPORT(s_y);
  ADREPORT(sigma);
  
  REPORT(G);
  REPORT(N_pred0);
  REPORT(N_pred);
  REPORT(inside);
  REPORT(rmat);
  REPORT(M);
  REPORT(pen);
  REPORT(pen_smooth);
  
  return nll;
}
