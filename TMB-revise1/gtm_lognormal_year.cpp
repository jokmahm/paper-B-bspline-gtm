#include <TMB.hpp>

template<class Type>
Type posfun_local(const Type &x, const Type &eps, Type &pen) {
  if (x >= eps) return x;
  Type y = eps / (Type(2) - x / eps);
  pen += (x - eps) * (x - eps);
  return y;
}

template<class Type>
Type invlogit_local(const Type& x) {
  return Type(1) / (Type(1) + exp(-x));
}

template<class Type>
Type map_interval(const Type& z, const Type& lo, const Type& hi) {
  return lo + (hi - lo) * invlogit_local(z);
}

template<class Type>
vector<Type> maturity_vec(const vector<Type>& x_mid,
                          const Type& p1,
                          const Type& p2) {
  int n = x_mid.size();
  vector<Type> r(n);
  for (int i = 0; i < n; i++) {
    r(i) = Type(1) / (Type(1) + exp(Type(4) * p1 * (p2 - x_mid(i))));
  }
  return r;
}

template<class Type>
vector<Type> mortality_vec_constant(const int n, const Type& m_fix) {
  vector<Type> M(n);
  for (int i = 0; i < n; i++) M(i) = m_fix;
  return M;
}

template<class Type>
Type plnorm_local(const Type& x, const Type& meanlog, const Type& sdlog) {
  if (x <= Type(0)) return Type(0);
  return pnorm((log(x) - meanlog) / sdlog);
}

template<class Type>
matrix<Type> build_G_lognormal(const vector<Type>& x_mid,
                               const vector<Type>& x_edges,
                               const Type& Linf,
                               const Type& K,
                               const Type& CV,
                               const Type& eps_mu,
                               Type& pen) {
  int n = x_mid.size();
  matrix<Type> G(n, n);
  G.setZero();
  
  Type eK = exp(-K);
  Type one_meK = Type(1) - eK;
  Type max_edge = x_edges(n);
  
  Type s2 = log(Type(1) + CV * CV);
  Type s  = sqrt(s2);
  
  for (int j = 0; j < n; j++) {
    Type lj = x_mid(j);
    
    Type mu = (Linf - lj) * one_meK;
    mu = posfun_local(mu, eps_mu, pen);
    
    Type mlog = log(mu) - Type(0.5) * s2;
    
    Type a_tr = Type(0);
    Type b_tr = max_edge - lj;
    
    Type Z = plnorm_local(b_tr, mlog, s) - plnorm_local(a_tr, mlog, s);
    
    if (Z <= Type(1e-14)) {
      G(j, j) = Type(1);
      continue;
    }
    
    for (int i = j; i < n; i++) {
      Type lo = x_edges(i)   - lj;
      Type hi = x_edges(i+1) - lj;
      
      if (lo < a_tr) lo = a_tr;
      if (hi > b_tr) hi = b_tr;
      
      Type pij = Type(0);
      if (hi > lo) {
        pij = (plnorm_local(hi, mlog, s) - plnorm_local(lo, mlog, s)) / Z;
      }
      G(i, j) = pij;
    }
    
    Type cs = Type(0);
    for (int i = 0; i < n; i++) cs += G(i, j);
    
    if (cs > Type(1e-14)) {
      for (int i = 0; i < n; i++) G(i, j) /= cs;
    } else {
      G(j, j) = Type(1);
    }
  }
  
  return G;
}

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_VECTOR(N_prev);
  DATA_VECTOR(N_obs);
  DATA_VECTOR(C_prev);
  DATA_VECTOR(x_mid);
  DATA_VECTOR(x_edges);
  DATA_VECTOR(w);
  DATA_SCALAR(obs_scale);
  DATA_SCALAR(eps_pos);
  DATA_SCALAR(eps_mu);
  
  DATA_SCALAR(p1_fix);
  DATA_SCALAR(p2_fix);
  DATA_SCALAR(m_fix);
  
  PARAMETER(log_sigma);
  PARAMETER(z_Linf);
  PARAMETER(z_K);
  PARAMETER(z_CV);
  
  int n = x_mid.size();
  Type pen = Type(0);
  
  Type sigma = obs_scale * exp(log_sigma);
  
  Type xmax = x_mid.maxCoeff();
  
  Type Linf = map_interval(z_Linf, xmax + Type(1e-6), Type(2.0) * xmax);
  Type K    = map_interval(z_K,    Type(0.01),       Type(2.0));
  Type CV   = map_interval(z_CV,   Type(0.02),       Type(1.0));
  
  vector<Type> rmat = maturity_vec(x_mid, p1_fix, p2_fix);
  vector<Type> M = mortality_vec_constant(n, m_fix);
  
  vector<Type> Surv(n);
  for (int i = 0; i < n; i++) {
    Surv(i) = N_prev(i) * (Type(1) - rmat(i));
  }
  
  matrix<Type> G = build_G_lognormal(x_mid, x_edges, Linf, K, CV, eps_mu, pen);
  
  vector<Type> inside(n);
  for (int i = 0; i < n; i++) {
    Type em = exp(-M(i));
    Type eh = exp(-Type(0.5) * M(i));
    inside(i) = em * Surv(i) - eh * C_prev(i);
  }
  
  vector<Type> N_pred0 = G * inside;
  for (int i = 0; i < n; i++) {
    N_pred0(i) = posfun_local(N_pred0(i), eps_pos, pen);
  }
  
  Type num = Type(0);
  Type den = Type(0);
  for (int i = 0; i < n; i++) {
    if (w(i) > Type(0)) {
      num += N_pred0(i) * N_obs(i);
      den += N_pred0(i) * N_pred0(i);
    }
  }
  
  Type s_y = Type(1);
  if (den > Type(1e-12)) s_y = num / den;
  if (s_y < Type(1e-12)) s_y = Type(1e-12);
  
  vector<Type> N_pred = s_y * N_pred0;
  
  Type nll = Type(0);
  for (int i = 0; i < n; i++) {
    if (w(i) > Type(0)) {
      nll -= dnorm(N_obs(i), N_pred(i), sigma, true);
    }
  }
  
  nll += Type(1e-4) * pen;
  
  REPORT(G);
  REPORT(N_pred0);
  REPORT(N_pred);
  REPORT(inside);
  REPORT(rmat);
  REPORT(M);
  REPORT(s_y);
  REPORT(sigma);
  REPORT(Linf);
  REPORT(K);
  REPORT(CV);
  REPORT(pen);
  
  ADREPORT(s_y);
  ADREPORT(sigma);
  ADREPORT(Linf);
  ADREPORT(K);
  ADREPORT(CV);
  
  return nll;
}
