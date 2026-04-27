#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
List sequential_kf_fast(NumericVector m_hat, 
                        NumericMatrix s_hat, 
                        NumericVector y, 
                        NumericMatrix Hm, 
                        NumericVector R_diag, 
                        IntegerVector NA_ind) {
  
  int m = m_hat.size();
  int n_stations = y.size();
  double log_lik = 0.0;
  const double log2pi = std::log(2.0 * M_PI);
  
  std::vector<bool> is_na(n_stations, false);
  for(int k = 0; k < NA_ind.size(); ++k) {
    if(NA_ind[k] > 0 && NA_ind[k] <= n_stations) is_na[NA_ind[k] - 1] = true;
  }

  for(int j = 0; j < n_stations; ++j) {
    if(is_na[j]) continue;

    double dy_j = y[j];
    for(int k = 0; k < m; ++k) {
      dy_j -= Hm(j, k) * m_hat[k];
    }

    NumericVector a(m);
    for(int r = 0; r < m; ++r) {
      double sum = 0;
      for(int c = 0; c < m; ++c) {
        sum += s_hat(r, c) * Hm(j, c);
      }
      a[r] = sum;
    }

    double sj = R_diag[j];
    for(int k = 0; k < m; ++k) {
      sj += Hm(j, k) * a[k];
    }

    log_lik -= 0.5 * (log2pi + std::log(sj) + (dy_j * dy_j) / sj);

    for(int k = 0; k < m; ++k) {
      double kj = a[k] / sj;
      m_hat[k] += kj * dy_j;
      
      for(int l = 0; l < m; ++l) {
        s_hat(k, l) -= (a[k] * a[l]) / sj;
      }
    }
  }

  return List::create(
    _["m_hat"]   = m_hat,
    _["s_hat"]   = s_hat,
    _["log_lik"] = log_lik
  );
}