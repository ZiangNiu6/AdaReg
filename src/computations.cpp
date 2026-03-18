#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double bdm_stat_cpp(NumericVector Y, IntegerVector A,
                    IntegerVector batch_id) {
  // Identify unique batches
  IntegerVector batches = sort_unique(batch_id);
  int B = batches.size();
  int n = Y.size();
  double bdm_sum = 0.0;
  double sqrtB = std::sqrt((double)B);

  for (int bi = 0; bi < B; bi++) {
    int b = batches[bi];
    // Compute batch-level means, counts, and variances in one pass
    double sum1 = 0.0, sum2 = 0.0;
    int N1 = 0, N2 = 0;
    for (int i = 0; i < n; i++) {
      if (batch_id[i] == b) {
        if (A[i] == 1) {
          sum1 += Y[i];
          N1++;
        } else {
          sum2 += Y[i];
          N2++;
        }
      }
    }
    if (N1 < 2 || N2 < 2) continue;

    double mean1 = sum1 / N1;
    double mean2 = sum2 / N2;
    double S_b = mean1 - mean2;

    // Compute variances in second pass
    double ss1 = 0.0, ss2 = 0.0;
    for (int i = 0; i < n; i++) {
      if (batch_id[i] == b) {
        if (A[i] == 1) {
          double d = Y[i] - mean1;
          ss1 += d * d;
        } else {
          double d = Y[i] - mean2;
          ss2 += d * d;
        }
      }
    }
    double var1 = ss1 / N1;
    double var2 = ss2 / N2;
    double se_b = std::sqrt(var1 / N1 + var2 / N2);
    if (se_b <= 0.0) continue;

    bdm_sum += S_b / (sqrtB * se_b);
  }

  return bdm_sum;
}

// [[Rcpp::export]]
List ipw_cpp(NumericVector Y, IntegerVector A,
             NumericVector sampling_prob, IntegerVector batch_id,
             bool one_batch = false, int batch_output = 1) {
  int n = Y.size();
  IntegerVector arm_list = sort_unique(A);
  IntegerVector batch_list = sort_unique(batch_id);
  int n_arms = arm_list.size();
  int n_batches = batch_list.size();

  // Build a map from (batch, arm) -> sampling probability
  // Use first occurrence
  std::map<std::pair<int,int>, double> prob_map;
  for (int i = 0; i < n; i++) {
    auto key = std::make_pair((int)batch_id[i], (int)A[i]);
    if (prob_map.find(key) == prob_map.end()) {
      prob_map[key] = sampling_prob[i];
    }
  }

  // Compute IPW per batch per arm
  // IPW_list[arm_idx][batch_idx]
  NumericMatrix IPW_list(n_arms, n_batches);
  for (int bi = 0; bi < n_batches; bi++) {
    int b = batch_list[bi];
    for (int ai = 0; ai < n_arms; ai++) {
      int a = arm_list[ai];
      double prob = prob_map[std::make_pair(b, a)];
      double sum_y = 0.0;
      for (int i = 0; i < n; i++) {
        if (A[i] == a && batch_id[i] == b) {
          sum_y += Y[i];
        }
      }
      IPW_list(ai, bi) = sum_y / (n * prob);
    }
  }

  NumericVector point(n_arms);
  NumericVector variance(n_arms);

  if (one_batch) {
    // Count observations in the target batch
    int n_batch = 0;
    for (int i = 0; i < n; i++) {
      if (batch_id[i] == batch_output) n_batch++;
    }

    // Compute point estimate using batch size (not total n)
    for (int ai = 0; ai < n_arms; ai++) {
      int a = arm_list[ai];
      int bo_idx = -1;
      for (int bi = 0; bi < n_batches; bi++) {
        if (batch_list[bi] == batch_output) { bo_idx = bi; break; }
      }
      double prob = prob_map[std::make_pair(batch_output, a)];
      double sum_y = 0.0;
      for (int i = 0; i < n; i++) {
        if (A[i] == a && batch_id[i] == batch_output) {
          sum_y += Y[i];
        }
      }
      point[ai] = sum_y / (n_batch * prob);
    }

    // Variance estimate
    for (int ai = 0; ai < n_arms; ai++) {
      int a = arm_list[ai];
      // Construct imputed_Y for this batch
      NumericVector imputed_Y(n_batch, 0.0);
      int j = 0;
      for (int i = 0; i < n; i++) {
        if (batch_id[i] == batch_output) {
          if (A[i] == a) {
            double prob = sampling_prob[i];
            imputed_Y[j] = Y[i] / prob;
          }
          j++;
        }
      }
      double ss = 0.0;
      for (int j2 = 0; j2 < n_batch; j2++) {
        double d = imputed_Y[j2] - point[ai];
        ss += d * d;
      }
      variance[ai] = ss / n_batch;
    }
  } else {
    // Sum across batches
    for (int ai = 0; ai < n_arms; ai++) {
      double s = 0.0;
      for (int bi = 0; bi < n_batches; bi++) {
        s += IPW_list(ai, bi);
      }
      point[ai] = s;
    }

    // Variance
    for (int ai = 0; ai < n_arms; ai++) {
      int a = arm_list[ai];
      // Construct imputed_Y
      NumericVector imputed_Y(n, 0.0);
      for (int i = 0; i < n; i++) {
        if (A[i] == a) {
          imputed_Y[i] = Y[i] / sampling_prob[i];
        }
      }
      double ss = 0.0;
      for (int i = 0; i < n; i++) {
        double d = imputed_Y[i] - point[ai];
        ss += d * d;
      }
      variance[ai] = ss / ((double)n * n);
    }
  }

  // Return as named matrix matching R format
  NumericMatrix output(2, n_arms);
  output(0, _) = point;
  output(1, _) = variance;
  CharacterVector rnames = CharacterVector::create("point", "variance");
  CharacterVector cnames(n_arms);
  for (int i = 0; i < n_arms; i++) cnames[i] = std::to_string(arm_list[i]);
  rownames(output) = rnames;
  colnames(output) = cnames;

  return List::create(Named("output") = output);
}

// [[Rcpp::export]]
NumericMatrix weighted_ipw_cpp(NumericVector Y, IntegerVector A,
                               NumericVector sampling_prob,
                               IntegerVector batch_id,
                               String weighting = "aw") {
  int n = Y.size();
  // Assume 2 arms: 1 and 2
  int n_arms = 2;

  // Build sampling_prob_mat: for each obs, prob of arm 1 and arm 2
  NumericMatrix sp_mat(n, 2);
  for (int i = 0; i < n; i++) {
    if (A[i] == 1) {
      sp_mat(i, 0) = sampling_prob[i];
      sp_mat(i, 1) = 1.0 - sampling_prob[i];
    } else {
      sp_mat(i, 0) = 1.0 - sampling_prob[i];
      sp_mat(i, 1) = sampling_prob[i];
    }
  }

  // Augmented term: Y/prob if assigned, 0 otherwise
  NumericMatrix aug(n, 2);
  for (int i = 0; i < n; i++) {
    if (A[i] == 1) {
      aug(i, 0) = Y[i] / sp_mat(i, 0);
    } else {
      aug(i, 1) = Y[i] / sp_mat(i, 1);
    }
  }

  // Weight constant
  NumericMatrix wc(n, 2);
  bool is_aw = (weighting == "aw");
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < 2; j++) {
      wc(i, j) = is_aw ? std::sqrt(sp_mat(i, j) / n) : sp_mat(i, j) / n;
    }
  }

  // Sum of weights
  NumericVector sum_wc(2, 0.0);
  for (int i = 0; i < n; i++) {
    sum_wc[0] += wc(i, 0);
    sum_wc[1] += wc(i, 1);
  }

  // Point estimate
  NumericVector point(2, 0.0);
  for (int i = 0; i < n; i++) {
    point[0] += aug(i, 0) * wc(i, 0);
    point[1] += aug(i, 1) * wc(i, 1);
  }
  point[0] /= sum_wc[0];
  point[1] /= sum_wc[1];

  // Variance estimate
  NumericVector var_est(2, 0.0);
  for (int i = 0; i < n; i++) {
    double d0 = (aug(i, 0) - point[0]) * wc(i, 0);
    double d1 = (aug(i, 1) - point[1]) * wc(i, 1);
    var_est[0] += d0 * d0;
    var_est[1] += d1 * d1;
  }
  var_est[0] /= (sum_wc[0] * sum_wc[0]);
  var_est[1] /= (sum_wc[1] * sum_wc[1]);

  // Build output matrix
  NumericMatrix output(2, 2);
  output(0, 0) = point[0];
  output(0, 1) = point[1];
  output(1, 0) = var_est[0];
  output(1, 1) = var_est[1];
  CharacterVector rn = CharacterVector::create("point", "variance");
  CharacterVector cn = CharacterVector::create("1", "2");
  rownames(output) = rn;
  colnames(output) = cn;

  return output;
}

// Helper: generate standard normal using Box-Muller
static void rnorm_pair(double &z1, double &z2) {
  double u1 = R::runif(0.0, 1.0);
  double u2 = R::runif(0.0, 1.0);
  double r = std::sqrt(-2.0 * std::log(u1));
  z1 = r * std::cos(2.0 * M_PI * u2);
  z2 = r * std::sin(2.0 * M_PI * u2);
}

// Generate bivariate normal with given correlation
static void rbvnorm(double rho, double &x1, double &x2) {
  double z1, z2;
  rnorm_pair(z1, z2);
  x1 = z1;
  x2 = rho * z1 + std::sqrt(1.0 - rho * rho) * z2;
}

// [[Rcpp::export]]
NumericVector plugin_bootstrap_cpp(NumericVector mean_estimate,
                                   NumericVector variance_estimate,
                                   double eps, double fs_p,
                                   int B = 2000,
                                   String normalization = "normalized",
                                   double h = 0.0,
                                   String weighting = "aw",
                                   String type = "eps_greedy") {
  double EY0 = mean_estimate[0], EY1 = mean_estimate[1];
  double VY0 = variance_estimate[0], VY1 = variance_estimate[1];
  double EY2_0 = EY0 * EY0 + VY0;
  double EY2_1 = EY1 * EY1 + VY1;

  double H1_0 = fs_p, H1_1 = 1.0 - fs_p;
  double V1_0 = EY2_0 - H1_0 * EY0 * EY0;
  double V1_1 = EY2_1 - H1_1 * EY1 * EY1;

  double cov1 = -std::sqrt((H1_0 * H1_1) / (V1_0 * V1_1)) * EY0 * EY1;
  // cov1 is the correlation since the marginals have variance 1

  // Determine weighting exponent m
  double m;
  if (weighting == "aw") m = 0.5;
  else if (weighting == "cw") m = 0.0;
  else m = 1.0; // dm

  bool is_normalized = (normalization == "normalized");

  NumericVector result(B);

  for (int b = 0; b < B; b++) {
    // Generate A_1 ~ bivariate normal with correlation cov1
    double A1_0, A1_1;
    rbvnorm(cov1, A1_0, A1_1);

    // Compute limiting sampling function for H_2
    double val0 = std::sqrt(V1_0 / H1_0) * A1_0;
    double val1 = std::sqrt(V1_1 / H1_1) * A1_1;
    double diff_value = val0 - val1 + h / std::sqrt(2.0);

    double H2_0;
    if (type == "eps_greedy") {
      H2_0 = (diff_value >= 0.0) ? (1.0 - eps / 2.0) : (eps / 2.0);
    } else { // thompson
      H2_0 = R::pnorm(diff_value, 0.0, 1.0, 1, 0);
      H2_0 = std::max(std::min(H2_0, 1.0 - eps), eps);
    }
    double H2_1 = 1.0 - H2_0;

    double V2_0 = EY2_0 - H2_0 * EY0 * EY0;
    double V2_1 = EY2_1 - H2_1 * EY1 * EY1;
    double cov2 = -std::sqrt((H2_0 * H2_1) / (V2_0 * V2_1)) * EY0 * EY1;

    // Generate A_2 ~ bivariate normal with correlation cov2
    double A2_0, A2_1;
    rbvnorm(cov2, A2_0, A2_1);

    // Compute core statistics
    double R1_0 = std::sqrt(H1_0 / V1_0);
    double R1_1 = std::sqrt(H1_1 / V1_1);
    double R2_0 = std::sqrt(H2_0 / V2_0);
    double R2_1 = std::sqrt(H2_1 / V2_1);

    double denom0 = std::pow(H1_0, m) + std::pow(H2_0, m);
    double denom1 = std::pow(H1_1, m) + std::pow(H2_1, m);

    double M1_0 = 0.5 * std::pow(std::pow(H1_0, m) / (0.5 * denom0), 2.0);
    double M1_1 = 0.5 * std::pow(std::pow(H1_1, m) / (0.5 * denom1), 2.0);
    double M2_0 = 0.5 * std::pow(std::pow(H2_0, m) / (0.5 * denom0), 2.0);
    double M2_1 = 0.5 * std::pow(std::pow(H2_1, m) / (0.5 * denom1), 2.0);

    double w1_0, w1_1, w2_0, w2_1;
    if (is_normalized) {
      double denom_n = std::sqrt(M1_0 / (R1_0*R1_0) + M1_1 / (R1_1*R1_1) +
                                 M2_0 / (R2_0*R2_0) + M2_1 / (R2_1*R2_1));
      w1_0 = std::sqrt(M1_0 / (R1_0*R1_0)) / denom_n;
      w1_1 = std::sqrt(M1_1 / (R1_1*R1_1)) / denom_n;
      w2_0 = std::sqrt(M2_0 / (R2_0*R2_0)) / denom_n;
      w2_1 = std::sqrt(M2_1 / (R2_1*R2_1)) / denom_n;
    } else {
      w1_0 = std::sqrt(M1_0 / (R1_0*R1_0));
      w1_1 = std::sqrt(M1_1 / (R1_1*R1_1));
      w2_0 = std::sqrt(M2_0 / (R2_0*R2_0));
      w2_1 = std::sqrt(M2_1 / (R2_1*R2_1));
    }

    // diff_operator = c(1, -1)
    result[b] = (A1_0 * w1_0 - A1_1 * w1_1) + (A2_0 * w2_0 - A2_1 * w2_1);
  }

  return result;
}
