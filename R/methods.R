# This is a Rscript including all the tests considered

#' Test with sample splitting and computed on first batch
#'
#' @param data Data to be analyzed
#'
#' @return left, right and both sided p-values
#' @export

first_batch_only <- function(data){

  # compute the test statistic
  data_to_analyze <- data$data_to_analyze
  IPW_pt <- IPW(data_to_analyze, one_batch = TRUE, batch_output = 1)
  point_estimate <- IPW_pt["point", 1] - IPW_pt["point", 2]
  sd_estimate <- sqrt(IPW_pt["variance", 1] + IPW_pt["variance", 2] + 2 * prod(IPW_pt["point", ]))

  # compute p-value
  n <- length(data_to_analyze$arm)
  p_value_left <- pnorm(sqrt(n / 2) * point_estimate / sd_estimate)
  p_value_right <- pnorm(sqrt(n / 2) * point_estimate / sd_estimate, lower.tail = FALSE)
  p_value_both <- 2 * min(p_value_left, p_value_right)

  # output
  return(list(
    left = p_value_left,
    right = p_value_right,
    both = p_value_both
  ))
}


#' Test with sample splitting
#'
#' @param data Data to be analyzed
#'
#' @return left, right and both sided p-values
#' @export

second_batch_only <- function(data){

  # compute the test statistic
  data_to_analyze <- data$data_to_analyze
  IPW_pt <- IPW(data_to_analyze, one_batch = TRUE, batch_output = 2)
  point_estimate <- IPW_pt["point", 1] - IPW_pt["point", 2]
  sd_estimate <- sqrt(IPW_pt["variance", 1] + IPW_pt["variance", 2] + 2 * prod(IPW_pt["point", ]))

  # compute p-value
  n <- length(data_to_analyze$arm)
  p_value_left <- pnorm(sqrt(n / 2) * point_estimate / sd_estimate)
  p_value_right <- pnorm(sqrt(n / 2) * point_estimate / sd_estimate, lower.tail = FALSE)
  p_value_both <- 2 * min(p_value_left, p_value_right)

  # output
  return(list(
    left = p_value_left,
    right = p_value_right,
    both = p_value_both
  ))
}

#' IPW test with normal quantile
#'
#' @param data Data to be analyzed
#'
#' @return left, right and both sided tests
#' @export

normal_test <- function(data){

  # compute the test statistic
  data_to_analyze <- data$data_to_analyze
  IPW_pt <- IPW(data_to_analyze)
  n <- length(data_to_analyze$arm)
  point_estimate <- IPW_pt["point", 1] - IPW_pt["point", 2]
  sd_estimate <- sqrt(IPW_pt["variance", 1] + IPW_pt["variance", 2] + 2 * prod(IPW_pt["point",]) / n)

  # compute p-value
  p_value_left <- pnorm(point_estimate / sd_estimate)
  p_value_right <- pnorm(point_estimate / sd_estimate, lower.tail = FALSE)
  p_value_both <- 2 * min(p_value_left, p_value_right)

  # output
  return(list(
    left = p_value_left,
    right = p_value_right,
    both = p_value_both
  ))
}

#' Unnormalized adaptive weighting test
#'
#' @param data Data to be analyzed
#'
#' @return left, right and both sided tests
#' @export

UA_test <- function(data){

  # compute the weighted test statistic
  data_to_analyze <- data$data_to_analyze
  IPW_pt <- weighted_IPW(data_to_analyze)
  point_estimate <- IPW_pt["point", 1] - IPW_pt["point", 2]

  # compute the variance and second moments of potential outcomes
  data_to_analyze$reward <- data_to_analyze$reward^2
  IPW_sq <- weighted_IPW(data_to_analyze)
  variance_Y <- c(IPW_sq["point", 1] - IPW_pt["point", 1]^2,
                  IPW_sq["point", 2] - IPW_pt["point", 2]^2)
  
  # impute the negative variance by an upper bound
  idx_negative <- which(variance_Y <= 0)
  variance_Y[idx_negative] <- IPW_sq["point", idx_negative]

  # if the estimated variance is negative, output the NA p-values
  if(any(variance_Y <= 0)){
    return(list(
      left = NA,
      right = NA,
      both = NA
    ))
  }

  # compute the bootstrap sample
  eps <- data$eps
  initial_prob <- data$initial_prob
  candidate_sample <- plugin_bootstrap(mean_estimate = IPW_pt["point", ],
                                       variance_estimate = variance_Y,
                                       eps = eps, fs_p = initial_prob[1], normalization = "unnormalized",
                                       weighting = "aw", type = data$type, B = 5000)

  # compute p-value
  n <- length(data_to_analyze$arm)
  p_value_left <- (length(which(candidate_sample <= (sqrt(n) * point_estimate))) + 1) / (length(candidate_sample) + 1)
  p_value_right <- (length(which(candidate_sample >= (sqrt(n) * point_estimate))) + 1) / (length(candidate_sample) + 1)
  p_value_both <- 2 * min(p_value_left, p_value_right)

  # output
  return(list(
    left = p_value_left,
    right = p_value_right,
    both = p_value_both
  ))
}

#' Normalized adaptive weighting test
#'
#' @param data Data to be analyzed
#'
#' @return left, right and both sided tests
#' @export

NA_test <- function(data){

  # compute the weighted test statistic
  data_to_analyze <- data$data_to_analyze
  IPW_pt <- weighted_IPW(data_to_analyze)
  point_estimate <- IPW_pt["point", 1] - IPW_pt["point", 2]
  sd_estimate <- sqrt(IPW_pt["variance", 1] + IPW_pt["variance", 2])

  # compute the variance and second moments of potential outcomes
  data_to_analyze$reward <- data_to_analyze$reward^2
  IPW_sq <- weighted_IPW(data_to_analyze)
  variance_Y <- c(IPW_sq["point", 1] - IPW_pt["point", 1]^2,
                  IPW_sq["point", 2] - IPW_pt["point", 2]^2)
  
  # impute the negative variance by an upper bound
  idx_negative <- which(variance_Y <= 0)
  variance_Y[idx_negative] <- IPW_sq["point", idx_negative]

  # if the estimated variance is negative, output the NA p-values
  if(any(variance_Y <= 0)){
    return(list(
      left = NA,
      right = NA,
      both = NA
    ))
  }

  # compute the bootstrap sample
  eps <- data$eps
  initial_prob <- data$initial_prob
  candidate_sample <- plugin_bootstrap(mean_estimate = IPW_pt["point", ],
                                       variance_estimate = variance_Y,
                                       eps = eps, fs_p = initial_prob[1], normalization = "normalized",
                                       weighting = "aw", type = data$type, B = 5000)

  # compute p-value
  p_value_left <- (length(which(candidate_sample <= (point_estimate / sd_estimate))) + 1) / (length(candidate_sample) + 1)
  p_value_right <- (length(which(candidate_sample >= (point_estimate / sd_estimate))) + 1) / (length(candidate_sample) + 1)
  p_value_both <- 2 * min(p_value_left, p_value_right)

  # output
  return(list(
    left = p_value_left,
    right = p_value_right,
    both = p_value_both
  ))
}


#' Unnormalized constant weighting test
#'
#' @param data Data to be analyed
#'
#' @return left, right and both sided tests
#' @export

UC_test <- function(data){

  # compute the test statistic
  data_to_analyze <- data$data_to_analyze
  IPW_pt <- IPW(data_to_analyze)
  point_estimate <- IPW_pt["point", 1] - IPW_pt["point", 2]

  # compute the variance and second moments of potential outcomes
  data_to_analyze$reward <- data_to_analyze$reward^2
  IPW_sq <- weighted_IPW(data_to_analyze)
  variance_Y <- c(IPW_sq["point", 1] - IPW_pt["point", 1]^2,
                  IPW_sq["point", 2] - IPW_pt["point", 2]^2)
  
  # impute the negative variance by an upper bound
  idx_negative <- which(variance_Y <= 0)
  variance_Y[idx_negative] <- IPW_sq["point", idx_negative]

  # compute the bootstrap sample for unnormalized test statistic
  eps <- data$eps
  initial_prob <- data$initial_prob
  candidate_sample <- plugin_bootstrap(mean_estimate = IPW_pt["point", ],
                                       variance_estimate = variance_Y,
                                       eps = eps, fs_p = initial_prob[1], normalization = "unnormalized",
                                       weighting = "cw", type = data$type, B = 5000)

  # compute unnormalized ci interval
  n <- length(data_to_analyze$arm)
  p_value_left <- (length(which(candidate_sample <= (sqrt(n) * point_estimate))) + 1) / (length(candidate_sample) + 1)
  p_value_right <- (length(which(candidate_sample >= (sqrt(n) * point_estimate))) + 1) / (length(candidate_sample) + 1)
  p_value_both <- 2 * min(p_value_left, p_value_right)

  # output
  return(list(
    left = p_value_left,
    right = p_value_right,
    both = p_value_both
  ))
}

#' Normalized constant weighting test
#'
#' @param data Data to be analyzed
#'
#' @return left, right and both sided tests
#' @export

NC_test <- function(data){

  # compute the test statistic
  data_to_analyze <- data$data_to_analyze
  IPW_pt <- IPW(data_to_analyze)
  point_estimate <- IPW_pt["point", 1] - IPW_pt["point", 2]
  sd_estimate <- sqrt(IPW_pt["variance", 1] + IPW_pt["variance", 2])

  # compute the variance and second moments of potential outcomes
  data_to_analyze$reward <- data_to_analyze$reward^2
  IPW_sq <- weighted_IPW(data_to_analyze)
  variance_Y <- c(IPW_sq["point", 1] - IPW_pt["point", 1]^2,
                  IPW_sq["point", 2] - IPW_pt["point", 2]^2)
  
  # impute the negative variance by an upper bound
  idx_negative <- which(variance_Y <= 0)
  variance_Y[idx_negative] <- IPW_sq["point", idx_negative]

  # compute the bootstrap sample
  eps <- data$eps
  initial_prob <- data$initial_prob
  candidate_sample <- plugin_bootstrap(mean_estimate = IPW_pt["point", ],
                                       variance_estimate = variance_Y,
                                       eps = eps, fs_p = initial_prob[1], normalization = "normalized",
                                       weighting = "cw", type = data$type, B = 5000)

  # compute two types of confidence interval
  p_value_left <- (length(which(candidate_sample <= (point_estimate / sd_estimate))) + 1) / (length(candidate_sample) + 1)
  p_value_right <- (length(which(candidate_sample >= (point_estimate / sd_estimate))) + 1) / (length(candidate_sample) + 1)
  p_value_both <- 2 * min(p_value_left, p_value_right)

  # output
  return(list(
    left = p_value_left,
    right = p_value_right,
    both = p_value_both
  ))
}


#' Test with unnormalized difference-in-means test statistic
#'
#' @param data Data to be analyzed
#'
#' @return left, right and both sided tests
#' @export

UDM_test <- function(data){
  
  # compute the weighted test statistic
  data_to_analyze <- data$data_to_analyze
  IPW_pt <- weighted_IPW(data_to_analyze, weighting = "dm")
  point_estimate <- IPW_pt["point", 1] - IPW_pt["point", 2]
  
  # compute the variance and second moments of potential outcomes
  data_to_analyze$reward <- data_to_analyze$reward^2
  IPW_sq <- weighted_IPW(data_to_analyze, weighting = "dm")
  variance_Y <- c(IPW_sq["point", 1] - IPW_pt["point", 1]^2,
                  IPW_sq["point", 2] - IPW_pt["point", 2]^2)
  
  # impute the negative variance by an upper bound
  idx_negative <- which(variance_Y <= 0)
  variance_Y[idx_negative] <- IPW_sq["point", idx_negative]
  
  # if the estimated variance is negative, output the NA p-values
  if(any(variance_Y <= 0)){
    return(list(
      left = NA,
      right = NA,
      both = NA
    ))
  }
  
  # compute the bootstrap sample
  eps <- data$eps
  initial_prob <- data$initial_prob
  candidate_sample <- plugin_bootstrap(mean_estimate = IPW_pt["point", ],
                                       variance_estimate = variance_Y,
                                       eps = eps, fs_p = initial_prob[1], normalization = "unnormalized",
                                       weighting = "dm", type = data$type, B = 5000)
  
  # compute p-value
  n <- length(data_to_analyze$arm)
  p_value_left <- (length(which(candidate_sample <= (sqrt(n) * point_estimate))) + 1) / (length(candidate_sample) + 1)
  p_value_right <- (length(which(candidate_sample >= (sqrt(n) * point_estimate))) + 1) / (length(candidate_sample) + 1)
  p_value_both <- 2 * min(p_value_left, p_value_right)
  
  # output
  return(list(
    left = p_value_left,
    right = p_value_right,
    both = p_value_both
  ))
}

#' Test with normalized difference-in-means test statistic
#'
#' @param data Data to be analyzed
#'
#' @return left, right and both sided tests
#' @export

NDM_test <- function(data){
  
  # compute the weighted test statistic
  data_to_analyze <- data$data_to_analyze
  IPW_pt <- weighted_IPW(data_to_analyze, weighting = "dm")
  point_estimate <- IPW_pt["point", 1] - IPW_pt["point", 2]
  sd_estimate <- sqrt(IPW_pt["variance", 1] + IPW_pt["variance", 2])
  
  # compute the variance and second moments of potential outcomes
  data_to_analyze$reward <- data_to_analyze$reward^2
  IPW_sq <- weighted_IPW(data_to_analyze, weighting = "dm")
  variance_Y <- c(IPW_sq["point", 1] - IPW_pt["point", 1]^2,
                  IPW_sq["point", 2] - IPW_pt["point", 2]^2)
  
  # impute the negative variance by an upper bound
  idx_negative <- which(variance_Y <= 0)
  variance_Y[idx_negative] <- IPW_sq["point", idx_negative]
  
  # if the estimated variance is negative, output the NA p-values
  if(any(variance_Y <= 0)){
    return(list(
      left = NA,
      right = NA,
      both = NA
    ))
  }
  
  # compute the bootstrap sample
  eps <- data$eps
  initial_prob <- data$initial_prob
  candidate_sample <- plugin_bootstrap(mean_estimate = IPW_pt["point", ],
                                       variance_estimate = variance_Y,
                                       eps = eps, fs_p = initial_prob[1], normalization = "normalized",
                                       weighting = "dm", type = data$type, B = 5000)
  
  # compute two types of confidence interval
  p_value_left <- (length(which(candidate_sample <= (point_estimate / sd_estimate))) + 1) / (length(candidate_sample) + 1)
  p_value_right <- (length(which(candidate_sample >= (point_estimate / sd_estimate))) + 1) / (length(candidate_sample) + 1)
  p_value_both <- 2 * min(p_value_left, p_value_right)

  # output
  return(list(
    left = p_value_left,
    right = p_value_right,
    both = p_value_both
  ))
}

#' Batchwise Difference-in-Means (BDM) test
#'
#' Normalizes the difference-in-means per batch and aggregates across batches.
#' Under the null, the test statistic is asymptotically standard normal.
#' See Hirano and Porter (2025), Appendix D.
#'
#' @param data Data to be analyzed
#'
#' @return left, right and both sided p-values
#' @export

BDM_test <- function(data){

  # extract data and compute test statistic via Rcpp
  data_to_analyze <- data$data_to_analyze
  bdm_sum <- bdm_stat_cpp(
    as.numeric(data_to_analyze$reward),
    as.integer(data_to_analyze$arm),
    as.integer(data_to_analyze$batch_id)
  )

  # compute p-values using standard normal
  p_value_left <- pnorm(bdm_sum)
  p_value_right <- pnorm(bdm_sum, lower.tail = FALSE)
  p_value_both <- 2 * min(p_value_left, p_value_right)

  # output
  return(list(
    left = p_value_left,
    right = p_value_right,
    both = p_value_both
  ))
}

#' Concentration-based inference test
#'
#' Non-asymptotic test using self-normalized martingale concentration bounds
#' from Abbasi-Yadkori, Pal, and Szepesvari (2011). P-values are obtained by
#' inverting the simultaneous armwise confidence intervals (paper Eq. 14).
#' Valid under adaptive sampling without requiring bootstrap or variance estimation.
#'
#' @param data Data to be analyzed
#'
#' @return left, right and both sided p-values
#' @export

concentration_test <- function(data){

  dat <- data$data_to_analyze
  A <- dat$arm
  Y <- dat$reward
  K <- 2

  N1 <- sum(A == 1)
  N2 <- sum(A == 2)
  Xbar1 <- mean(Y[A == 1])
  Xbar2 <- mean(Y[A == 2])
  Dhat <- Xbar1 - Xbar2

  # armwise confidence radius from Eq. (14) with K arms
  radius <- function(N, alpha) {
    sqrt(((1 + N) / N^2) * (1 + 2 * log(K * sqrt(1 + N) / alpha)))
  }

  # total radius for the gap
  R <- function(alpha) {
    radius(N1, alpha) + radius(N2, alpha)
  }

  # invert R(alpha) to find p-value via bisection
  invert_alpha <- function(target) {
    if (target <= 0) return(1)
    if (target <= R(1 - 1e-12)) return(1)
    lo <- 1e-12
    hi <- 1 - 1e-12
    for (i in 1:60) {
      mid <- (lo + hi) / 2
      if (target > R(mid)) {
        hi <- mid
      } else {
        lo <- mid
      }
    }
    return(hi)
  }

  p_value_right <- invert_alpha(Dhat)
  p_value_left <- invert_alpha(-Dhat)
  p_value_both <- invert_alpha(abs(Dhat))

  return(list(
    left = p_value_left,
    right = p_value_right,
    both = p_value_both
  ))
}

