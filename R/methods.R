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

  # extract data
  data_to_analyze <- data$data_to_analyze
  Y <- data_to_analyze$reward
  A <- data_to_analyze$arm
  batch_id <- data_to_analyze$batch_id
  batch_list <- sort(unique(batch_id))
  B <- length(batch_list)

  # compute per-batch statistics and aggregate
  bdm_sum <- 0
  for(b in batch_list){
    idx <- which(batch_id == b)
    idx_1 <- idx[A[idx] == 1]
    idx_2 <- idx[A[idx] == 2]
    N_b1 <- length(idx_1)
    N_b2 <- length(idx_2)

    # skip batch if either arm has fewer than 2 observations
    if(N_b1 < 2 || N_b2 < 2) next

    # batch-level means
    beta_b1 <- mean(Y[idx_1])
    beta_b2 <- mean(Y[idx_2])

    # batch-level difference in means
    S_b <- beta_b1 - beta_b2

    # batch-level variances
    sigma2_b1 <- sum((Y[idx_1] - beta_b1)^2) / N_b1
    sigma2_b2 <- sum((Y[idx_2] - beta_b2)^2) / N_b2

    # batch-level weight
    se_b <- sqrt(sigma2_b1 / N_b1 + sigma2_b2 / N_b2)
    if(se_b <= 0) next

    omega_b <- 1 / (sqrt(B) * se_b)

    bdm_sum <- bdm_sum + omega_b * S_b
  }

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

