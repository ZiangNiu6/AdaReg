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

