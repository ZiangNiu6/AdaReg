#' Sampling algorithm
#'
#' @param eps A small probability that is used to draw random sample from arms
#' @param reward A vector of aggregated rewards to determine which arm to sample
#' @param n Number of data that is needed to sample
#' @param dist_fam Distribution that is used for sampling
#' @param dist_params Parameters for the distribution, e.g., mean or variance
#' @param data_generation A logical value determing if we want to generate the data or not
#' @param sample_size sample size used to calculate expected reward estimate; can equal n
#' @param type sampling algorithm
#'
#' @return the random sampled and random index with length n
#' @importFrom stats rnorm rbinom pnorm rt rpois
#' @export
sampling_function <- function(eps, reward, n, dist_fam = NULL, dist_params = NULL, data_generation = TRUE,
                              sample_size = NULL, type = "eps_greedy"){

  # specify an empty vector of result
  sample_id <- numeric(n)

  # compute the number of distinct arms
  num_arm <- length(reward)

  # decide the sampling scheme
  switch (type,
    eps_greedy = {
      # hard threshold
      # compute the random probability given the reward
      max_id <- rbinom(1, 1, 0.5) + 1
      prob_estimate <- sapply(1:num_arm,
                              function(x){
                                if(all(reward == max(reward))){
                                  return(dplyr::if_else(x == max_id, 1, 0)*(1-eps) + eps / num_arm)
                                }else{
                                  dplyr::if_else(reward[x] == max(reward), 1, 0)*(1-eps) + eps / num_arm
                                }
                              })

      # decide the index of random sample
      random_coin <- rbinom(n, 1, eps)
      num_random_sample <- sum(random_coin)
      sample_id[which(random_coin == 1)] <- sapply(1:num_random_sample,
                                                   function(x) sample(1:num_arm, 1))
      sample_id[which(random_coin == 0)] <- dplyr::if_else(all(reward == max(reward)), 
                                                           max_id, 
                                                           min(which(reward == max(reward))))
    },
    thompson = {

      # alternative smooth version
      scaled_reward <- sqrt(sample_size) * reward
      # compute the random probability given the reward
      prob_estimate <- sapply(1:num_arm,
                              function(x){
                                dplyr::if_else(scaled_reward[x] == max(scaled_reward),
                                               min(pnorm(max(scaled_reward) - min(scaled_reward)), 1 - eps),
                                               max(pnorm(min(scaled_reward) - max(scaled_reward)), eps))
                              })
      # decide the index of random sample
      sample_id <- rbinom(n, size = 1, prob = prob_estimate[2]) + 1
    }
  )

  # return the result if outcome is not needed
  if(!data_generation){
    return(list(
      sample_id = sample_id,
      prob_estimate = prob_estimate
    ))
  }else{
    if(is.null(dist_fam) | is.null(dist_params)){
      stop("Distribution family and parameters should be specified given the data generation is required!")
    }
    
    # sample according to dist_fam
    switch (dist_fam,
            gaussian = {
              if(ncol(dist_params) != 2){
                stop("Both mean and variance should be provided to sample!")
              }
              mu <- dist_params[, 1]
              sigma <- dist_params[, 2]
              random_sample <- sapply(sample_id,
                                      function(x) rnorm(n = 1,
                                                        mean = mu[x],
                                                        sd = sigma[x]))
            },
            Student = {
              if(ncol(dist_params) != 2){
                stop("Both mean and df should be provided to sample!")
              }
              mu <- dist_params[, 1]
              df <- dist_params[, 2]
              random_sample <- sapply(sample_id,
                                      function(x) {rt(n = 1, df = df[x]) + mu[x]})
            },
            mix_normal = {
              if(length(dist_params) != 3){
                stop("All mean, sd and mixture prob should be provided to sample!")
              }
              # extract distribution information
              mu <- dist_params$mu
              sigma <- dist_params$sd
              weights <- dist_params$prob
              
              # generate second stage outcome
              random_sample <- sapply(sample_id,
                                      function(x) {
                                        rmvnormmix(n = 1, 
                                                   lambda = weights[, x], 
                                                   mu = mu[, x], 
                                                   sigma = sigma[, x])
                                      })
            },
            Bernoulli = {
              if(ncol(dist_params) != 1){
                stop("Mean should be provided to sample!")
              }
              mu <- dist_params[, 1]
              random_sample <- sapply(sample_id,
                                      function(x) rbinom(n = 1, size = 1,
                                                         prob = mu[x]))
            },
            Poisson = {
              if(ncol(dist_params) != 1){
                stop("Lambda should be provided to sample!")
              }
              lambda <- dist_params[, 1]
              random_sample <- sapply(sample_id,
                                      function(x) rpois(n = 1, lambda = lambda[x]))
            }
    )
    return(list(
      random_sample = random_sample,
      random_index = sample_id,
      prob_estimate = prob_estimate
    ))
  }
}



#' Function fitting GLM regression
#'
#' @param data Data containing covariate and reward for the specific arm
#' @param hyperparams Hyperparmeters containing the family
#'
#' @return Same output as in fit_lasso function
#' @importFrom stats glm
#' @export
fit_glm <- function(data, hyperparams){
  # extract the covariate and reward
  X <- data$covariate
  Y <- data$reward

  # extract the family hyperparameter
  family <- hyperparams$family

  # check if there is data distribution there or not
  if(is.null(family)){
    stop("Specification of GLM should be announced with name family!")
  }

  # use the cross-validation for model selection
  if(is.null(X) | is.null(ncol(X))){
    glm_fit <- glm(Y ~ 1, family = family)
  }else{
    glm_fit <- glm(Y ~ X, family = family)
  }

  # obtain the optimal coefficient, linear predictions and conditional mean
  optimal_coef <- unname(coef(glm_fit))
  linear_pred <- unname(predict(glm_fit, newx = X, type = "link"))
  transformed_pred <- unname(predict(glm_fit, newx = X, type = "response"))

  # obtain the empirical l2 error
  insample_err <- sum(linear_pred^2) / length(linear_pred)

  # return the output
  return(list(
    conditional_mean = transformed_pred,
    coefficient = optimal_coef,
    insampled_err = insample_err,
    linear_pred = linear_pred
  ))
}


#' lasso fit
#'
#' @param data Data containing covariate and reward
#' @param hyperparams Hyperparameter for lasso regression
#'
#' @return A list of output containing the estimated coef, conditional mean etc
#' @importFrom stats coef predict
#' @importFrom glmnet cv.glmnet
#' @export
fit_lasso <- function(data, hyperparams){

  # extract the covariate and reward
  X <- data$covariate
  Y <- data$reward

  # extract the family and lambda hyperparameters
  family <- hyperparams$family
  lambda <- hyperparams$lambda
  # check if there is data distribution there or not
  if(is.null(family)){
    stop("Specification of GLM should be announced with name family!")
  }
  # check if there is choice of regularization parameter
  if(is.null(lambda)){
    stop("Regularization type should be announced with name lambda!")
  }

  # use the cross-validation for model selection
  cv_fit <- glmnet::cv.glmnet(X, Y, family = family)

  # obtain the optimal coefficient, linear predictions and conditional mean
  optimal_coef <- unname(coef(cv_fit, s = lambda))
  linear_pred <- unname(predict(cv_fit, newx = X,
                                s = lambda, type = "link"))
  transformed_pred <- unname(predict(cv_fit, newx = X,
                                     s = lambda, type = "response"))

  # obtain the empirical l2 error
  insample_err <- sum(linear_pred^2) / length(linear_pred)

  # return the output
  return(list(
    conditional_mean = transformed_pred,
    coefficient = optimal_coef,
    insampled_err = insample_err,
    linear_pred = linear_pred
  ))
}

#' Variance estimate
#'
#' @param data Data including covariate, arm and reward
#' @param conditional_mean Conditional mean estimate for this specific arm
#' @param mean_estimate The estimate for mean of potential outcome
#' @param arm The arm to be evaluated
#'
#' @return the variance estimate
#' @export
variance_estimate <- function(data,
                              conditional_mean,
                              mean_estimate, arm){
  if(!is.numeric(arm)){
    stop("Variance for which arm should be specified!")
  }

  Y <- data$reward
  A <- data$arm
  X <- data$covariate
  n <- length(A)

  # extract the arm id
  arm_id <- which(A == arm)

  # compute the var_F
  var_F <- mean(((Y - conditional_mean)^2)[arm_id])

  # compute the second part of the variance estimate
  var_B <- mean((conditional_mean - mean_estimate)^2)

  # return the variance estimate
  return(list(var_F = var_F, var_B = var_B))
}

#' Covariance estimate
#'
#' @param conditional_mean Conditional mean vector of length n
#' @param mean_estimate Vector of mean estimate of two arms
#'
#' @return Covariance estimate
#' @export
covariance_estimate <- function(conditional_mean, mean_estimate){
  # Check if two arm id is specified or not
  if(length(mean_estimate) != 2 | ncol(conditional_mean) != 2){
    stop("Only two arms should be used to compute the covariance!")
  }

  n <- nrow(conditional_mean)

  # compute the covariance
  residual_mat <- conditional_mean - matrix(rep(mean_estimate, n),
                                            nrow = n, byrow = TRUE)
  covariance_estimate <- mean(apply(residual_mat, 1, prod))
  return(covariance_estimate)
}



#' Data generation function
#'
#' @param eps The input for sampling algorithm
#' @param dist_fam A vector storing the type of sampling distribution
#' @param n_1 Number of sample in stage 1
#' @param n_2 Number of sample in stage 2
#' @param p Dimension of the covariate
#' @param coef_list List of coefficients
#' @param deciding_scheme IPW or AIPW
#' @param reg_type NULL as default choice
#' @param hyperparams Null as default choice
#' @param dist_params A matrix storing the hyperparameters for sampling
#' @param initial_prob A vector storing the initial sampling probability for each arm
#' @param link When it is linear, linear model is considered; when it is nonlinear, nonlinear model is considered
#' @param type sampling algorithm choice; either epsilon greedy or Thompson sampling
#'
#' @return Data list containing reward, arm, covariate and batch ids
#' @importFrom stats rnorm rmultinom
#' @importFrom mixtools rmvnormmix
#' @import dplyr
#' @export
data_generate_online <- function(eps, dist_fam, n_1, n_2, p, coef_list = NULL,
                                 deciding_scheme = "IPW",
                                 reg_type = NULL, hyperparams = NULL,
                                 dist_params, initial_prob, link = "linear", type = "eps_greedy"){

  if(any(sum(initial_prob) != 1 | initial_prob <= 0)){
    stop("The sum of prob vector is one and prob should be positive!")
  }
  n <- n_1 + n_2
  batch_1 <- 1:n_1
  batch_2 <- setdiff(1:n, batch_1)

  # generating the data
  if(!is.null(coef_list)){
    X <- matrix(rnorm(n*p), nrow = n, ncol = p)
  }
  A <- numeric(n)
  Y <- numeric(n)
  if(is.null(coef_list)){
    linear_pred <- matrix(0, nrow = n, ncol = 2)
  }else{
    if(link == "linear"){
      linear_pred <- unname(sapply(coef_list,
                                   function(b){
                                     X %*% b
                                   }))
    }else{
      linear_pred <- unname(sapply(coef_list,
                                   function(b){
                                     X %*% b
                                   })) + exp(X[, 1])
    }
  }
  
  # assign treatment in the first stage
  A[batch_1] <- apply(rmultinom(n_1, 1, initial_prob), 2, function(x) which(x != 0))
  
  # swith among different distributions
  switch (dist_fam,
    gaussian = {
      
      # extract distribution information
      mu <- dist_params[, 1]
      sigma <- dist_params[, 2]
      
      # first batch
      Y[batch_1] <- sapply(batch_1, function(x){
        rnorm(n = 1,
              mean = mu[A[x]] + linear_pred[x, A[x]],
              sd = sigma[A[x]])
      })
    },
    mix_normal = {
      
      # extract distribution information
      mu <- dist_params$mu
      sigma <- dist_params$sd
      weights <- dist_params$prob
      
      # first batch
      Y[batch_1] <- sapply(batch_1, function(x){
        rmvnormmix(n = 1, 
                   lambda = weights[, A[x]], 
                   mu = mu[, A[x]] + linear_pred[x, A[x]], 
                   sigma = sigma[, A[x]])
      })

    },
    Student = {
      
      # extract distribution information
      mu <- dist_params[, 1]
      df <- dist_params[, 2]
      
      # first batch
      Y[batch_1] <- sapply(batch_1, function(x){
        rt(n = 1, df = df[A[x]]) + mu[A[x]] + linear_pred[x, A[x]]
      })
    },
    Bernoulli = {
      
      # extract distribution information
      mu <- dist_params[, 1]

      # first batch
      Y[batch_1] <- sapply(batch_1, function(x){
        rbinom(n = 1, size = 1, prob = mu[A[x]])
      })
    },
    Poisson = {
      
      # extract distribution information
      lambda <- dist_params[, 1]
      
      # first batch
      Y[batch_1] <- sapply(batch_1, function(x){
        rpois(n = 1, lambda = lambda[A[x]])
      })
    }
  )
  
  # compute the selection rule using data from first stage
  batch_1_estimate <- switch (deciding_scheme,
                              IPW = {
                                # compute the IPW estimate
                                IPW(
                                  data = list(
                                    reward = Y[batch_1],
                                    arm = A[batch_1],
                                    sampling_prob = sapply(batch_1, function(x) initial_prob[A[x]]),
                                    batch_id = rep(1, length(batch_1))
                                  ))["point", ]
                              },
                              AIPW = {
                                if(all(c(is.null(reg_type), is.null(hyperparams)))){
                                  stop("Both reg_type and hyperparams should be specified when AIPW is considered!")
                                }
                                # compute the IPW estimate
                                AIPW(
                                  data = list(
                                    reward = Y[batch_1],
                                    arm = A[batch_1],
                                    covariate = X[batch_1, ],
                                    sampling_prob = sapply(batch_1, function(x) initial_prob[A[x]]),
                                    batch_id = rep(1, length(batch_1))
                                  ), hyperparams = hyperparams, reg_type = reg_type)
                              }
  )
  
  # pass the result to sampling algorithm and get an estimate for the sampling probability
  second_stage_data <- sampling_function(eps = eps, reward = batch_1_estimate,
                                         n = n_2, dist_fam = dist_fam,
                                         sample_size = n_1,
                                         dist_params = dist_params, type = type)
  prob_estimate <- second_stage_data$prob_estimate
  
  # second batch
  A[batch_2] <- second_stage_data$random_index
  Y[batch_2] <- second_stage_data$random_sample + sapply(batch_2,
                                                         function(x){
                                                           linear_pred[x, A[x]]
                                                         })

  # store the probability vector
  sampling_prob <- sapply(1:n, function(x){
    dplyr::if_else(x <= n_1, initial_prob[A[x]], prob_estimate[A[x]])
  })

  # store the data
  if(p == 0){
    # store the data
    data <- list(
      arm = A,
      reward = Y,
      batch_id = c(rep(1, n_1), rep(2, n_2)),
      sampling_prob = sampling_prob
    )
  }else{
    # store the data
    data <- list(
      covariate = X,
      arm = A,
      reward = Y,
      batch_id = c(rep(1, n_1), rep(2, n_2)),
      sampling_prob = sampling_prob
    )
  }

  # return the result
  return(data)
}

#' Generating random variable in the limiting distribution
#'
#' @param A_1 First stage bivariate limiting normal distribution
#' @param V_1 Variance estimate in the first stage
#' @param h Signal (c)
#' @param H_1 Limiting weight vector in the first stage
#' @param eps Epsilon parameter in epsilon-greedy algorithm
#' @param type either epsilon greedy algorithm or Thompson sampling
#' @param weight_vec A weight vector depending how ATE is computed in sampling function (encourage or discourage)
#'
#' @return Probability of choosing specific arm
#' @export
#' @importFrom stats pnorm

limiting_sampling_function <- function(A_1, V_1, h, H_1, eps, type = "eps_greedy",
                                       weight_vec = c(1, -1)){

  # compute the core quantity
  value_vec <- c(sqrt(V_1[1] / H_1[1])*A_1[1], sqrt(V_1[2] / H_1[2])*A_1[2])
  diff_value <-  sum(value_vec * weight_vec)  + h / sqrt(2)

  # denpending on which sampling algorithm is used
  prob_first_arm <-   switch (type,
                              eps_greedy = {
                                if_else(diff_value >= 0, 1 - eps/2, eps / 2)
                              },
                              thompson = {
                                max(min(pnorm(diff_value), 1 - eps), eps)
                              }
  )

  # return the output
  return(c(prob_first_arm, 1 - prob_first_arm))
}


#' Bootstrap function for constant weighting unnormalized test statistic
#'
#' @param mean_estimate Estimate for the expected outcome (bivariate vector)
#' @param variance_estimate Variance estimate (bivariate vector)
#' @param fs_p Treatment assignment probability in the first stage for the first arm
#' @param B Number of bootstrap sample required
#' @param normalization A hyperparameter to use normalization or not in the test statistic
#' @param weighting A hyperparameter to use adaptive or constant weighting
#' @param eps Epsilon tolerance in epsilon-greedy algorithm
#' @param h Signal strength present in the sampling function
#' @param type Sampling algorithm
#' @param weight_vec A weight vector depending how ATE is computed in sampling function (encourage or discourage)
#'
#' @return A vector of bootstrap sample of length B
#' @export
#' @importFrom katlabutils fast_generate_mvn
#' @importFrom dplyr case_when

plugin_bootstrap <- function(mean_estimate,
                             variance_estimate, eps,
                             fs_p, B = 2000, normalization = "normalized",
                             h = 0,
                             weighting = "aw", type = "eps_greedy",
                             weight_vec = c(1, -1)){

  # specify the expectation and variance
  EY_vec <- mean_estimate
  VY_vec <- variance_estimate
  EY2_vec <- EY_vec^2 + VY_vec
  num_arm <- 2

  # append two hyperparameter together
  testing_method <- paste(normalization, weighting, sep = "_")

  # specify the first stage sampling probability
  fs_prob <- c(fs_p, 1 - fs_p)

  # specify the parameter in the limiting distribution from the first stage
  H_1 <- fs_prob
  V_1 <- EY2_vec - H_1*EY_vec^2
  cov_1 <- -sqrt(prod(H_1) / prod(V_1)) * prod(EY_vec)
  covariance_1 <- matrix(c(1, cov_1, cov_1, 1), num_arm, num_arm)

  # generate B sample from the limiting distribution for each h

  ## generate sample from first batch
  A_1 <- fast_generate_mvn(mean = rep(0, num_arm),
                           covariance = covariance_1,
                           num_samples = B)

  # generate sample from the limiting distribution
  final_sample <- sapply(1:B, function(y){
    H_2 <- limiting_sampling_function(A_1[y, ], V_1, h = h, H_1, eps = eps, type = type, weight_vec = weight_vec)
    V_2 <- EY2_vec - H_2*EY_vec^2
    cov_2 <- -sqrt(prod(H_2) / prod(V_2)) * prod(EY_vec)
    covariance_2 <- matrix(c(1, cov_2, cov_2, 1), num_arm, num_arm)

    # generate second stage sample
    A_2 <- fast_generate_mvn(mean = rep(0, num_arm),
                             covariance = covariance_2,
                             num_samples = 1)

    # extract which weighing method is used
    m <- case_when(
      weighting == "aw" ~ 0.5,
      weighting == "cw" ~ 0,
      TRUE ~ 1
    )

    # compute core statistics
    R_1 <- sqrt(H_1 / V_1)
    R_2 <- sqrt(H_2 / V_2)
    M_1 <- 0.5 * (H_1^{m} / (0.5 * sum(H_1^{m})))^2
    M_2 <- 0.5 * (H_2^{m} / (0.5 * sum(H_2^{m})))^2

    # compute the weighting depending on the normalization or not
    switch (normalization,
      normalized = {
        denominator <- sqrt(sum(M_1 / R_1^2 + M_2 / R_2^2))
        w_1 <- sqrt(M_1 / R_1^2) / denominator
        w_2 <- sqrt(M_2 / R_2^2) / denominator
      },
      unnormalized = {
        w_1 <- sqrt(M_1 / R_1^2)
        w_2 <- sqrt(M_2 / R_2^2)
      }
    )


    # # compute the denominator in the weighting vector
    # switch (testing_method,
    #   normalized_aw = {
    #     R_1 <- 1 / (sqrt(H_1[1]) + sqrt(H_2[1]))
    #     R_2 <- 1 / (sqrt(H_1[2]) + sqrt(H_2[2]))
    #     denominator_total <-  sqrt(R_1^2 * (V_1[1] + V_2[1]) + R_2^2 * (V_2[2] + V_1[2]))
    #     denominator_1 <- denominator_2 <- c(denominator_total / R_1, denominator_total / R_2)
    #   },
    #   normalized_cw = {
    #     denominator_1 <- sqrt(sum(V_1 / H_1) + sum(V_2 / H_2)) * sqrt(H_1)
    #     denominator_2 <- sqrt(sum(V_1 / H_1) + sum(V_2 / H_2)) * sqrt(H_2)
    #   },
    #   unnormalized_aw = {
    #     R_1 <- (sqrt(H_1[1]) + sqrt(H_2[1]))
    #     R_2 <- (sqrt(H_1[2]) + sqrt(H_2[2]))
    #     denominator_1 <- denominator_2 <- c(R_1, R_2)
    #   },
    #   unnormalized_cw = {
    #     denominator_1 <- sqrt(H_1)
    #     denominator_2 <- sqrt(H_2)
    #   }
    # )
    #
    # # generate the weight w_1 and w_2
    # w_1 <- sqrt(V_1) / denominator_1
    # w_2 <- sqrt(V_2) / denominator_2

    # define the diff_operator
    diff_operator <- c(1, -1)

    # generate final sample
    return(sum(A_1[y, ]*w_1*diff_operator) + sum(A_2*w_2*diff_operator))
  })

  # return final sample
  return(final_sample)
}


#' AIPW estimator
#'
#' @param data Data containing the reward, arm and covariate
#' @param hyperparams The hyperparameters associated to the regression, an input to fit_glm and fit_lasso
#' @param reg_type Regression methods, e.g., lasso
#' @param return_reg Return the regression coefficient or not
#' @param first_batch TRUE if only data in first batch is sued; FALSE if full data is sued
#'
#' @return A scalar used for downstream inference
#' @export
AIPW <- function(data, hyperparams, reg_type = "lasso", return_reg = TRUE, first_batch = FALSE){

  # extract the arm list
  arm_list <- sort(data$arm |> unique())

  # construct an empty matrix
  output <- matrix(data = NA,
                   nrow = length(arm_list),
                   ncol = 2,
                   dimnames = list(estimate = c("point", "variance"),
                                   arm = arm_list))

  # construct the AIPW estimator
  ## switch depending on the logic value first_batch
  if(first_batch){

    # only consider the first_batch data
    batch_id <- data$batch_id
    data <- lapply(data, function(x) as.matrix(x)[which(batch_id == 1), ])

    # extract the covariate, reward, arm and sampling probability
    X <- data$covariate
    Y <- data$reward
    A <- data$arm
    n <- length(A)

    # construct regression estimator
    reg_result <- sapply(arm_list, function(x){
      arm_id <- which(A == x)
      data_to_reg <- lapply(data, function(x) as.matrix(x)[arm_id, ])
      switch (reg_type,
              lasso = {
                fit_lasso(data = data_to_reg, hyperparams = hyperparams)
              },
              MLE = {
                fit_glm(data = data_to_reg, hyperparams = hyperparams)
              }
      )
    })

    # impute conditional mean
    conditional_mean_imputed <-   switch (hyperparams$family,
                                          gaussian = {
                                            sapply(arm_list,
                                                   function(x){
                                                     coef_fit <- reg_result["coefficient", x]$coefficient
                                                     if(is.null(X)){
                                                       return(rep(as.vector(coef_fit), n))
                                                     }else{
                                                       return(as.vector(cbind(1, X) %*% coef_fit))
                                                     }
                                                   })
                                          }
    )

    # G-estimation term
    g_estimate <- apply(conditional_mean_imputed, 2, mean)

    # augmented term
    residual_vec <- sapply(1:length(A), function(x){
      Y[x] - conditional_mean_imputed[x, A[x]]
    })
    data_aug <- data
    data_aug$reward <- residual_vec
    augmented_term <- IPW(data = data_aug)["point", ]

    # store the mean estimate
    output["point", ] <- augmented_term + g_estimate

    # compute the variance term
    output["variance", ] <- sapply(arm_list,
                                   function(y){
                                     prob_to_use <- data$sampling_prob[which(data$arm == y)]
                                     imputed_Y <- numeric(n)
                                     imputed_Y[which(A == y)] <- residual_vec[which(A == y)] / prob_to_use
                                     sum((imputed_Y)^2) / n^2
                                   })
  }else{

    # extract the covariate, reward, arm and sampling probability
    X <- data$covariate
    Y <- data$reward
    A <- data$arm
    n <- length(A)

    # construct regression estimator
    reg_result <- sapply(arm_list, function(x){
      arm_id <- which(A == x)
      data_to_reg <- lapply(data, function(x) as.matrix(x)[arm_id, ])
      switch (reg_type,
              lasso = {
                fit_lasso(data = data_to_reg, hyperparams = hyperparams)
              },
              MLE = {
                fit_glm(data = data_to_reg, hyperparams = hyperparams)
              }
      )
    })

    # impute conditional mean
    conditional_mean_imputed <-   switch (hyperparams$family,
                                          gaussian = {
                                            sapply(arm_list,
                                                   function(x){
                                                     coef_fit <- reg_result["coefficient", x]$coefficient
                                                     if(is.null(X)){
                                                       rep(as.vector(coef_fit), n)
                                                     }else{
                                                       as.vector(cbind(1, X) %*% coef_fit)
                                                     }
                                                   })
                                          }
    )

    # G-estimation term
    g_estimate <- apply(conditional_mean_imputed, 2, mean)

    # augmented term
    residual_vec <- sapply(1:length(A), function(x){
      Y[x] - conditional_mean_imputed[x, A[x]]
    })
    data_aug <- data
    data_aug$reward <- residual_vec
    augmented_term <- IPW(data = data_aug)["point", ]

    # store the mean estimate
    output["point", ] <- augmented_term + g_estimate

    # compute the variance term
    output["variance", ] <- sapply(arm_list,
                                   function(y){
                                     prob_to_use <- data$sampling_prob[which(data$arm == y)]
                                     imputed_Y <- numeric(n)
                                     imputed_Y[which(A == y)] <- residual_vec[which(A == y)] / prob_to_use
                                     mean((imputed_Y)^2)
                                   })
  }

  if(return_reg){
    # return the final value together with regression estimate
    return(list(output = output,
                coef_estimate = sapply(arm_list, function(x){
                  as.vector(reg_result["coefficient", x]$coefficient)
                })))
  }else{
    # return the final value
    return(output)
  }
}


#' This is a function for IPW point estimate
#'
#' @param data Data containing reward, arm, batch_id and sampling probability
#' @param one_batch A logical variable deciding if only one of the batches is only used
#' @param batch_output Which batch result will be output
#'
#' @return A vector of estimated mean and variance for each arm
#' @export
IPW <- function(data, one_batch = FALSE, batch_output = 1){

  # extract the data list
  Y <- as.vector(data$reward)
  A <- as.vector(data$arm)
  n <- length(A)
  sampling_prob <- data$sampling_prob
  batch_id <- data$batch_id
  batch_list <- unique(batch_id)
  arm_list <- sort(unique(A))

  # construct prob_tibble
  prob_tibble <- as_tibble(cbind(A, batch_id, sampling_prob)) |>
    dplyr::distinct()

  # construct an empty matrix
  output <- matrix(data = NA,
                   nrow = length(arm_list),
                   ncol = 2,
                   dimnames = list(estimate = c("point", "variance"),
                                   arm = arm_list))

  # compute the IPW estimator
  IPW_list <- sapply(batch_list, function(x){
    sapply(arm_list,
           function(y){
             prob_to_use <- prob_tibble |>
               filter(batch_id == x & A == y) |>
               dplyr::select(sampling_prob) |>
               pull()
             sum(Y[which(A == y & batch_id == x)]) / (n * prob_to_use)
           })
  })

  # depend on if only first batch is used
  if(one_batch){

    # compute the corresponding IPW with batch_output
    output["point", ] <- 2 * IPW_list[, batch_output]
    output["variance", ] <- sapply(arm_list,
                                   function(y){
                                     prob_to_use <- data$sampling_prob[which((data$arm == y) & (data$batch_id == batch_output))]
                                     imputed_Y <- numeric(n / 2)
                                     imputed_Y[which((A == y) & (batch_id == batch_output)) - (n/2)*(batch_output - 1)] <- Y[which((A == y) & (batch_id == batch_output))] / prob_to_use
                                     mean((imputed_Y  - output["point", y])^2)
                                   })

    # return the value
    return(output)
  }else{

    # store the point estimate
    output["point", ] <- apply(IPW_list, 1, sum)

    # compute the variance estimator
    output["variance", ] <- sapply(arm_list,
                                   function(y){
                                     prob_to_use <- data$sampling_prob[which(data$arm == y)]
                                     imputed_Y <- numeric(n)
                                     imputed_Y[which(A == y)] <- Y[which(A == y)] / prob_to_use
                                     sum((imputed_Y  - output["point", y])^2) / n^2
                                   })
  }
  return(output)
}

#' Tnis is a function building up the point and variance estimates wight weighting strategy
#'
#' @param data Same as AIPW function
#' @param hyperparams Same as AIPW function
#' @param reg_type Same as AIPW function
#'
#' @return An output matrix with point and variance estimate for each arm
#' @export
weighted_AIPW <- function(data, hyperparams, reg_type = "lasso"){

  # extract the covariate, reward, arm and sampling probability
  X <- data$covariate
  Y <- as.vector(data$reward)
  A <- as.vector(data$arm)
  n <- length(A)
  sampling_prob <- data$sampling_prob
  batch_id <- data$batch_id
  batch_list <- unique(batch_id)
  arm_list <- sort(unique(A))

  # construct regression estimator
  reg_result <- sapply(arm_list, function(x){
    arm_id <- which(A == x)
    data_to_reg <- lapply(data, function(x) as.matrix(x)[arm_id, ])
    switch (reg_type,
            lasso = {
              fit_lasso(data = data_to_reg, hyperparams = hyperparams)
            },
            MLE = {
              fit_glm(data = data_to_reg, hyperparams = hyperparams)
            }
    )
  })

  # construct prob_tibble
  prob_tibble <- as_tibble(cbind(A, batch_id, sampling_prob)) |>
    dplyr::distinct()

  # construct the sampling prob mat
  sampling_prob_mat <- t(
    sapply(1:length(A), function(x){
      as.numeric(
        sapply(arm_list, function(y){
          batch_id_rename <- batch_id
          prob_tibble |>
            filter((A == y) & (batch_id == batch_id_rename[x])) |>
            dplyr::select(sampling_prob) |> pull()
        })
      )
    })
  )

  # construct augmented term
  if(is.null(X)){

    # compute the mean via IPW
    conditional_mean_imputed <- matrix(c(rep(2, n), rep(6, n)),
                                       nrow = n,
                                       byrow = FALSE)

    augmented_term <- t(
      sapply(1:length(A), function(x){
        residual <- Y[x] - conditional_mean_imputed[x, A[x]]
        sapply(arm_list, function(y){
          dplyr::if_else(y == A[x], 1, 0) * residual / sampling_prob_mat[x, y]
        })
      })
    )
  }else{
    # construct g-estimate term by imputing the conditional mean
    conditional_mean_imputed <-   switch (hyperparams$family,
                                          gaussian = {
                                            sapply(arm_list,
                                                   function(x){
                                                     coef_fit <- reg_result["coefficient", x]$coefficient
                                                     as.vector(cbind(1, X) %*% coef_fit)
                                                   })
                                          }
    )

    # compute the augmented term
    augmented_term <- t(
      sapply(1:length(A), function(x){
        residual <- Y[x] - conditional_mean_imputed[x, A[x]]
        sapply(arm_list, function(y){
          dplyr::if_else(y == A[x], 1, 0) * residual / sampling_prob_mat[x, y]
        })
      })
    )
  }

  # final matrix
  Gamma <- conditional_mean_imputed + augmented_term

  # constant allocation
  weight_constant <- sqrt(sampling_prob_mat / length(batch_id))

  # compute the sum of weights
  sum_weight <- apply(weight_constant, 2, sum)

  # construct the point estimate
  weighted_AIPW <- apply(Gamma * weight_constant, 2, sum) / sum_weight

  # compute the difference of the matrix
  diff_mat <- Gamma - matrix(weighted_AIPW, length(A), length(arm_list), byrow = TRUE)

  # construct variance estimator
  variance_estimate <-  apply((diff_mat*weight_constant)^2, 2, sum) / sum_weight^2

  # construct an empty matrix
  output <- matrix(data = NA,
                   nrow = length(arm_list),
                   ncol = 2,
                   dimnames = list(estimate = c("point", "variance"),
                                   arm = arm_list))

  # impute the final output
  output["point", ] <- weighted_AIPW
  output["variance", ] <- variance_estimate

  # return the final confidence interval
  return(output)
}

#' Tnis is a function building up the point and variance estimates wight weighting strategy
#'
#' @param data Same as AIPW function
#' @param weighting Weighting method including aw, dm
#'
#' @return An output matrix with point and variance estimate for each arm
#' @export
weighted_IPW <- function(data, weighting = "aw") {
  # Extract the reward, arm, sampling probability, and batch id
  Y <- data$reward
  A <- data$arm
  sampling_prob <- data$sampling_prob
  batch_id <- data$batch_id
  arm_list <- sort(unique(A))

  # Construct prob_tibble using distinct values
  prob_tibble <- as_tibble(cbind(A, batch_id, sampling_prob)) %>%
    dplyr::distinct()

  # Construct the sampling probability matrix using vectorized operations
  combined_data <- data.frame(A = A, sampling_prob = sampling_prob)
  sampling_prob_mat <- combined_data |>
    dplyr::mutate(arm_1 = if_else(A == 1, sampling_prob, 1 - sampling_prob),
                  arm_2 = if_else(A == 2, sampling_prob, 1 - sampling_prob)) |>
    dplyr::select(arm_1, arm_2)

  # Construct augmented term using vectorized operations
  augmented_term <- matrix(0, nrow = length(A), ncol = length(arm_list))
  for (i in seq_along(arm_list)) {
    augmented_term[, i] <- ifelse(A == arm_list[i], Y / sampling_prob_mat[, i], 0)
  }

  # Constant allocation
  if(weighting == "aw"){
    weight_constant <- sqrt(sampling_prob_mat / length(batch_id))
  }else{
    weight_constant <- sampling_prob_mat / length(batch_id)
  }

  # Compute the sum of weights
  sum_weight <- colSums(weight_constant)

  # Construct the point estimate
  weighted_IPW <- colSums(augmented_term * weight_constant) / sum_weight

  # Compute the difference of the matrix
  diff_mat <- augmented_term - matrix(weighted_IPW, nrow = length(A), ncol = length(arm_list), byrow = TRUE)

  # Construct variance estimator
  variance_estimate <- colSums((diff_mat * weight_constant) ^ 2) / sum_weight ^ 2

  # Construct an empty matrix for the final output
  output <- matrix(NA, nrow = 2, ncol = length(arm_list), dimnames = list(c("point", "variance"), arm_list))

  # Impute the final output
  output["point", ] <- weighted_IPW
  output["variance", ] <- variance_estimate

  # Return the final confidence interval
  return(output)
}

#' A function providing confidence interval for weighted_AIPW
#'
#' @param data Data as input same as in CI_AIPW function
#' @param hyperparams Same as CI_AIPW function
#' @param reg_type Same as CI_AIPW function
#' @param estimand Same as CI_AIPW function
#' @param arm Same as CI_AIPW function
#' @param sig_level Same as CI_AIPW function
#'
#' @return Same as CI_AIPW function
#' @export
#' @importFrom stats qnorm
CI_weighted_AIPW <- function(data, hyperparams, reg_type, estimand = "ATE",
                             arm, sig_level = 0.05){

  weight_AIPW_result <-   weighted_AIPW(data = data,
                                        hyperparams = hyperparams,
                                        reg_type = reg_type)
  switch (estimand,
          ATE = {

            # chekc if a pair of arm is specified or not
            if(length(arm) != 2){
              stop("A pair of arm should be provided to compute the confidence interval!")
            }

            # compute the variance estimate
            variance_estimate <- sqrt(sum(weight_AIPW_result["variance", arm]))
            point_estimate <- weight_AIPW_result["point", arm[1]] - weight_AIPW_result["point", arm[2]]

            # compute the lower and upper confidence interval
            confidence_interval_lower <- point_estimate + qnorm(sig_level / 2)*variance_estimate
            confidence_interval_upper <- point_estimate + qnorm(1 - sig_level / 2)*variance_estimate
          },
          arm = {

            # check if arm is specified or not
            if(is.null(arm)){
              stop("Specific arm index should be provided to compute the confidence interval!")
            }

            # compute the variance estimate
            variance_estimate <- sqrt(sum(weight_AIPW_result["variance", arm]))
            point_estimate <- weight_AIPW_result["point", arm]

            # compute the lower and upper confidence interval
            confidence_interval_lower <- point_estimate + qnorm(sig_level / 2)*variance_estimate
            confidence_interval_upper <- point_estimate + qnorm(1 - sig_level / 2)*variance_estimate
          }
  )
  return(c(confidence_interval_lower, confidence_interval_upper))
}

#' A function providing confidence interval for weighted_IPW
#'
#' @param data Data as input same as in CI_AIPW function
#' @param estimand Same as CI_AIPW function
#' @param arm Same as CI_AIPW function
#' @param sig_level Same as CI_AIPW function
#'
#' @return Same as CI_AIPW function
#' @export
CI_weighted_IPW <- function(data, estimand = "ATE",
                            arm, sig_level = 0.05){

  weight_IPW_result <-   weighted_IPW(data = data)
  switch (estimand,
          ATE = {

            # chekc if a pair of arm is specified or not
            if(length(arm) != 2){
              stop("A pair of arm should be provided to compute the confidence interval!")
            }

            # compute the variance estimate
            sd_estimate <- sqrt(sum(weight_IPW_result["variance", arm]))
            point_estimate <- weight_IPW_result["point", arm[1]] - weight_IPW_result["point", arm[2]]

            # compute the lower and upper confidence interval
            confidence_interval_lower <- point_estimate + qnorm(sig_level / 2)*sd_estimate
            confidence_interval_upper <- point_estimate + qnorm(1 - sig_level / 2)*sd_estimate
          },
          arm = {

            # check if arm is specified or not
            if(is.null(arm)){
              stop("Specific arm index should be provided to compute the confidence interval!")
            }

            # compute the variance estimate
            sd_estimate <- sqrt(sum(weight_IPW_result["variance", arm]))
            point_estimate <- weight_IPW_result["point", arm]

            # compute the lower and upper confidence interval
            confidence_interval_lower <- point_estimate + qnorm(sig_level / 2)*sd_estimate
            confidence_interval_upper <- point_estimate + qnorm(1 - sig_level / 2)*sd_estimate
          }
  )
  return(list(
    confidence_interval = c(confidence_interval_lower, confidence_interval_upper),
    point_estimate = point_estimate,
    sd_estimate = sd_estimate
  ))
}

