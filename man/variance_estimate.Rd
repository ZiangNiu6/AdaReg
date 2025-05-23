% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/help.R
\name{variance_estimate}
\alias{variance_estimate}
\title{Variance estimate}
\usage{
variance_estimate(data, conditional_mean, mean_estimate, arm)
}
\arguments{
\item{data}{Data including covariate, arm and reward}

\item{conditional_mean}{Conditional mean estimate for this specific arm}

\item{mean_estimate}{The estimate for mean of potential outcome}

\item{arm}{The arm to be evaluated}
}
\value{
the variance estimate
}
\description{
Variance estimate
}
