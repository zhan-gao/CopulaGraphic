% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/copula_graphic.R
\name{rw}
\alias{rw}
\title{Copula Graphic Estimator}
\usage{
rw(
  tau,
  d,
  rho,
  alpha,
  copula = "joe",
  t_eval = NULL,
  quant_point = c(0.25, 0.5, 0.75),
  num_b = 10000,
  seed = NULL
)
}
\arguments{
\item{tau}{min(T, C, Z)}

\item{d}{I(T < C); it is observed only when rho = 1}

\item{rho}{indicator for not censored by Z}

\item{alpha}{Joe Copula parameter (scalar >= 1)}

\item{copula}{c("joe", "clayton")}

\item{num_b}{number of Bootstrap sampling}

\item{seed}{}
}
\value{
A list
\item{surv}{Estimated survival function}
\item{quant_time}{Estimated quantiles of survival times}
\item{surv_eval}{Estimated survival probability at time _eval}
\item{quant_time_bt}{Bootsrapped estimated quantiles}
\item{surv_eval_bt}{Bootstrapped estimated survival probability at time_eval}
}
\description{
Estimate the survival function under both dependent and independent censoring
by a modified Rivest-Well (2001) estimator
}
