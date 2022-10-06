#' Bootstrapped Kaplan-Meier
#'
#' @param tau min(T, C)
#' @param d I(T < C)
#'
#' @import survival
#' @import doParallel
#' @import foreach
#' @import doRNG
#'
#' @return A list
#' \item{quant_time}{Estimated quantiles of survival times}
#' \item{quant_time_bt}{Bootsrapped estimated quantiles}
#' \item{CI}{Confidence Intervals}
#'
#' @export

km <- function(tau, d, sig_level = 0.05,
               quant_point = c(0.25, 0.5, 0.75), num_b = 10000, seed = NULL) {

    # Estimation
    n <- length(tau)
    num_q <- length(quant_point)

    # Fit
    km_fit <- km_est(tau, d)
    quant_hat <- quantile_surv_time(km_fit$surv, km_fit$time, quant_point)

    # Bootstrap
    km_fit_c <- km_est(tau, 1 - d)
    if(! is.null(seed)) {set.seed(seed)}

    if(!foreach::getDoParRegistered()) {
        num_cores <- parallel::detectCores()
        cltr <- parallel::makeCluster(num_cores)
        doParallel::registerDoParallel(cltr, cores = num_cores)
    }

    bt_out <- foreach(
        b = 1:num_b,
        .combine = rbind,
        .packages = c("survival"),
        .export = c("inverse_cdf", "km_est",
                    "quantile_surv_time")
    ) %dorng% {

        t_b <- inverse_cdf(1 - km_fit$surv,
                           km_fit$time,
                           n)
        c_b <- inverse_cdf(1 - km_fit_c$surv,
                           km_fit_c$time,
                           n)
        tau_b <- apply(cbind(t_b, c_b), 1, min)
        d_b <- as.numeric(t_b < c_b)
        km_fit_b <- km_est(tau_b, d_b)
        quantile_surv_time(km_fit_b$surv, km_fit_b$time, quant_point)
    }

    CI <- apply(bt_out, 2, quantile, probs = c(sig_level / 2,
                                               1 - (sig_level / 2)))

    return(list(
        surv = km_fit$surv,
        time = km_fit$time,
        quant_time = quant_hat,
        quant_tiime_bt = bt_out,
        CI = CI
    ))

}


#' generate data by inverse a estimated CDF function
inverse_cdf <- function(fi, xi, n, seed = NULL, round_digits = 5) {

    # Generate N uniformly distributed samples between 0 and 1.
    if(!is.null(seed)) {set.seed(seed)}
    u <- runif(n)
    # Map these to the points on the empirical CDF.
    x <- approx(fi, xi, u, method = 'linear', rule = 2)$y

    return(round(x, digits = round_digits))
}

#' Kaplan-Meier Estimation
#'
#' Mimicing Matlab ecdf function.
#'
#' @param tau min(T, C)
#' @param d I(T < C)
#'
#' @import dplyr
#'
#' @return survival probability and survival time for those with d == 1 with
#'   repeated entries removed
#'
#' @export
#'
km_est <- function(tau, d) {


    # The survifit outputs surv for each unique time in tau
    surv_fit <- survival::survfit(survival::Surv(tau, d) ~ 1)
    # Add the repeated time entry back
    pi_t <- rep(surv_fit$surv, as.numeric(table(tau)))
    # Extract those with Delta_T = 1 with a leading 1
    pi_t <- c(1, pi_t[as.logical(d[order(tau)])])

    surv_time <- sort(tau[as.logical(d)])
    surv_time <- c(surv_time[1], surv_time)

    dplyr::distinct(
        data.frame(surv = pi_t, time = surv_time)
    )
}
