#' Copula Graphic Estimator
#'
#' Estimate the survival function under both dependent and independent censoring
#' by a modified Rivest-Well (2001) estimator
#'
#' @inheritParams rw_est
#' @param t_eval
#' @param quant_point
#' @param num_b number of Bootstrap sampling
#' @param seed
#'
#'
#' @return A list
#' \item{surv}{Estimated survival function}
#' \item{quant_time}{Estimated quantiles of survival times}
#' \item{surv_eval}{Estimated survival probability at time _eval}
#' \item{quant_time_bt}{Bootsrapped estimated quantiles}
#' \item{surv_eval_bt}{Bootstrapped estimated survival probability at time_eval}
#'
#' @export
#'
rw <- function(tau, d, rho, alpha, copula = "joe", t_eval = NULL,
               quant_point = c(0.25, 0.5, 0.75), num_b = 10000, seed = NULL) {

    # Estimation
    rw_out <- rw_est(tau, d, rho, alpha, copula)
    quant_time_vec <- quantile_surv_time(rw_out$surv, rw_out$time, quant_point)
    surv_eval_vec <- surv_func_eval(rw_out$surv, rw_out$time, t_eval)

    # Bootstrap
    rw_bt_out <- rw_bootstrap(tau, d, rho, alpha, copula,
                              t_eval, quant_point, num_b, seed)
    quant_time_bt_mat <- rw_bt_out$quantile
    surv_eval_bt_mat <- rw_bt_out$surv_eval

    return(list(
        surv = rw_out$surv,
        quant_time = quant_time_vec,
        surv_eval = surv_eval_vec,
        quant_time_bt = quant_time_bt_mat,
        surv_eval_bt = surv_eval_bt_mat
    ))
}

#' Copula Graphic Estimator
#'
#' Modified Rivest and Wells (2001)'s estimator
#' using data (t_i, d_i) of individuals with rho_i = 1 to estimate the
#' survival function of the event at points t_(1), t_(2), ..., t_(n_1)
#'
#' @import survival
#'
#' @param tau min(T, C, Z)
#' @param d I(T < C); it is observed only when rho = 1
#' @param rho indicator for not censored by Z
#' @param alpha Joe Copula parameter (scalar >= 1)
#' @param copula c("joe", "clayton")
#'
#' @return A data.frame
#' \item{surv}{n_1-by-1 vector estimated survival function}
#' \item{time}{survival time}
#'
#' @export
#'
rw_est <- function(tau, d, rho, alpha, copula = "joe") {

    if (copula == "joe") {
        phi = joe_phi
        phi_inverse = joe_phi_inverse
    } else if (copula == "clayton") {
        phi = clayton_phi
        phi_inverse = clayton_phi_inverse
    } else {
        stop("WRONG COPULA INPUT")
    }

    # Keep noncensored individuals (rho_i = 1) in estimation
    index_rho <- as.logical(rho)
    t <- tau[index_rho]
    d <- d[index_rho]

    # Order t
    surv_time <- sort(t)
    t_order <- order(t)
    d_est <- d[t_order]

    # Kaplan-Meier estimator of the survival function of T at tau(rho = 1)
    km_fit <- km_est(tau, rho)
    pi_t <- km_fit$surv[-1]
    # # resolve the tied survival times
    pi_t <- rep(pi_t, as.numeric(table(t)))

    # km_fit <- survival::survfit(survival::Surv(tau, rho) ~ 1)

    # Rivst Wells Estimation
    len_pi <- length(pi_t)
    phi_inv_input <- -cumsum(
        (phi(pi_t[1:(len_pi - 1)], alpha) - phi(pi_t[2:len_pi], alpha)) *
            d_est[2:len_pi]
    )
    phi_inv_input[is.nan(phi_inv_input)] <- Inf
    s_y <- c(1, phi_inverse(phi_inv_input, alpha))

    return(data.frame(
        surv = s_y,
        time = surv_time
    ))

}


#' Bootstrap Copula Graphic Estimator
#'
#' @inheritParams rw
#'
#' @import doParallel
#' @import foreach
#' @import doRNG
#'
#' @return A list
#' \item{quantile}{A matrix of bootstrapped quantiles}
#' \item{surv_eval}{A matrix of bootstrapped survival fucntion evaluated at t_eval}
#' @export
#'

rw_bootstrap <- function(tau, d, rho, alpha, copula = "joe", t_eval = NULL,
                         quant_point = c(0.25, 0.5, 0.75),
                         num_b = 10000, seed = NULL) {
    n <- length(tau)
    num_q <- length(quant_point)

    if(!foreach::getDoParRegistered()) {
        num_cores <- parallel::detectCores()
        cltr <- parallel::makeCluster(num_cores)
        doParallel::registerDoParallel(cltr, cores = num_cores)
    }

    if(! is.null(seed)) {set.seed(seed)}
    bt_out <- foreach(
        b = 1:num_b,
        .combine = rbind,
        .packages = c("survival", "dplyr"),
        .export = c("rw_est", "quantile_surv_time", "surv_func_eval",
                    "km_est", "joe_phi", "joe_phi_inverse",
                    "clayton_phi", "clayton_phi_inverse")
    ) %dorng% {

        resample_id <- sample(n, n, replace = TRUE)
        tau_b <- tau[resample_id]
        d_b <- d[resample_id]
        rho_b <- rho[resample_id]

        rw_out_b <- rw_est(tau_b, d_b, rho_b, alpha, copula)
        c(quantile_surv_time(rw_out_b$surv, rw_out_b$time, quant_point),
          surv_func_eval(rw_out_b$surv, rw_out_b$time, t_eval))
    }

    list(quantile = bt_out[, 1:num_q],
         surv_eval = bt_out[, -(1:num_q)])
}
