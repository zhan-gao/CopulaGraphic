# Implement the MLE estimator in Czado Keilegom (2021)


#' MLE estimator in Czado and van Keilegom (2021) with bootstrap inference
#'
#' @export
#'
#'


ck_mle_weibull <- function(tau,
                           d,
                           rho,
                           copula = "joe",
                           x0 = NULL,
                           lb = NULL,
                           ub = NULL,
                           xtol_rel = 1e-8,
                           maxeval = 20000,
                           quant_point = c(0.75, 0.5, 0.25),
                           num_b = 10000,
                           seed = NULL) {

    #Estimation
    mle_out <- ck_mle_weibull_est(
        tau,
        d,
        rho,
        copula,
        x0 = x0,
        lb = lb,
        ub = ub,
        xtol_rel = xtol_rel,
        maxeval = maxeval
    )


    mle_bt_out <- ck_mle_weibull_bootstrap(tau,
                                           d,
                                           rho,
                                           copula,
                                           x0 = mle_out$para_hat,
                                           lb = lb,
                                           ub = ub,
                                           xtol_rel = xtol_rel,
                                           maxeval = maxeval,
                                           quant_point = quant_point,
                                           num_b = num_b,
                                           seed = seed)

    return(list(
        para_hat = mle_out$para_hat,
        quant_time = qweibull(quant_point, mle_out$para_hat[1], mle_out$para_hat[2]),
        para_hat_mat = mle_bt_out$para_hat_bt,
        alpha_hat = mle_bt_out$alpha_hat,
        quant_time_bt = mle_bt_out$quant
    ))

}




#' MLE Estimator in Czado Keilegom (2021) with weibull distributions
#'
#' @import nloptr
#'
#' @param tau min(T, C, Z)
#' @param d I(T < C); it is observed only when rho = 1
#' @param rho indicator for not censored by Z
#' @param copula c("joe", "clayton")
#' @param xtol_rel parameter for nloptr sovler
#' @param maxeval parameter for nloptr sovler
#'
#' @export


ck_mle_weibull_est <- function(tau,
                               d,
                               rho,
                               copula = "joe",
                               x0 = NULL,
                               lb = NULL,
                               ub = NULL,
                               xtol_rel = 1e-8,
                               maxeval = 20000,
                               quant_point = c(0.75, 0.5, 0.25)) {
    if (is.null(x0)) {
        para_t <- weibull_mle(tau[as.logical(d)])
        para_c <- weibull_mle(tau[!as.logical(d)])
        a <- joe_k_tau_to_alpha(0.3)
        x0 <- c(para_t, para_c, a)
    }

    if (is.null(lb)) {
        lb <- c(0, 0, 0, 0, 1)
    }

    if (is.null(ub)) {
        ub <- c(x0[-5] * 2, 15)
    }

    opts <- list(
        "algorithm" = "NLOPT_LN_COBYLA",
        "xtol_rel" = xtol_rel,
        "maxeval" = maxeval
    )

    nlopt_sol <- nloptr::nloptr(
        x0 = x0,
        eval_f = weibull_loglik,
        lb = lb,
        ub = ub,
        opts = opts,
        tau = tau,
        d = d,
        rho = rho,
        copula = copula
    )

    return(list(
        para_hat = nlopt_sol$solution,
        quant = qweibull(quant_point, nlopt_sol$solution[1], nlopt_sol$solution[2])
    ))
}

#' Likelihood function
weibull_loglik <- function(para, tau, d, rho, copula = "joe") {
    shape_T <- para[1]
    scale_T <- para[2]
    shape_C <- para[3]
    scale_C <- para[4]
    a <- para[5]

    if (copula == "joe") {
        phi <- joe_phi
        phi_inv <- joe_phi_inverse
        phi_deriv <- joe_phi_deriv
    } else if (copula == "copula") {
        phi <- clayton_phi
        phi_inv <- clayton_phi_inverse
        phi_deriv <- clayton_phi_deriv
    } else {
        stop("WRONG COPULA INPUT")
    }

    surv_t <- 1 - pweibull(tau, shape = shape_T, scale = scale_T)
    surv_c <- 1 - pweibull(tau, shape = shape_C, scale = scale_C)
    partial_copula <-
        log(partial_archimedean_copula(surv_t, surv_c, phi, phi_inv, phi_deriv, a))
    log_f_t <-
        dweibull(tau,
                 shape = shape_T,
                 scale = scale_T,
                 log = TRUE)
    log_f_c <-
        dweibull(tau,
                 shape = shape_C,
                 scale = scale_C,
                 log = TRUE)

    I_1 = sum(rho * d * (log_f_t + partial_copula[, 1]))
    I_2 = sum(rho * (1 - d) * (log_f_c + partial_copula[, 2]))
    I_3 = sum((1 - rho) * log(archimedean_copula(surv_t, surv_c, phi, phi_inv, a)))

    - (I_1 + I_2 + I_3)
}


#' Bootstrap Weibull MLE in Czado and van Keilegom (2021)
#'
#' @import doParallel
#' @import foreach
#' @import doRNG
#'
#' @export
ck_mle_weibull_bootstrap <- function (tau,
                                      d,
                                      rho,
                                      copula = "joe",
                                      x0 = NULL,
                                      lb = NULL,
                                      ub = NULL,
                                      xtol_rel = 1e-8,
                                      maxeval = 20000,
                                      quant_point = c(0.75, 0.5, 0.25),
                                      num_b = 10000,
                                      seed = NULL) {
    n <- length(tau)

    if (!foreach::getDoParRegistered()) {
        num_cores <- parallel::detectCores()
        cltr <- parallel::makeCluster(num_cores)
        doParallel::registerDoParallel(cltr, cores = num_cores)
    }

    if (!is.null(seed)) {
        set.seed(seed)
    }
    bt_out <- foreach(
        b = 1:num_b,
        .combine = rbind,
        .packages = c("nloptr"),
        .export = c(
            "ck_mle_weibull_est",
            "weibull_loglik",
            "joe_phi",
            "joe_phi_inverse",
            "joe_phi_deriv",
            "clayton_phi",
            "clayton_phi_inverse",
            "clayton_phi_deriv",
            "archimedean_copula",
            "partial_archimedean_copula"
        )
    ) %dorng% {
        resample_id <- sample(n, n, replace = TRUE)
        tau_b <- tau[resample_id]
        d_b <- d[resample_id]
        rho_b <- rho[resample_id]

        mle_out_b <- ck_mle_weibull_est(
            tau_b,
            d_b,
            rho_b,
            copula,
            x0 = x0,
            lb = lb,
            ub = ub,
            xtol_rel = xtol_rel,
            maxeval = maxeval
        )$para_hat
        mle_out_b
    }

    quant_mat <- t(apply(bt_out, 1, function(x) {
        qweibull(quant_point, x[1], x[2])
    }))

    return(list(
        para_hat_bt = bt_out[, -5],
        alpha_hat = bt_out[, 5],
        quant = quant_mat
    ))

}



#' Weibull MLE
#'
#' @export
#'
weibull_mle <- function(x,
                        x0 = c(1, mean(x)),
                        lb = c(0, 0),
                        ub = c(10, max(x)),
                        xtol_rel = 1e-8,
                        maxeval = 20000) {
    loglik <- function(parm, tt) {
        gamma <- parm[1]
        lambda <- parm[2]
        loglik <-
            sum(dweibull(
                tt,
                shape = gamma,
                scale = lambda,
                log = TRUE
            ))
        return(-loglik)
    }
    opts <- list(
        "algorithm" = "NLOPT_LN_COBYLA",
        "xtol_rel" = 1e-8,
        "maxeval" = 5000
    )
    nlopt_sol <- nloptr::nloptr(
        x0 = x0,
        eval_f = loglik,
        lb = lb,
        ub = ub,
        opts = opts,
        tt = x,
    )

    return(nlopt_sol$solution)
}
