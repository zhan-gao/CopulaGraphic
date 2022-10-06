


get_data_dp <- function() {
    data(tau, d, p, d_p, rho, envir=environment())
    ind_non_random_censor <- as.logical(rho)
    tau <- tau[ind_non_random_censor]
    d <- d[ind_non_random_censor]
    p <- p[ind_non_random_censor]
    ind_dp <- as.logical(d_p[ind_non_random_censor])
    tau <- tau[ind_dp]
    d <- d[ind_dp]
    p <- p[ind_dp]

    return(data.frame(
        tau = tau,
        d = d,
        p = p
    ))
}

#' Compute the quantiles of survival time
#'
#' Caveat: Note that the estimated survival function and time are sorted
#'
#'
#' @param surv
#' @param time
#' @param quant_point
#'
#' @export
#'
quantile_surv_time <- function(surv, time, quant_point = c(0.25, 0.5, 0.75)){
    n <- length(surv)
    num_q <- length(quant_point)
    q <- rep(0, num_q)

    for (i in 1:num_q){
        quant_i = quant_point[i]
        if (surv[n] >= quant_i) {
            q[i] <- time[n]
        } else {
            high_indx <- ((surv - quant_i) >= 0);
            low_indx <- ((quant_i - surv) > 0);
            q[i] <- (time[sum(high_indx)] + time[n - sum(low_indx) + 1]) / 2
        }
    }

    return(q)
}


#' Evaluate Survival Function
#'
#' @param surv_est
#' @param time_est
#' @param time_eval
#'
#'
#' @export
#'

surv_func_eval <- function(surv_est, time_est, time_eval) {

    if(is.null(time_eval)) return(NULL)

    surv_est <- c(1, surv_est)
    time_est <- c(0, time_est)

    ind <- sapply(time_eval, function(x) {
        sum(x >= time_est)
    })

    surv_eval <- surv_est[ind]

    return(surv_eval)
}
