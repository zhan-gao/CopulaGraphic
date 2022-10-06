#' Pseudo Model Method of Moments to choose the depedent parameter alpha
#'
#' @import nloptr
#'
#' @param tau min{T, C} without random censoring and Delta_P = 1
#' @param d I(T < C)
#' @param p disease progress
#' @param moment_option_vec
#' @param S = 20000
#' @param lb
#' @param ub
#' @param para_init
#' @param range_factor
#' @param seed
#' @param num_rand_perturb
#' @param two_stage
#' @param min_sig
#' @param x_tol_rel
#' @param maxeval
#'
#' @export

mm_alpha <- function(tau, d, p, moment_options_vec = 1:6, S = 20000,
                     lb = c(-10, 0, 0, -10, 0, 0, 0),
                     ub = c(10, 5, 5, 10, 5, 5, 1),
                     para_init = c(0.5, 1, 0.4, 0.5, 1, 0.4, 0.0055),
                     range_factor = c(10, 2, 2, 10, 2, 2, 500),
                     seed = 99,
                     num_rand_perturb = 100,
                     two_stage = FALSE,
                     min_sig = 0.25,
                     x_tol_rel = 1.0e-8,
                     maxeval = 1.0e4) {

    num_m <- length(moment_options_vec)


    theta_hat_1 <- matrix(0, 6, num_m)
    gamma_hat_1 <- rep(0, num_m)

    ktau_vec <- rep(0, num_m)
    for (i in 1:num_m) {
        moment_sel <- moment_comb(moment_options_vec[i])

        # Compute the sample moments
        m_sample <- sample_moments_log(tau, d, p, moment_sel)
        md_out <- md_est(m_sample, S, moment_sel, lb, ub, para_init, range_factor, seed,
                         num_rand_perturb, x_tol_rel, maxeval)
        para_hat <- md_out$para_hat
        if (two_stage) {
            para_init <-  md_out$para_hat
            if (para_init[3] < min_sig) {
                para_init[3] <- min_sig
            }
            if (para_init[6] < min_sig) {
                para_init[6] <- min_sig
            }

            md_out_2 <- md_est(
                m_sample,
                para_init,
                range_factor,
                seed,
                num_rand_perturb,
                x_tol_rel = x_tol_rel,
                maxeval = maxeval
            )

            if (md_out_2$obj < md_out$obj) {
                para_hat <- md_out_2$para_hat
            }
        }

        p <- rexp(S, para_hat[length(para_hat)])
        data_tc <- generate_tc_lognormal(para_hat[-length(para_hat)], p)
        t <- data_tc$t
        c <- data_tc$c

        ktau_vec[i] <- cor(t, c, method = "kendall")

    }

    return (ktau_vec)
}


#' Minimum Distance Estimation
#'
md_est <- function(m_sample,
                   S,
                   moment_sel,
                   lb,
                   ub,
                   para_init,
                   range_factor,
                   seed,
                   num_rand_perturb,
                   x_tol_rel = 1.0e-8,
                   maxeval = 1.0e4) {

    para_hat_mat <- matrix(0, num_rand_perturb, length(para_init))
    obj_value_vec <- rep(Inf, num_rand_perturb)
    set.seed(seed)



    for (j in 1:num_rand_perturb) {
        para_init_j <- para_init +
            (runif(length(para_init)) - 0.5) / range_factor
        t_0 <- Sys.time()
        opts <- list(
            "algorithm" = "NLOPT_LN_COBYLA",
            "xtol_rel" = x_tol_rel,
            "maxeval" = maxeval
        )
        nlopt_sol <- nloptr::nloptr(
            x0 = para_init_j,
            eval_f = min_dist_obj,
            lb = lb,
            ub = ub,
            opts = opts,
            S = S,
            m_sample = m_sample,
            moment_sel = moment_sel
        )
        para_hat_mat[j, ] <- nlopt_sol$solution
        obj_value_vec[j] <- nlopt_sol$objective
        t_1 <- Sys.time()
        print(difftime(t_1, t_0, units = "secs"))
    }
    j_min <- which.min(obj_value_vec)

    return(
        list(
            para_hat = para_hat_mat[j_min, ],
            obj = obj_value_vec[j_min]
        )
    )
}

#' moment combination indexes
#'
moment_comb <- function(opt) {
    if (opt == 1) {
        moment_sel = c(2:3, 6:7, 11, 12:13) # 7 moments
    } else if (opt == 2) {
        moment_sel = c(2:4, 7:9, 12:13); # 8 moments
    } else if (opt == 3) {
        moment_sel = c(2:3, 5, 7:8, 10, 12:13) # 8 moments
    } else if (opt == 4) {
        moment_sel = c(1:4, 6:9) # 8 moments
    } else if (opt == 5) {
        moment_sel = c(1:3, 5, 6:8, 10) # 8 moments
    } else if (opt == 6) {
        moment_sel = 1:10 # 10 moments
    }

    return(moment_sel)
}

#' Compute the sample moments of log variables
#'
#' @param tau min{T, C} without random censoring and Delta_P = 1
#' @param d I(T < C)
#' @param p disease progress
#' @moment_sel index of moments in use
#'
sample_moments_log <- function(tau, d, p, moment_sel = NULL) {

    log_P <- log(p)
    log_tau <- log(tau)
    Delta_T <- as.logical(d)
    Delta_C <- !Delta_T

    E_TP <- mean(log_tau[Delta_T] * log_P[Delta_T])
    E_T <- mean(log_tau[Delta_T])
    E_P_T <- mean(log_P[Delta_T]);
    E_P2_T <- mean(log_P[Delta_T]^2)
    E_T2 <- mean(log_tau[Delta_T]^2);
    frac_T_P <- mean(Delta_T)
    Corr_TP <- (E_TP - E_T * E_P_T) / (sqrt(E_T2 - E_T^2) * sqrt(E_P2_T - E_P_T^2))

    E_CP <- mean(log_tau[Delta_C] * log_P[Delta_C])
    E_C <- mean(log_tau[Delta_C])
    E_P_C <- mean(log_P[Delta_C])
    E_P2_C <- mean(log_P[Delta_C]^2)
    E_C2 <- mean(log_tau[Delta_C]^2)
    Corr_CP <- (E_CP - E_C * E_P_C) / (sqrt(E_C2 - E_C^2) * sqrt(E_P2_C - E_P_C^2))

    m_sample <- c(E_TP, E_T, E_P_T, E_P2_T, E_T2,
                  E_CP, E_C, E_P_C, E_P2_C, E_C2,
                  frac_T_P, Corr_TP, Corr_CP)
    if (is.null(moment_sel)) {
        m_out <- m_sample
    } else {
        m_out <- m_sample[moment_sel]
    }
}

#' Compute the moments of log variables based on the pseudo model
#'
#' @param theta parameters of the log-normal distributions of T and C
#' @param gamma parameter of the exponential distribution of P
#' @param S number of sample draws
#' @param moment_sel index of moments in use
#'
eval_moment_log <- function(theta, gamma, S, moment_sel = NULL, seed = NULL) {

    if(!is.null(seed)) {set.seed(seed)}

    p <- rexp(S, gamma)
    data_tc <- generate_tc_lognormal(theta, p)
    t <- data_tc$t
    c <- data_tc$c

    Delta_TP <- as.logical((t <= c) * (p <= t))
    Dleta_CP <- as.logical((c <= t) * (p <= c))
    Delta_P <- as.logical((p <= c) * (p <= t))


    # Take only Delta_P == 1
    tau <- apply(cbind(t, c), 1, min)
    tau <- tau[Delta_P]
    p_obs <- p[Delta_P]
    d <- as.logical((t <= c))
    d <- d[Delta_P]

    m_theta <- sample_moments_log(tau, d, p_obs, moment_sel)

    return(m_theta)
}

#' Objective function of minimum distance (Log normal T and C and log moments)
#'
#'

min_dist_obj <- function(para, S, m_sample, moment_sel) {
    theta <- para[-length(para)]
    gamma <- para[length(para)]
    m_theta <- eval_moment_log(theta, gamma, S, moment_sel)
    sum((m_sample - m_theta)^2)
}


#' Generate T and C based log-normal parameters and p
#'
#' @param theta
#' @param p
#'
generate_tc_lognormal <- function(theta, p, min_sig = 0.0001) {

    S = length(p)

    phi_0T <- theta[1]
    phi_1T <- theta[2]
    sigma_T <- max(theta[3], min_sig)
    phi_0C <- theta[4]
    phi_1C <- theta[5]
    sigma_C <- max(theta[6], min_sig)

    tc_mat <- sapply(p, function(x){
        c(rlnorm(1, phi_0T + phi_1T * log(x), sigma_T),
          rlnorm(1, phi_0C + phi_1C * log(x), sigma_C))
    })


    return(data.frame(t = tc_mat[1, ],
                     c = tc_mat[2, ]))

}
