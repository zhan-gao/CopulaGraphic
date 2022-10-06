#' Joe copula generator
#'
#'
#' @export
joe_phi <- function(u, a) {
    -log(1 - (1 - u)^a)
}

#'First order derivative of Joe copula generator
joe_phi_deriv <- function(u, a) {
    - (a * (1 - u)^(a - 1)) / (1 - (1 - u)^a)
}

#' Joe copula generator inverse
#'
#'
#'@export
joe_phi_inverse <- function(u, a) {
    1 - (1 - exp(-u))^(1/a)
}


#' Clayton Copula generator
#'
#' @export
#'
clayton_phi <- function(u, a) {
    (u ^ (-a) - 1) / a
}

#'First order derivative of clayton copula generator
clayton_phi_deriv <- function(u, a) {
    -u^(-a - 1)
}

#' Clayton Copula generator inverse
#'
#' @export
#'
clayton_phi_inverse <- function(u, a) {
    (1 + a * u) ^ (-1 / a)
}


#' Convert Kendall's tau to the dependence parameter alpha (Joe)
#'
#' @export
#'
joe_k_tau_to_alpha <- function(k_tau) {

    data("kendall_alpha_joe", envir=environment())
    ind_alpha <- sapply(k_tau, function(k) {
        which.min(abs(k - kendall_alpha_joe$ktau))
    })
    return(kendall_alpha_joe$alpha[ind_alpha])

}


#' Convert dependence parameter alpha (Joe) to Kendall's tau
#'
#' @export
#'
joe_alpha_to_k_tau <- function(alpha) {

    data("kendall_alpha_joe", envir=environment())
    ind_ktau <- sapply(alpha, function(a) {
        which.min(abs(a - kendall_alpha_joe$alpha))
    })
    return(kendall_alpha_joe$ktau[ind_ktau])
}



#' Convert Kendall's tau to the dependence parameter alpha (Clayton)
#'
#' @export
#'
clayton_k_tau_to_alpha <- function(k_tau) {
    (2 * k_tau) /(1 - k_tau)
}


archimedean_copula <- function(u, v, phi, phi_inv, a) {

    phi_inv(phi(u, a) + phi(v, a), a)

}


partial_archimedean_copula <- function(u, v, phi, phi_inv, phi_deriv, a) {

    num <- cbind(phi_deriv(u, a), phi_deriv(v, a))
    denom <- phi_deriv(
        phi_inv(phi(u, a) + phi(v, a), a), a
    )

    return(num / denom)

}

