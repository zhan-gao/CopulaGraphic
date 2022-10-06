data(tau, d, rho)
t_tr <- tau[as.logical(d)]
weibull_mle(t_tr)
ck_mle_weibull(tau, d, rho)

ck_mle_lognormal(tau, d, rho)



plot(density(t_tr))
lines(
    0:1500,
    dweibull(
        x = 0:1500,
        shape = 1.385965,
        scale = 468.052
    )
)
