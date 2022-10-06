# CopulaGraphic
A package for "*Copula Graphic Estimation of Survival Function with Dependent Censoring and its Application to an Analysis of Pancreatic Cancer Clinical Trial*" by Jung Hyun Jo, Zhan Gao, Inkyung Jung, Si Young Song, Geert Ridder and Hyungsik Roger Moon. 

### Installation

```R
devtools::install_github("zhan-gao/CopulaGraphic")
```

### Example

```R
library(CopulaGraphic)
data(tau, d, rho)
rw_out <- rw(
    tau,
    d,
    rho,
    alpha = 3.83,
    copula = "joe",
    quant_point = c(0.25, 0.5, 0.75),
    num_b = 10000,
    seed = 100
)
```

```R
# Estimated quantiles of survival times
print(rw_out$quant_time)
```

```
[1] 845.8924 366.0999 102.8447

```

```R
surv_plot(tau[rho == 1], rw_out$surv)
```

![surv_plot](surv_plot.png)
