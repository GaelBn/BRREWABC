
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BRREWABC <img src="man/figures/icon.png" align="right" width="90" />

<!-- badges: start -->
<!-- badges: end -->

BRREWABC (Batched Resilient and Rapid Estimation Workflow through
Approximate Bayesian Computation) : an R package designed to facilitate
inference through a parallelized Approximate Bayesian Computation
Sequential Monte Carlo (ABC SMC) algorithm. This package streamlines the
process of conducting Bayesian inference for complex models by
implementing efficient parallelization techniques.

## Overview

The algorithm used corresponds to an Approximate Bayesian Computation
approach using a Sequential Monte Carlo sampler. This iterative
algorithm enhances the basic ABC algorithm by incorporating two main
steps: weighted resampling of simulated particles and a gradual
reduction in tolerance. Similar to the ABC rejection approach, a prior
distribution is defined, aiming to estimate a posterior distribution. In
ABC-SMC, this estimation is achieved sequentially by constructing
intermediate distributions in each iteration, converging towards the
posterior distribution.

The specific implementation used in the package improves upon Del Moral
et al. (2006) original algorithm[^1] in three ways: - (1) an adaptive
threshold schedule selection based on quantiles of distances between
simulated and observed data [^2] [^3] - (2) an adaptive perturbation
kernel width during the sampling step, dependent on the previous
intermediate posterior distribution [^4] [^5] - and (3) the capability
to use multiple criteria simultaneously.

The algorithm used corresponds to an Approximate Bayesian Computation
approach using a Sequential Monte Carlo sampler (Del Moral et al.,
2006). This iterative algorithm enhances the basic ABC algorithm by
incorporating two main steps: weighted resampling of simulated particles
and a gradual reduction in tolerance. Similar to the ABC rejection
approach, a prior distribution is defined, aiming to estimate a
posterior distribution. In ABC-SMC, this estimation is achieved
sequentially by constructing intermediate distributions in each
iteration, converging towards the posterior distribution.

The specific implementation used in the package improves upon Del Moral
et al. (2006) original algorithm in three ways: - (1) an adaptive
threshold schedule selection based on quantiles of distances between
simulated and observed data (Del Moral et al., 2011; Drovandi and
Pettit, 2011) - (2) an adaptive perturbation kernel width during the
sampling step, dependent on the previous intermediate posterior
distribution (Beaumont et al., 2009; Toni et al., 2009), - and (3) the
capability to use multiple criteria simultaneously.

Del Moral, P., Doucet, A., and Jasra, A. (2006). Sequential Monte Carlo
samplers. Journal of the Royal Statistical Society : Series B
(Statistical Methodology), 68 :411–436.

Del Moral, P., A. Doucet and A. Jasra (2011): “An adaptive sequential
monte carlo method for approximate bayesian computation,” Stat. Comput.,
22, 1–12.

Drovandi, C. C. and A. N. Pettitt (2011): “Estimation of parameters for
macroparasite population evolution using approximate bayesian
computation” Biometrics, 67(1), 225–233.

Beaumont, M. A., Cornuet, J.-M., Marin, J.-M., and Robert, C. P. (2009).
Adaptive approximate Bayesian computation. Biometrika, 96(4) :983–990.

Toni, T., Welch, D., Strelkowa, N., Ipsen, A., and Stumpf, M. P. H.
(2009). Approximate Bayesian computation scheme for parameter inference
and model selection in dynamical systems. Journal of The Royal Society
Interface, 6(31) :187–202.

## Features

- **Parallelized ABC SMC**: Conduct inference using an Approximate
  Bayesian Computation Sequential Monte Carlo algorithm, parallelized
  for enhanced computational efficiency.
- **Flexible Model Specification**: Easily specify complex models.
- **Customizable Settings**: Fine-tune algorithm parameters to suit
  specific modeling needs and computational resources.
- **Scalable**: Utilize parallel computing capabilities to handle large
  datasets and complex models with ease.
- **Comprehensive Documentation**: Detailed documentation and examples
  to guide users through package functionality and usage.

## Installation

You can install the development version of BRREWABC from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("GaelBn/BRREWABC")
```

## Usage

<!-- For a basic example, see the [Get Started](`vignette("BRREWABC")`) section. -->

For a basic example, see the [Get
Started](https://gaelbn.github.io/BRREWABC/articles/BRREWABC.html)
section.

## Getting help

``` r
# Access package documentation
help(package = "BRREWABC")
```

## Contributing

Contributions to BRREW-ABC are welcome! If you encounter any issues,
have feature requests, or would like to contribute enhancements, please
feel free to contact us.

## Project status

This project is constantly evolving, according to needs and suggestions,
at a pace that depends on the time that can be devoted to it.

## Authors and acknowledgment

BRREW-ABC was developed by Gaël Beaunée thanks to the advice of a number
of people, many thanks to them.

## Troubleshooting

If you encounter any issues, please feel free to contact us.

[^1]: Del Moral, P., Doucet, A., and Jasra, A. (2006). Sequential Monte
    Carlo samplers. Journal of the Royal Statistical Society : Series B
    (Statistical Methodology), 68 :411–436

[^2]: Del Moral, P., A. Doucet and A. Jasra (2011): “An adaptive
    sequential monte carlo method for approximate bayesian computation,”
    Stat. Comput., 22, 1–12.

[^3]: Drovandi, C. C. and A. N. Pettitt (2011): “Estimation of
    parameters for macroparasite population evolution using approximate
    bayesian computation” Biometrics, 67(1), 225–233.

[^4]: Beaumont, M. A., Cornuet, J.-M., Marin, J.-M., and Robert, C. P.
    (2009). Adaptive approximate Bayesian computation. Biometrika, 96(4)
    :983–990.

[^5]: Toni, T., Welch, D., Strelkowa, N., Ipsen, A., and Stumpf, M. P.
    H. (2009). Approximate Bayesian computation scheme for parameter
    inference and model selection in dynamical systems. Journal of The
    Royal Society Interface, 6(31) :187–202.
