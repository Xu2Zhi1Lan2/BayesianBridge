# BBHD: Bayesian Bridge with High Dimensions in R

## Overview

The `BBHD` package provides a robust framework for Bayesian Bridge Regression, leveraging advanced Bayesian statistical methods to deliver both shrinkage and selection capabilities in regression analysis. This package includes specialized functions for fitting Bayesian models under bridge penalties, offering efficient algorithms tailored for high-dimensional data.

## Installation

To install the latest version of the BBHD package, use the following command in R:

```r
# Install directly from GitHub repository
# install.packages("devtools")
devtools::install_github("https://github.com/Xu2Zhi1Lan2/BayesianBridge")
```

## Package Structure

The `BBHD` package mainly includes:

- `R/`: Contains the core implementation of the package functions.
  - `BBLR.R`: Implements the Bayesian Bridge Linear Regression.
  - `BBNM.R`: Implements a Bayesian Bridge Nonlinear Model.
- `vignettes/`: Includes a detailed tutorial on Bayesian Bridge Regression (`Bayesian-Bridge-tutorial.Rmd`).

## Key Features and Functions

### Functions

- **`BBLR()`**: Performs Bayesian Bridge Linear Regression.
- **`BBNM()`**: Extends the Bayesian Bridge Normal Mean.

### Key Parameters

- `y`: Response variable vector.
- `X`: Predictor matrix.
- `step_sizes_tuning_iterations`: The number of iterations dedicated to tuning the MCMC sampling algorithm for optimal step sizes.
- `burn`: Burn-in period for the sampling algorithm.
- `nmc`: Number of Monte Carlo samples to draw post burn-in period.
- `thin`: Thinning factor for the samples.
- `method.alpha`: Method for updating the `alpha` parameter during sampling.

Both functions utilize advanced MCMC techniques, specifically Metropolis-Hastings within Gibbs sampling, to efficiently explore the posterior distributions of the model parameters in high dimensions case.

## Examples

You can find practical examples in the vignettes folder. To access the vignette after installation, you can run:

```r
vignette("Bayesian-Bridge-tutorial", package = "BBHD")
```

This vignette includes a comprehensive guide on using `BBHD` for Bayesian Bridge Regression, including setting up data, running the models.

## Contributing

Contributions to improve `BBHD` or extend its capabilities are welcome. Please submit pull requests or report issues on the repository page.

## References

[1] Nicholas G. Polson, James G. Scott, Jesse Windle, The Bayesian Bridge, *Journal of the Royal Statistical Society Series B: Statistical Methodology*, Volume 76, Issue 4, September 2014, Pages 713–733, https://doi.org/10.1111/rssb.12042

[2] Mallick, H., & Yi, N. (2018). Bayesian Bridge Regression. *Journal of applied statistics*, *45*(6), 988–1008. https://doi.org/10.1080/02664763.2017.1324565