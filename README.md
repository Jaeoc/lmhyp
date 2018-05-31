
<!-- README.md is generated from README.Rmd. Please edit that file -->
lmhyp
=====

This package provides an easy way to test hypotheses about continuous predictors in multiple regression. Simply fit a linear model in R with `lm`, specify your hypotheses and test them with `test_hyp`. The package enables formal testing of hypotheses using a default Bayesian prior and is particularly useful for testing multiple contradicting hypotheses. A tutorial paper for the package is forthcoming.

Basic example
-------------

``` r
library(lmhyp)

###Fit a linear model with lm
dt <- as.data.frame(scale(mtcars[, c(1, 3:4, 6)]))
fit <- lm(mpg ~ disp + hp + wt, data = dt)

###Define hypotheses and test them
hyp <- "wt > disp; wt < disp; wt = disp"
test_hyp(fit, hyp)
#> Hypotheses:
#> 
#>   H1:   "wt>disp"
#>   H2:   "wt<disp"
#>   H3:   "wt=disp"
#> 
#> Posterior probability of each hypothesis (rounded):
#> 
#>   H1:   0.0302
#>   H2:   0.5021
#>   H3:   0.4677
```

Installation
------------

You can install lmhyp from github with:

``` r
# install.packages("devtools")
devtools::install_github("Jaeoc/lmhyp")
```

Usage
-----

This function is based on a Bayes factor method by Mulder (2014). In essence, it uses a number of observations equal to the number of predictors in the linear model to construct a minimally informative prior, and the remainder of the observations are then used to test the hypotheses.

The function `test_hyp` takes as primary input a linear model object and a string vector specifying hypotheses. Simply specify the linear model as you normally would in R using `lm`, although for interpretation purposes it is preferable to standardize variables first.

``` r
fit <- lm(mpg ~ disp + hp + wt, data = dt)
```

Hypotheses are then specified using the variable names of the fitted object. Spaces do not matter. For example, if we think that both `wt` and `disp` have negative associations with the DV, this can be specified as:

``` r
Hyp <- "wt < 0, disp < 0"
```

When comparing several variables with the same value or variable they can be grouped with parentheses. The following specification is equivalent to the above:

``` r
Hyp <- "(wt, disp) < 0"
```

Parentheses become more important when defining hypotheses that involve more variables. For example, a hypothesis might be that the third variable `hp` also has a negative association with the DV, but that `wt` has a stronger negative association than either of the other variables. In that case it becomes necessary to use parentheses because the function does not permit repeating variables when specifying a hypothesis. This hypothesis would be then be specified as follows:

``` r
H1 <- "wt < (disp, hp) < 0"
```

If we have several contradicting hypotheses, as in the basic example, we can compare these directly by including them all in the same string vector separated with semicolons. For example, an alternative hypothesis might be that `wt` has a negative association, but that `disp` and `hp` have no association with the DV:

``` r
H2 <- "wt < (disp, hp) = 0"
```

Since H1 and H2 make predictions regarding the same variables we can test them directly against each other. Note the semicolons that separate each hypothesis:

``` r
H1v2 <- "wt < (disp, hp) < 0; wt < (disp, hp) = 0"
```

Testing the hypotheses is then as simple as inputting the `lm` object and hypothesis object into the `test_hyp` function:

``` r
result <- test_hyp(fit, H1v2)
```

When specified hypotheses are overlapping as in this case (in both hypotheses wt &lt; 0) the function takes slightly longer to finish since it runs 1e4 Monte Carlo iterations to check the extent of the overlap.

The function gives as output which hypotheses were specified and the posterior probability of each hypothesis. If the input hypotheses are not exhaustive (i.e., do not cover the full parameter space), the posterior probability of the complement is also output. The complement is the hypothesis that neither of the input hypotheses is true. In the current case we get the following output:

``` r
result
#> Hypotheses:
#> 
#>   H1:   "wt<(disp,hp)<0"
#>   H2:   "wt<(disp,hp)=0"
#>   Hc:   "Not H1-H2"
#> 
#> Posterior probability of each hypothesis (rounded):
#> 
#>   H1:   0.9915
#>   H2:   0.0027
#>   Hc:   0.0058
```

We see that there is strong evidence in the data for H1 and very weak evidence for H2 and the complement Hc. Exactly how the hypotheses compare to each other we can find out by printing the Bayes factor matrix:

``` r
result$BF_matrix
#>          H1          H2          Hc
#> H1   1.0000 0.002741001 0.005866392
#> H2 364.8302 1.000000000 2.140237177
#> Hc 170.4625 0.467237935 1.000000000
```

By definition, the Bayes Factor is interpreted as the likelihood of the data under a hypothesis compared to another. However, because the prior is the same for all hypotheses we can directly interpret the BFs as the likelihood of a hypothesis compared to another. Thus, given the data, H1 is 365 times as likely as H2 (BF = 364.83) and 170 times as likely as the complement (BF = 170.46).

References
----------

Mulder, J. (2014). Prior adjusted default Bayes factors for testing (in) equality constrained hypotheses. Computational Statistics & Data Analysis, 71, 448-463.
