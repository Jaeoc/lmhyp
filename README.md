
<!-- README.md is generated from README.Rmd. Please edit that file -->
lmhyp
=====

This package provides an easy way to test hypotheses about coefficients in multiple regression. Simply fit a linear model in R with `lm`, specify your hypotheses and test them with `test_hyp`. The package enables formal testing of hypotheses using a default Bayesian prior and is particularly useful for testing multiple contradicting hypotheses. A tutorial paper for the package is forthcoming.

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
#>   H1:   0.030
#>   H2:   0.502
#>   H3:   0.468
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

The function `test_hyp` takes as primary input a linear model object and a string vector specifying hypotheses. Simply specify the linear model as you normally would in R using `lm`. If hypotheses involve multiple variables it is typically preferable to standardize the variables first to facilitate interpretation.

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

The function gives as primary output which hypotheses were specified and the posterior probability of each hypothesis. If the input hypotheses are not exhaustive (i.e., do not cover the full parameter space), the posterior probability of the complement is also output. The complement is the hypothesis that neither of the input hypotheses is true. In the current case we get the following output:

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
#>   H1:   0.991
#>   H2:   0.003
#>   Hc:   0.006
```

We see that there is strong evidence in the data for H1 and very weak evidence for H2 and the complement Hc. Exactly how the hypotheses compare to each other we can find out by printing the Bayes factor matrix:

``` r
result$BF_matrix
#>       H1      H2      Hc
#> H1 1.000 356.994 166.581
#> H2 0.003   1.000   0.467
#> Hc 0.006   2.143   1.000
```

By definition, the Bayes Factor is interpreted as the likelihood of the data under a hypothesis compared to another. However, because the prior by default is the same for all hypotheses we can directly interpret the BFs as the likelihood of a hypothesis compared to another. Thus, given the data, H1 is 357 times as likely as H2 (*B*<sub>12</sub> = 356.99) and 167 times as likely as the complement (*B*<sub>1*c*</sub> = 166.58).

It is also possible to specify one's own prior probabilites as a vector of numbers. These must then be specified for all hypotheses, including the complement if one exists. The function will throw an error if too few or too many probabilites have been specified. They can be specified both as probabilites, e.g, c(0.5, 0.2, 0.3) or as relative weights, e.g., c(5, 2, 3). In the latter case the function will normalize the input values. So for our example we might get:

``` r
test_hyp(fit, H1v2, c(5, 2, 3))
#> Warning in test_hyp(fit, H1v2, c(5, 2, 3)): priorprobs did not sum to 1 and
#> have been normalized. Used priorprobs: 0.5 0.2 0.3
#> Hypotheses:
#> 
#>   H1:   "wt<(disp,hp)<0"
#>   H2:   "wt<(disp,hp)=0"
#>   Hc:   "Not H1-H2"
#> 
#> Posterior probability of each hypothesis (rounded):
#> 
#>   H1:   0.995
#>   H2:   0.001
#>   Hc:   0.004
```

As can be seen the posterior probabilites are now different from when we used the default equal prior probabilites.

An alternative to specifying explicit hypotheses is to specify the option "exploratory". This will print the posterior probabilites for the hypotheses "X &lt; 0; X = 0; X &gt; 0" for all variables, including the intercept, and assumes equal prior probabilities.

``` r
test_hyp(fit, "exploratory")
#> Hypotheses:
#> 
#>   H1:   "X < 0"
#>   H2:   "X = 0"
#>   H3:   "X > 0"
#> 
#> Posterior probabilities for each variable (rounded), 
#> assuming equal prior probabilities:
#> 
#>                H1    H2    H3
#>             X < 0 X = 0 X > 0
#> (Intercept) 0.117 0.767 0.117
#> disp        0.125 0.766 0.109
#> hp          0.897 0.098 0.005
#> wt          0.985 0.014 0.001
```

Also here it would be possible to ask for the BF\_matrix for each variable if we had saved the test as an object.

References
----------

Mulder, J. (2014). Prior adjusted default Bayes factors for testing (in) equality constrained hypotheses. Computational Statistics & Data Analysis, 71, 448-463.
