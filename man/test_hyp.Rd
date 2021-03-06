% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_hyp.R
\name{test_hyp}
\alias{test_hyp}
\title{Testing Informed Hypotheses}
\usage{
test_hyp(object, hyp, priorprob = 1, mcrep = 1e+06)
}
\arguments{
\item{object}{A regression model object fit using the \code{lm} function.}

\item{hyp}{A string specifying hypotheses to be tested using the variable
names of the \code{lm} object, or the string \dQuote{exploratory}, see
details.}

\item{priorprob}{Vector of prior probabilities for the input hypotheses, by
default equal.}

\item{mcrep}{Integer specifying the number of iterations if no analytical
solutions is possible. This is rare and only the case if the rank of the
constraint matrix is less than its number of rows.}
}
\description{
Test competing hypotheses about coefficients in an \code{\link{lm}}-object.
}
\details{
This function is based on a method by Mulder (2014), a modification of the
fractional Bayes factor approach. In essence, it uses a number of observations
equal to the number of predictors in the \code{lm} model to construct a
minimally informative prior, and the remainder of the observations are then
used to test the hypotheses.

Hypotheses are specified using the variable names from the \code{lm} object.
If a hypothesis involves multiple variables it is usually preferable to
standardize relevant variables before fitting the model with \code{lm} to
facilitate interpretation. This is done simply by substracting the mean of a
variable from each observation and dividing by the standard deviation. A
simple option for achieving this is to use the \code{\link{scale}} function.

Multiple hypotheses can be specified at the same time by separating them with
a semicolon. It is advisable to only specify competing hypotheses in this way,
that is, hypotheses regarding the same variables, e.g., \dQuote{X1 > 0; X1 <
0; X1 = 0}. If specifying multiple hypotheses and comparing against a value it
is currently only possible to compare against the same value, i.e., \dQuote{X1
= 0; X1 = 2} is not functional input. This is because the prior is centered
around the input value (or zero if no input value), which is not possible in
the case of several comparison values.

Parentheses can be used to compare multiple variables with the same variable
or value. For example, \dQuote{(X1, X2) > X3} is read as \dQuote{X1 > X3 and
X2 > X3}. Each variable should only be specified once in a single hypothesis.

An alternative to specifying hypotheses is to input the string
\dQuote{exploratory}. This will compare the hypotheses \dQuote{X < 0; X = 0; X
> 0} for all independent variables in the regression object, including the
intercept.

For each specified hypothesis the posterior probability is output. If the
hypotheses are not exhaustive (i.e., do not cover the entire parameter space)
this includes the posterior probability of the complement to the input
hypotheses. The complement is the hypothesis that neither of the input
hypotheses is true. For example, inputting \dQuote{X1 > 0; X1 < 0} gives
posterior probabilities for only for these hypotheses, whereas inputting
\dQuote{(X1, X2) > 0} gives posterior probabilities for \dQuote{(X1, X2) > 0}
and \dQuote{not (X1, X2) > 0}.

If not using the \dQuote{exploratory} option, it is possible to specify prior
probabilities for the input hypotheses. By default these are equal (priorprob
= 1). Prior probabilites can both be input as probabilites, e.g., c(0.2, 0.3,
0.5) or relative weights of each hypothesis, e.g, c(2, 3, 5). If the input
probabilites do not sum to 1 they will simply be normalized. Prior
probabilities must be specified for all hypotheses, including the complement
if one exists.

By saving the test as an object it is also possible to access the
\code{BF_matrix} which compares the hypotheses directly against each other
(see examples). This matrix divides the row hypothesis by each column
hypothesis and, assuming equal prior probabilities, can be interpreted as
\dQuote{given the data, [row hypothesis] is [value] times as likely as [column
hypothesis]}.

Supplementary output intended to provide a deeper understanding of the
underlying method and primary output can be printed when the test has been
saved as an object. Calling \code{BF_computation} illustrates the computation
of the Bayes factor of each hypothesis against the unconstrained hypothesis.
In the output \code{c(E)} is the prior density, \code{c(I|E)} the prior
probability, \code{c} the product of these two, and columns prefixed with "f"
the equivalent for the posterior. \code{B(t,u)} is the Bayes factor of
hypothesis t against the unconstrained (u) and \code{PP(t)} is the
posterior probability of t. Finally, calling \code{BFu_CI} provides 90\%
credibility intervals for those few cases where the Bayes factor was
calculated numerically.
}
\section{References}{
 Mulder, J. (2014). Prior adjusted default Bayes factors
 for testing (in) equality constrained hypotheses. Computational Statistics &
 Data Analysis, 71, 448-463.
}

\examples{
###Standardize variables and fit the linear model
dt <- as.data.frame(scale(mtcars[, c(1, 3:4, 6)]))
fit <- lm(mpg ~ disp + hp + wt, data = dt)

###Exploratory analysis
test_hyp(fit, "exploratory")

###Define hypotheses based on theory and test them
hyp <- "(wt, hp) > disp > 0; (wt, hp) > disp = 0"
res <- test_hyp(fit, hyp)
res

###Bayes factor comparison of hypotheses
res$BF_matrix

}
