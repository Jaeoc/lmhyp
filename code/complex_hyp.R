#Internship 3 - Joris Mulder, complex hypothesis testing with Bayes Factor
#Code: Anton Ohlsson Collentine

#*************************************
#Function to simulate regression data----
#*************************************
sim_reg_data <- function(betas, intercept = 0,  sigma2 = 1, n = 100){
  beta <- c(intercept, betas)
  Xmat <- matrix(NA, nrow = n, ncol = length(beta)) #prepare the data matrix for predictors
  Xmat[,1] <- 1 #intercept

  #generate predictor data
  Xmat[,-1] <- rnorm((length(beta)-1)*n)
  
  #Generate error
  error <- rnorm(n,sd=sqrt(sigma2))

  #generate outcome variables as a linear function
  #of the predictors
  y <-  Xmat %*% beta + error #matrix multiplicaton

  dat <- data.frame(Xmat, y) #Put into dataframe because lm function requires dataframe format
  names(dat)[1:length(beta)] <- paste0("X", 0:(length(beta)-1)) #for clarity change names so that intercept is "X0"
  dat
}

##Possible refinements: Allow to state means and SD of rnorm for data
### If so, generate data with the follwing code: Xmat[-1,] <- rmvnorm(n,mean=c(1,2,3),sigma=diag(rep(1,3)))
##possible refinement: Allow to specify distribution of data (for loop then becomes necessary)
##Possible refinements: Other defaults?

#*************************************
#Function to check if beta1 == 0----
#*************************************
#Requires mvtnorm
if(!require("mvtnorm")){install.packages("mvtnorm")}
library(mvtnorm)


test_equality <- function(object, R_e = c(0, 1, 0), r_e = 0){

#setup
betahat <- object$coefficients # ML estimates for betas
varnames <- all.vars(object$terms) #all.vars(q$terms) provides the names of the object-formula objects, y is the first
y <- varnames[1]

k <- length(varnames) #varnames has DV but not intercept, but the length is the same as the number of parameters
n <- length(object$fitted.values) # df posterior = n - k
b <- (k + 1) / n #df prior = nb - k

R_e <- t(matrix(R_e, ncol = length(r_e))) #Input required as vector for now, should be in the shape of a matrix specifying which coefficients we are testing and how, see notes
r_e <- as.matrix(r_e) #should be in the shape of a 1-col matrix specifying what values we are testing the coefficients against

delta <- R_e %*% betahat #This is the value we're actually checking now
delta_zero <- R_e %*% as.matrix(rep(0, k)) #What we're comparing against

X <- model.matrix(object) #X-values including intercept

RX <- as.vector(R_e %*% solve((t(X) %*% X)) %*% t(R_e)) #Needs to be vector for later calculations

s2 <- sum((model.frame(object)$y - X %*% betahat)^2) #Using this gives correct result

#Scale matrix for posterior t-distribution
scale_m <- s2 * RX / (n - k)

#Scale matrix for prior t-distribution
scale_star <- s2 * RX / (n*b - k)

#because the code of dmvt requires that ncol(x) == ncol(sigma) I had to convert these to matrices
S1 <- as.matrix(scale_m) #posterior scale matrix value for beta1
S2 <- as.matrix(scale_star) #prior scale matrix value for beta1

#Hypothesis test, betahat[2] is beta1
log_BF <- dmvt(x = r_e, delta = delta, sigma = S1, df = n - k, log = TRUE) - #using logs and backtransforming is more robust
  dmvt(x = r_e, delta = delta_zero, sigma = S2, df = n*b - k, log = TRUE) 

BF <- exp(log_BF)
names(BF) <- "BF of 'beta1 = 0' versus 'beta1 != 0'"

BF
}


##I need to find S^2 to be able to extract (X'X)^-1 from the vcov(object). DONE, better way?
#methods(class = "lm") tells us which functions work. 
# model.matrix(q) gives us the X-values.

##Fixed my S2 and the RX calculations, now effect is consistent

##Problem, ncol(x) needs to == ncol(sigma). previously I just used scale_m[2,2]
##but what should I do now? Say I wish the check if beta1 == beta2, what do I use for sigma?

#*************************************
#Function to check if beta1 > 0----
#*************************************
#Requires mvtnorm
if(!require("mvtnorm")){install.packages("mvtnorm")}
library(mvtnorm)


test_area <- function(object, R_e = c(0, 1, 0), r_e = 0){
  
  #setup
  betahat <- object$coefficients # ML estimates for betas
  varnames <- all.vars(object$terms) #all.vars(q$terms) provides the names of the object-formula objects, y is the first
  y <- varnames[1]
  
  k <- length(varnames) #varnames has DV but not intercept, but the length is the same as the number of parameters
  n <- length(object$fitted.values) # df posterior = n - k
  b <- (k + 1) / n #df prior = nb - k
  
  R_e <- t(matrix(R_e, ncol = length(r_e))) #Input required as vector for now, should be in the shape of a matrix specifying which coefficients we are testing and how, see notes
  r_e <- r_e #For pmvt must be a vector contrary to for dmvt. Specifies our null hypotheses
  
  delta <- as.vector(R_e %*% betahat) #This is the value we're actually checking now
  delta_zero <- as.vector(R_e %*% as.matrix(rep(0, k))) #What we're comparing against
  
  X <- model.matrix(object) #X-values including intercept
  
  RX <- as.vector(R_e %*% solve((t(X) %*% X)) %*% t(R_e)) #Needs to be vector for later calculations

  s2 <- sum((model.frame(object)$y - X %*% betahat)^2) #Using this gives correct result

  #Scale matrix for posterior t-distribution
  scale_m <- s2 * RX / (n - k)
  
  #Scale matrix for prior t-distribution
  scale_star <- s2 * RX / (n*b - k)

  #Hypothesis test, sigma[2,2] refers to beta1. Pmvt requires delta, r_e and sigma to be vectors, pmvt also has different default type
  BF_more_than <- pmvt(lower = r_e, upper = Inf, delta = delta, sigma = scale_m, df = n - k, type = "shifted") / #
    pmvt(lower = r_e, upper = Inf, delta = delta_zero, sigma = scale_star, df = n*b - k, type = "shifted") 
  
  BF_less_than <- pmvt(lower = -Inf, upper = r_e, delta = delta, sigma = scale_m, df = n - k, type = "shifted") / #
    pmvt(lower = -Inf, upper = r_e, delta = delta_zero, sigma = scale_star, df = n*b - k, type = "shifted") 
  
  BF <- as.vector(BF_more_than / BF_less_than) #remove attributes cluttering output by turning into vector
  
  #Alternative method using monte carlo draws, length(delta) == nrow(sigma) -> sigma must be matrix
  #because the code of rmvt requires that ncol(x) == ncol(sigma) I had to convert these to matrices
  S1 <- as.matrix(scale_m) #posterior scale matrix value for beta1
  S2 <- as.matrix(scale_star) #prior scale matrix value for beta1
  
  
  set.seed(1)
  draws_post <- rmvt(n = 1e4, delta = delta, sigma = S1, df = n - k)
  draws_pre <- rmvt(n = 1e4, delta = delta_zero, sigma = S2, df = n*b - k)
  
  BF2 <- mean(draws_post > 0) / mean(draws_pre > 0)
  
  names(BF) <- "BF of 'beta1 > 0' versus 'beta1 < 0'"
  
  c(BF, BF2)
}

#currently function is testing if beta1 > 0 vs. beta1 < 0

pmvt(lower = r_e, upper = Inf, delta = delta, corr = 1, df = n - k, type = "shifted") 

#*************************************
#Testing----
#*************************************
d <- sim_reg_data(c(2, 0.7))
q <- lm(y ~ X1 + X2, data = d)

test_equality(q)
test_area(q)

R_e <- c(0, 0, 1)
r_e <- c(0)
