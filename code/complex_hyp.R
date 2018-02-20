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
  
  #This function currently requires that R_e and r_e be input as vectors
  #It can check hypotheses of the types "beta1 = 0", "beta1 = beta2" and "beta1 = beta2 = 0"
  #These correspond to {R_e = c(0, 1, 0), r_e = 0}, {R_e = c(0, 1, -1), r_e = 0} and {R_e = c(0, 1, -1, 0, 1, 0), r_e = c(0, 0)}
  #Also works with more or fewer variables

#setup
betahat <- object$coefficients # ML estimates for betas
varnames <- all.vars(object$terms) #all.vars(q$terms) provides the names of the object-formula objects, y is the first
y <- varnames[1]

k <- length(varnames) #varnames has DV but not intercept, but the length is the same as the number of parameters
n <- length(object$fitted.values) # df posterior = n - k
b <- (k + 1) / n #df prior = nb - k

R_e <- matrix(R_e, nrow = length(r_e), byrow = TRUE) #should be in the shape of a matrix specifying which coefficients we are testing against what
r_e <- matrix(r_e, ncol = length(r_e)) #specifying what values we are testing the coefficients against. NB! ncol(r_e) has to equl nrow(R_e)

delta <- R_e %*% betahat #Posterior values we want to check
delta_zero <- R_e %*% rep(0, k) #Prior values

#Scale matrix components, had to separate them to calculate RX
X <- model.matrix(object) #X-values including intercept
RX <- R_e %*% solve((t(X) %*% X)) %*% t(R_e) 
s2 <- sum((model.frame(object)$y - X %*% betahat)^2) #Tried extracting s2 from vcov(object) but made some mistake, this gives correct result

#Scale matrix for posterior t-distribution
scale_m <- matrix(s2 * RX / (n - k), ncol = length(delta)) #ncol = number of effects, needs to be in matrix for dmvt

#Scale matrix for prior t-distribution
scale_star <- matrix(s2 * RX / (n*b - k), length(delta)) #ncol = number of effects, needs to be in matrix for dmvt

#Hypothesis test
log_BF <- dmvt(x = r_e, delta = delta, sigma = scale_m, df = n - k, log = TRUE) - #using logs and backtransforming is more robust
  dmvt(x = r_e, delta = delta_zero, sigma = scale_star, df = n*b - k, log = TRUE) 

BF <- exp(log_BF)
names(BF) <- "BF of 'beta1 = 0' versus 'beta1 != 0'"

BF
}

#*************************************
#Function to check if beta1 > 0----
#*************************************
#Requires mvtnorm
if(!require("mvtnorm")){install.packages("mvtnorm")}
library(mvtnorm)


test_area <- function(object, R_e = c(0, 1, 0), r_e = 0){
  
  #This function currently requires that R_e and r_e be input as vectors
  #The dmvt function is not working, but the monte carlo draws seem to be working
  #It can check hypotheses of the types "beta1 > 0"
  #This corresponds to {R_e = c(0, 1, 0), r_e = 0}

  #setup
  betahat <- object$coefficients # ML estimates for betas
  varnames <- all.vars(object$terms) #all.vars(q$terms) provides the names of the object-formula objects, y is the first
  y <- varnames[1]
  
  k <- length(varnames) #varnames has DV but not intercept, but the length is the same as the number of parameters
  n <- length(object$fitted.values) # df posterior = n - k
  b <- (k + 1) / n #df prior = nb - k
  
  R_e <- matrix(R_e, nrow = length(r_e), byrow = TRUE) #Input required as vector for now, should be in the shape of a matrix specifying which coefficients we are testing and how, see notes
  r_e <- r_e #For pmvt must be a vector contrary to for dmvt.

  
  delta <- as.vector(R_e %*% betahat) #Posterior values we want to check
  delta_zero <- as.vector(R_e %*% rep(0, k)) #Prior values
  
  #Scale matrix components
  X <- model.matrix(object) #X-values including intercept
  RX <- as.vector(R_e %*% solve((t(X) %*% X)) %*% t(R_e)) #Needs to be vector for later calculation
  s2 <- sum((model.frame(object)$y - X %*% betahat)^2) #Using this gives correct result

  #Scale matrix for posterior t-distribution
  scale_m <- matrix(s2 * RX / (n - k), ncol = length(delta)) #must be matrix for monte carlo draws
  
  #Scale matrix for prior t-distribution
  scale_star <- matrix(s2 * RX / (n*b - k), ncol = length(delta))

  #Hypothesis test. [Gives error when testing only one coefficient, can't figure out why. Strange result when testing several] Pmvt requires delta and r_e to be vectors. Default type for pmvt is not "shifted" as we want so must be specified.
  BF_more_than <- pmvt(lower = r_e, upper = Inf, delta = delta, sigma = scale_m, df = n - k, type = "shifted") / #These probabilities are not complementary to those below (testing several betas)
    pmvt(lower = r_e, upper = Inf, delta = delta_zero, sigma = scale_star, df = n*b - k, type = "shifted") #so obviously something is incorrect
  
  BF_less_than <- pmvt(lower = -Inf, upper = r_e, delta = delta, sigma = scale_m, df = n - k, type = "shifted") / #
    pmvt(lower = -Inf, upper = r_e, delta = delta_zero, sigma = scale_star, df = n*b - k, type = "shifted") 
  
  BF <- as.vector(BF_more_than / BF_less_than) #remove attributes cluttering output by turning into vector
  
  #Alternative method using monte carlo draws,[seems to work] length(delta) == nrow(sigma) -> sigma must be matrix
  set.seed(56)
  draws_post <- rmvt(n = 100, delta = delta, sigma = scale_m, df = n - k)
  draws_pre <- rmvt(n = 100, delta = delta_zero, sigma = scale_star, df = n*b - k)
  
  BF2 <- mean(draws_post > r_e) / mean(draws_pre > r_e)
  
  names(BF) <- names(BF2) <-  "BF of 'beta1 > 0' versus 'beta1 < 0'"
  
  c(BF, BF2)
}

#currently function is testing if beta1 > 0 vs. beta1 < 0

#*************************************
#Testing----
#*************************************
d <- sim_reg_data(c(0, 1))
q <- lm(y ~ X1 + X2, data = d)

test_equality(q, c(0, 1, 0))
test_area(q)

R_e <- c(0, 1, 0, 0, 0, 1)
r_e <- c(0, 0 )
