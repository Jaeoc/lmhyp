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
##possible refinement: Allow to specify distribution of data
##Possible refinements: Other defaults?

#*************************************
#Function to check if beta1 == 0----
#*************************************
#Requires mvtnorm
if(!require("mvtnorm")){install.packages("mvtnorm")}
library(mvtnorm)


test_null <- function(object){

#setup
betahat<- object$coefficients # ML estimates for betas
varnames <- all.vars(object$terms) #all.vars(q$terms) provides the names of the object-formula objects, y is the first
y <- varnames[1]

k <- length(varnames) #varnames has DV but not intercept, but the length is the same as the number of parameters
n <- length(object$fitted.values) # df posterior = n - k
b <- (k + 1) / n #df prior = nb - k

#Scale matrix for posterior t-distribution
scale_m2 <- vcov(object)

#Scale matrix for prior t-distribution
scale_star <- vcov(object) * (n - k) / (n*b - k)

#because the code of dmvt requires that ncol(x) == ncol(sigma) I had to convert these to matrices
S1 <- as.matrix(scale_m[2,2]) #posterior scale matrix value for beta1
S2 <- as.matrix(scale_star[2,2]) #prior scale matrix value for beta1

#Hypothesis test, betahat[2] is beta1
log_BF <- dmvt(x = 0, delta = betahat[2], sigma = S1, df = n - k, log = TRUE) - #using logs and backtransforming is more robust
  dmvt(x = 0, delta = 0, sigma = S2, df = n*b - k, log = TRUE) 

BF <- exp(log_BF)
names(BF) <- "BF of 'beta1 = 0' versus 'beta1 != 0'"

BF
}


#*************************************
#Testing----
#*************************************
d <- sim_reg_data(c(0, 0.7))
q <- lm(y ~ X1 + X2, data = d)

test_null(q)


