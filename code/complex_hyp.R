#Internship 3 - Joris Mulder, complex hypothesis testing with Bayes Factor

#*************************************
#Function to simulate regression data----
#*************************************
sim_reg_data <- function(betas, intercept = 0,  sigma2 = 1, n = 10000){
  beta <- c(intercept, betas)
  Xmat <- matrix(NA, nrow = n, ncol = length(beta)) #prepare the data matrix for predictors
  Xmat[,1] <- 1 #intercept

  #generate predictor data
  for(i in seq_along(beta)[-1]){
    Xmat[,i] <- rnorm(n)
  }
  
  #Generate error
  error <- rnorm(n,sd=sqrt(sigma2))

  #generate outcome variables as a linear function
  #of the predictors
  y <-  rowSums(Xmat %*% diag(beta)) + error #matrix multiplication (row * column), each rowsum in matrix + error is y-variable
  
  dat <- data.frame(Xmat, y) #Put into dataframe because lm function requires dataframe format
  names(dat)[1:length(beta)] <- paste0("X", 0:(length(beta)-1)) #for clarity change names so that intercept is "X0"
  dat
}

##Possible refinements: Allow to state means and SD of rnorm for data
##possible refinement: Allow to specify distribution of data
##Possible refinements: Better defaults?

#*************************************
#Function to check if beta1 == 0----
#*************************************
#Requires mvtnorm
if(!require("mvtnorm")){install.packages("mvtnorm")}
library(mvtnorm)


test_null <- function(lm, data){

#setup
betahat<- lm$coefficients # ML estimates for betas
varnames <- all.vars(lm$terms) #all.vars(q$terms) provides the names of the lm-formula objects, y is the first
y <- varnames[1]

k <- length(varnames) #varnames has DV but not intercept, but the length is the same as the number of parameters
n <- length(lm$fitted.values) # df posterior = n - k
b <- (k + 1) / n #df prior = nb - k

dat <- as.matrix(data[, colnames(data) %in% varnames]) #needed to make sure we only use the variables in the dataset that were used in lm (if there are more)
dat <- cbind(dat, X0 = rep(1, n)) #Add intercept column

#sums of squares
s2 <- sum((dat[, y] - dat[,colnames(dat) != y] %*% betahat)^2) # (y - X %*% betahat)^2
#colnames rather than names function necessary because matrix. varnames[1] == DV.

#Scale matrix for posterior t-distribution
scale_m <- s2*solve(t(dat[,colnames(dat) != y]) %*% dat[, colnames(dat) != y]) / (n - k) 

#Scale matrix for prior t-distribution
scale_star <- s2*solve(t(dat[,colnames(dat) != y]) %*% dat[, colnames(dat) != y]) / (n*b - k)

#Hypothesis test----

#because the code of dmvt requires that ncol(x) == ncol(sigma) I had to convert these to matrices
zero <- matrix(0, ncol = 1)
S1 <- matrix(diag(scale_m)[1], ncol = 1) #diagonal scale matrix value for beta1
S2 <- matrix(diag(scale_star)[1], ncol = 1) #diagonal scale matrix value for beta1

#betahat[2] is beta1
BF2 <- dmvt(x = zero, delta = betahat[2], sigma = S1, df = n - k, log = FALSE) / #for now I've set log = FALSE
  dmvt(x = zero, delta = zero, sigma = S2, df = n*b - k, log = FALSE) #How do I check if  log needs to be true?

names(BF2) <- "BF beta1 = 0"
BF2
}


#*************************************
#Testing----
#*************************************
d <- sim_reg_data(c(1, 0.7), n = 100)
q <- lm(y ~ X1 + X2, data = d)

test_null(q, d)

##Not sure if working properly












