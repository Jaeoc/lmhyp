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
#Function to test equality constraints----
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
varnames <- variable.names(object) #provides the variable names of the object, including intercept

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
scale_m <- matrix(s2 * RX / (n - k), ncol = nrow(R_e)) #ncol = number of effects, needs to be in matrix for dmvt

#Scale matrix for prior t-distribution
scale_star <- matrix(s2 * RX / (n*b - k), nrow(R_e)) #ncol = number of effects, needs to be in matrix for dmvt

#Hypothesis test
log_BF <- dmvt(x = r_e, delta = delta, sigma = scale_m, df = n - k, log = TRUE) - #using logs and backtransforming is more robust
  dmvt(x = r_e, delta = delta_zero, sigma = scale_star, df = n*b - k, log = TRUE) 

BF <- exp(log_BF)

# names(BF) <- "BF of 'beta1 = 0' versus 'beta1 != 0'" #needs to be generalized depending on input R_e
names(BF) <- "BF" #temporary

BF
}

#*************************************
#Function to test inequality constraints----
#*************************************
#Requires mvtnorm
if(!require("mvtnorm")){install.packages("mvtnorm")}
library(mvtnorm)


test_inequality <- function(object, R_i = c(0, 1, 0), r_i = 0){
  
  #This function currently requires that R_i and r_i be input as vectors
  #It can check hypotheses of the types "beta1 > 0" and "beta1 > 0, beta2 > 0"
  #This corresponds to {R_i = c(0, 1, 0), r_i = 0} and {R_i = c(0, 1, 0, 0, 0, 1), r_i = c(0, 0)}
  #Also works with more variables

  #setup
  betahat <- object$coefficients # ML estimates for betas
  varnames <- variable.names(object) #provides the variable names of the object, including intercept
  
  k <- length(varnames) #varnames has DV but not intercept, but the length is the same as the number of parameters
  n <- length(object$fitted.values) # df posterior = n - k
  b <- (k + 1) / n #df prior = nb - k
  
  R_i <- matrix(R_i, nrow = length(r_i), byrow = TRUE) #Input required as vector for now, should be in the shape of a matrix specifying which coefficients we are testing and how, see notes
  r_i <- r_i #For pmvt must be a vector contrary to for dmvt.

  
  delta <- as.vector(R_i %*% betahat) #Posterior values we want to check
  delta_zero <- as.vector(R_i %*% rep(0, k)) #Prior values
  
  #Scale matrix components
  X <- model.matrix(object) #X-values including intercept
  RX <- as.vector(R_i %*% solve((t(X) %*% X)) %*% t(R_i)) #Needs to be vector for later calculation
  s2 <- sum((model.frame(object)$y - X %*% betahat)^2) #Using this gives correct result

  #Scale matrix for posterior t-distribution
  scale_m <- matrix(s2 * RX / (n - k), ncol = nrow(R_i)) #must be matrix for monte carlo draws
  
  #Scale matrix for prior t-distribution
  scale_star <- matrix(s2 * RX / (n*b - k), ncol = nrow(R_i))

  #Hypothesis test using exact values
  if(nrow(scale_m) == 1){ #If univariate
    BF <- pt((r_i - delta) / sqrt(scale_m), df = n - k, lower.tail = FALSE)[1] / #posterior
      pt((r_i - delta_zero) / sqrt(scale_star), df = n*b - k, lower.tail = FALSE)[1] #prior
  } else { #if multivariate
  BF <- pmvt(lower = r_i, upper = Inf, delta = delta, sigma = scale_m, df = n - k, type = "shifted")[1] / #posterior
    pmvt(lower = r_i, upper = Inf, delta = delta_zero, sigma = scale_star, df = n*b - k, type = "shifted")[1] #prior
  }
  
  #Alternative method using monte carlo draws
  draws_post <- rmvt(n = 1e6, delta = delta, sigma = scale_m, df = n - k) #posterior draws
  satisfied_post2 <- apply(draws_post > r_i, 1, prod) #checks which posterior draws satisfy constraints

  draws_pre <- rmvt(n = 1e6, delta = delta_zero, sigma = scale_star, df = n*b - k) #prior draws
  satisfied_pre <- apply(draws_pre > r_i, 1, prod) #checks which prior draws satisfy constraints
  
  BF2 <- mean(satisfied_post) / mean(satisfied_pre) #proportion posterior draws satisfying all constrains / prior draws satisfying all constraints
  
  names(BF) <-  "BF pmvt" #needs to be generalized depending on input R_i
  names(BF2) <-  "BF Monte carlo" #needs to be generalized depending on input R_i
  
  list("beta1 > 0, beta2 > 0", BF, BF2)
}

#*************************************
#Function to test both equality and inequality----
#*************************************
test_hyp <- function(object, R_e = NULL, r_e = NULL, R_i = NULL, r_i = NULL){

  ##Function assuming both R_e and R_i are filled in
  
  #setup [Is common in all cases]
  betahat <- object$coefficients # ML estimates for betas
  varnames <- variable.names(object) #provides the variable names of the object, including intercept
  
  k <- length(varnames) #varnames has DV but not intercept, but the length is the same as the number of parameters
  n <- length(object$fitted.values) # df posterior = n - k
  b <- (k + 1) / n #df prior = nb - k
  
  #Scale matrix for posterior t-distribution
  scale_m <- vcov(object) #ncol = number of effects, needs to be in matrix for dmvt
  
  #Scale matrix for prior t-distribution
  scale_star <- vcov(object) * (n - k) / (n*b - k)#ncol = number of effects, needs to be in matrix for dmvt
  
  #Input formatting
  R_e <- matrix(R_e, nrow = length(r_e), byrow = TRUE) #should be in the shape of a matrix specifying which coefficients we are testing against what
  r_e <- matrix(r_e, ncol = length(r_e)) #specifying what values we are testing the coefficients against. NB! ncol(r_e) has to equl nrow(R_e)
  
  R_i <- matrix(R_i, nrow = length(r_i), byrow = TRUE) #should be in the shape of a matrix specifying which coefficients we are testing against what
  r_i <- matrix(r_i, ncol = length(r_i)) #specifying what values we are testing the coefficients against. NB! ncol(r_e) has to equl nrow(R_e)
  
  ###New stuff!
  
  #a)Transformation matrix
  D <- diag(k) - t(R_e) %*% solve(R_e %*% t(R_e)) %*% R_e 
  D2 <- unique(D) #Unique, must take unique first or else if only one row treats as vector
  D2 <- D2[as.logical(rowSums(D2 != 0)),] #Remove if only zeroes, this version keeps also rows where sum (+ -) ends up being zero
  Tm <- rbind(R_e, D2) #Transformation matrix, T is an object in base already (TRUE) so using Tm

  
  #b)
  w_post <- Tm %*% betahat
  w_prior <- Tm %*% rep(0, k) 
  K_post <- Tm %*% scale_m %*% t(Tm) 
  K_prior <- Tm %*% scale_star %*% t(Tm)
  
  #Equality BF - 
  log_BF <- dmvt(x = r_e, delta = w_post[1:nrow(R_e)], sigma = matrix(K_post[1:nrow(R_e), 1:nrow(R_e)], ncol = nrow(R_e)), df = n - k, log = TRUE) - #sigmas must be matrices due to code of dmvt
    dmvt(x = r_e, delta = w_prior[1:nrow(R_e)], sigma = matrix(K_prior[1:nrow(R_e), 1:nrow(R_e)], ncol = nrow(R_e)), df = n*b - k, log = TRUE)  #using logs and backtransforming is more robust
  
  BF1 <- exp(log_BF)

  #Inequality--------------------------------------------------------
  R_i2 <- R_i %*% solve(Tm) #R_i with tilde
  
  #Partitioning----
  #Makes inequality computations more understandable
  q_e <- nrow(R_e) #Used a lot in below calculations
  
  #Posterior parameters
  w_1_post <- w_post[1:q_e]
  w_2_post <- w_post[(q_e + 1):k]
  
  K_11_post <- K_post[1:q_e, 1:q_e]
  K_12_post <- K_post[1:q_e, (q_e + 1):k]
  K_21_post <- K_post[(q_e + 1):k, 1:q_e]
  K_22_post <- K_post[(q_e + 1):k, (q_e + 1):k]
  
  #prior parameters
  w_1_prior <- w_prior[1:q_e]
  w_2_prior <- w_prior[(q_e + 1):k]
  
  K_11_prior <- K_prior[1:q_e, 1:q_e]
  K_12_prior <- K_prior[1:q_e, (q_e + 1):k]
  K_21_prior <- K_prior[(q_e + 1):k, 1:q_e]
  K_22_prior <- K_prior[(q_e + 1):k, (q_e + 1):k]

  
  #Conditional parameters----

  #posterior
  w_2g1_post <- w_2_post + K_21_post %*% solve(K_11_post) %*% matrix(r_e - w_1_post) #w_2 given theta1, last part needs to be transposed to function as a vector in matrix calc.
  
  K_2g1_post <- as.vector(n - k + (t(matrix(r_e - w_1_post)) %*% solve(K_11_post) %*% matrix(r_e - w_1_post)) / #Scalar needs to be vector for multiplication
    (n - k + q_e)) * (K_22_post - K_21_post %*% solve(K_11_post) %*% t(K_21_post)) #K_2 given theta1

  #prior
  w_2g1_prior <- w_2_prior + K_21_prior %*% solve(K_11_prior) %*% matrix(r_e - w_1_prior) #w_2 given theta1
  
  K_2g1_prior <- as.vector(n*b - k + (t(matrix(r_e - w_1_prior)) %*% solve(K_11_prior) %*% matrix(r_e - w_1_prior)) / #Scalar needs to be vector for multiplication
                            (n*b - k + q_e)) * (K_22_prior - K_21_prior %*% solve(K_11_prior) %*% t(K_21_prior)) #K_2 given theta1
                  
  #Draws----
  r_e_matrix <- matrix(rep(r_e,1e6),nrow=1e6) #matrix used for extending draws
  r_i_matrix <- matrix(rep(r_i,1e6),nrow=1e6) #matrix used for checking inequality constraints
  
  #posterior
  draws_post <- rmvt(n = 1e6, delta = w_2g1_post, sigma = K_2g1_post, df = n - k + q_e) #posterior draws
  draws_post2 <- cbind(r_e_matrix, draws_post) #combine r_e + draws
  satisfied_post <- apply(draws_post2 %*% t(R_i2) > r_i_matrix,1,prod) #Check which posterior draws satisfy the inequality constraints
  
  #prior
  draws_prior <- rmvt(n = 1e6, delta = w_2g1_prior, sigma = K_2g1_prior, df = n*b - k + q_e) #prior draws
  draws_prior2 <- cbind(r_e_matrix, draws_prior) #combine r_e + draws
  satisfied_prior <- apply(draws_prior2 %*% t(R_i2) > r_i_matrix, 1, prod) #Check which prior draws satisfy the inequality constraints
  
  #Results----

  #Inequality BF
  BF2 <- mean(satisfied_post) / mean(satisfied_prior)
  
  #Total BF
  BF <- BF1 * BF2
  
  list(BF = BF, BF1 = BF1, BF2 = BF2)

  
}


#*************************************
#Testing functions----
#*************************************
d <- sim_reg_data(c(0.2, 0.1))
q <- lm(y ~ X1 + X2, data = d)

test_equality(q, R_e = c(0, 1, -1, 0, 0, 1), r_e = c(0, 0)) #beta1 == beta2 == 0?
test_inequality(q, R_i = c(0, 0, 1, 0, 1, 0), r_i = c(0, 0)) #beta1 > 0, beta2 > 0?

#for troubleshooting if necessary
object <- q
R_e <- c(0, 1, 0) 
r_e <- c(0)
R_i <- c(0, 1, 0) 
r_i <- c(0)

