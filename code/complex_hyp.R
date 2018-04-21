#Internship 3 - Joris Mulder, complex hypothesis testing with Bayes Factor
#Script purpose: Function to examine complex hypothesis for lm objects with a reasonable prior and BF as output
#Code: Anton Ohlsson Collentine

#Content----
#a) Function to simulate regression data
#b) hypothesis testing function
#c) code for trying out b)


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
#Hypothesis testing function----
#*************************************
#Requires mvtnorm
if(!require("mvtnorm")){install.packages("mvtnorm")}
library(mvtnorm)
#reqires Matrix for function 'rankMatrix' for choice between pmvt and montecarlo in only inequality option 
if(!require("Matrix")){install.packages("Matrix")}
library(Matrix)
if(!require("pracma")){install.packages("pracma")} #for function rref (reduced row echelon)
library(pracma)
if(!require("MASS")){install.packages("MASS")} #for function ginv (general inverse)
library(MASS)


#Parts of function
#1) Initial setup and checks, [gives warning if improper lm object input]
#2) convert input into matrices, [gives warnings if hypothesis input improperly]
#3) Check which (if any) of the constraints matrices are NULL and choose computation option based on that
#**3.1) Only equality comparisons
#**3.2) Only inequality comparisons
#**3.3) both
#END

#Note: If several hypotheses specified in input (separated by semicolons) function loops over parts 2-3 for each hypothesis


hyp_test <- function(object, hyp){
  
  #1) initial setup and checks of input----
  varnames <- variable.names(object) #provides the variable names of the object, including intercept
  if(is.null(varnames)) stop("Please input proper linear model object")
  betahat <- object$coefficients # ML estimates for betas
  
  k <- length(varnames) #varnames length is the same as the number of parameters
  n <- length(object$fitted.values) # df posterior = n - k
  b <- (k + 1) / n #df prior = nb - k
  
  hyp2 <- gsub("[ ()]", "", hyp) #removes all whitespace and parentheses
  if(!grepl("^[0-9a-zA-Z><=,;]+$", hyp2)) stop("Impermissable characters in hypotheses. Letters, numbers and > < = (),; permitted") #Self-explanatory
  if(grepl("[><=]{2,}", hyp2)) stop("Do not use combined comparison signs e.g., '>=' or '=='")
  
  step1 <- unlist(strsplit(hyp2, split = "[<>=,;]")) #split by comparison signs and unlist
  input_vars <- step1[grep("[a-zA-Z]+", step1)] #extract subunits that contain at least one letter
  if(!all(input_vars %in% varnames)) stop("Hypothesis variable(s) not in object, check spelling") #Checks if input variables exist in lm-object
  
  hyp <- unlist(strsplit(hyp, split = ";")) #For returning specified hypotheses with outcome
  hyps <- unlist(strsplit(hyp2, split = ";")) #Separated hypotheses (if several hypotheses) for use in computations
  out <- vector("list", length = length(hyps) + 1) #list for final output of each hypothesis, + 1 for product (overall BF)
  out[[length(out)]] <- 1 #First value for computing overall BF at end (in case of several hypotheses)
  
  for(h in seq_along(hyps)){ #loops over the rest of the function until penultimate }
    hyp2 <- hyps[[h]] #for each hypothesis, go through the rest of the function
    
    #2)hyp-to-matrices----
    pos_comparisons <- unlist(gregexpr("[<>=]", hyp2)) #Gives the positions of all comparison signs
    leftside <- rep(NA, length(pos_comparisons) + 1) #empty vector for loop below
    rightside <- rep(NA, length(pos_comparisons) + 1) #empty vector for loop below
    pos1 <- c(-1, pos_comparisons) #positions to extract data to the leftside of comparisons
    pos2 <- c(pos_comparisons, nchar(hyp2) + 1) #positions to extract data to the rightside of comparisons
    for(i in seq_along(pos1)){
      leftside[i] <- substring(hyp2, pos1[i] + 1, pos1[i+1] - 1) #Extract all variables or outcomes to the leftside of a comparison sign
      rightside[i] <- substring(hyp2, pos2[i] + 1, pos2[i+1] - 1) #Extract all variables or outcomes to the rightside of a comparison sign
    }
    leftside <- leftside[-length(leftside)] #remove last element which is a NA due to loop formatting
    rightside <- rightside[-length(rightside)] #remove last element which is a NA due to loop formatting
    comparisons <- substring(hyp2, pos_comparisons, pos_comparisons) #Extract comparison signs
    framed <- data.frame(left = leftside, comp = comparisons, right = rightside, stringsAsFactors = FALSE) #hypotheses as a dataframe
    
    commas <- unique(c(grep(",", framed$left), grep(",", framed$right))) #Gives us the unique rows that contain commas (multiple comparisons) from left or right columns
    if(length(commas) > 0){ #If there are any multiple comparisons e.g., (X1, X2) separate these
    multiples <- vector("list", length = length(commas)) #Empty vector to store results for each row in loop below
    for(r in seq_along(commas)){ #for each row containing commas
      several <- framed[commas,][r, ] #select row r
      leftvars <- unlist(strsplit(several$left, split = ",")) #separate left hand var
      rightvars <- unlist(strsplit(several$right, split = ",")) #separate right hand vars
      if(any(grepl("^$", leftvars)) || any(grepl("^$", rightvars))) stop("Misplaced comma in hypothesis") #if empty element after strsplit
      
      left <- rep(leftvars, each = length(rightvars)) #repeat each leftvars the number of rightvars
      right <- rep(rightvars, each = length(leftvars)) #complement for rightvars
      comp <- rep(several$comp, length(left)) #repeat the comparison a corresponding number of times
      
      multiples[[r]] <- data.frame(left = left, comp = comp, right = right, stringsAsFactors = FALSE) #save as df to be able to combine with 'framed'
    }
    
    framed <- framed[-commas,] #remove old unfixed rows with commas
    multiples <- do.call(rbind, multiples) #make list into dataframe
    framed <- rbind(multiples, framed) #recombine into one dataframe
    }
    
    equality <- framed[framed$comp == "=",]
    inequality <- framed[!framed$comp == "=",]
    
    #****Equality part string-to-matrix
    if(nrow(equality) == 0) { #If there are no '=' comparisons set to NULL
      list_equality <- NULL
    } else{
      outcomes <- suppressWarnings(apply(equality[, -2], 2, as.numeric)) #Convert left/right to numeric, non-numeric values (variables) coerced to NA 
      outcomes <- matrix(outcomes, ncol = 2, byrow = TRUE) #Conversion to matrix in case there was only one row in outcomes
      if(any(rowSums(is.na(outcomes)) == 0)) stop("Value compared with value rather than variable, e.g., '2 = 2', check hypotheses")
      cols <- which(rowSums(is.na(outcomes)) < 2) #which columns contain a numeric value (comparing variable to value), that is not two NA-values
      specified <- t(outcomes[cols,]) #transpose so that specified comparison values are extracted in correct order below
      specified <- specified[!is.na(specified)] #extract specified comparison values
      r_e <- ifelse(rowSums(is.na(outcomes)) == 2, 0, specified) #If variable = variable -> 0, if variable = value -> value
      r_e <- matrix(r_e, ncol = length(r_e)) #convert to matrix
      
      var_locations <- t(apply(equality[, -2], 1, function(x) ifelse(x %in% varnames, which(varnames %in% x), 0))) #convert non-variables to NA, only worked with rows but gets transposed
      var_locations <- matrix(var_locations, ncol = 2) #Necessary if only one comparison row
      
      R_e <- matrix(rep(0, nrow(equality)*length(varnames)), ncol = length(varnames)) #Create empty variable matrix
      
      for(i in seq_along(r_e)){ # for each row i in R_e, replace the columns specified in var_locations row i
        if(!all(var_locations[i, ] > 0)){ #If only one variable is specified (i.e., other one is set to zero)
          R_e[i, var_locations[i,]] <- 1 #Set this variable to 1 in R_e row i
        } else{ #If two variables specified
          R_e[i, var_locations[i,]] <- c(1, -1) #Set one column to 1 and the other to -1 in R_e row i
        }
      }
      list_equality <- list(R_e = R_e, r_e = r_e) #Note column 1 in R_e is for intercept
    }
    
    
    #****Inequality part string-to-matrix
  if(nrow(inequality) == 0) { #If there are no '>' or '<' comparisons set to NULL
    list_inequality <- NULL 
    } else{
      outcomes <- suppressWarnings(apply(inequality[, -2], 2, as.numeric)) #Convert left/right to numeric, non-numeric values (variables) coerced to NA 
      outcomes <- matrix(outcomes, ncol = 2, byrow = TRUE) #Conversion to matrix in case there was only one row in outcomes
      if(any(rowSums(is.na(outcomes)) == 0)) stop("Value compared with value rather than variable, e.g., '2 > 2', check hypotheses")
      cols <- which(rowSums(is.na(outcomes)) < 2) #which columns contain a numeric value (comparing variable to value), that is not two NA-values
      specified <- t(outcomes[cols,]) #transpose so that specified comparison values are extracted in correct order below
      specified <- specified[!is.na(specified)] #extract specified comparison values
      r_i <- ifelse(rowSums(is.na(outcomes)) == 2, 0, specified) #If variable = variable -> 0, if variable = value -> value
      r_i <- matrix(r_i, ncol = length(r_i)) #convert to matrix
  
      leq <- which(inequality$comp == "<") #gives the rows that contain '<' comparisons
      var_locations <- t(apply(inequality[, -2], 1, function(x) ifelse(x %in% varnames, which(varnames %in% x), 0))) #convert non-variables to NA, only worked with rows but gets transposed
      var_locations <- matrix(var_locations, ncol = 2) #Necessary if only one comparison row
      
      R_i <- matrix(rep(0, nrow(inequality)*length(varnames)), ncol = length(varnames)) #Create empty variable matrix
  
      for(i in seq_along(r_i)){ # for each row i in R_i, replace the columns specified in var_locations row i
        if(!all(var_locations[i, ] > 0)){ #If only one variable is specified (i.e., other one is set to zero)
          
          value <- if(i %in% leq) -1 else 1 #If comparison is 'lesser or equal' set to -1, if 'larger or equal' set to 1
          R_i[i, var_locations[i,]] <- value #Set this variable to 1 in R_i row i
          
          } else{ #If two variables specified
            value <- if(i %in% leq) c(-1, 1) else c(1, -1) #If comparison is 'leq' take var2 - var1, if 'larger or equal' take var1 - var2
            R_i[i, var_locations[i,]] <- value #Set one column to 1 and the other to -1 in R_i row i
          }
        }
  
      list_inequality<- list(R_i = R_i, r_i = r_i) #Note column 1 in R_i is for intercept
      }
    
    matrices <- list(equality = list_equality, inequality = list_inequality) #List with matrices of input hypothesis
    
    
    #3)check comparisons----------------
    if(is.null(matrices$inequality)){
      comparisons <- "only equality"
    } else if(is.null(matrices$equality)){
      comparisons <- "only inequality"
    } else{
      comparisons <- "both comparisons"
    }
    
    #set prior mean
    R_ei <- rbind(R_e,R_i) #Sets prior mean around the specified boundary point instead of zero, eg., in case b1 = 2
    r_ei <- rbind(r_e,r_i) #NEW:NB! if one of these is null (no problem, but then the object will not exist at all). Remove matrices listing solves this. 
	  Rr_ei <- cbind(R_ei,r_ei) #Creates adjusted matrix
    beta_zero <- ginv(R_ei)%*%r_ei #NEW: ginv comes from MASS package
    
    rref_ei <- rref(Rr_ei) #NEW: reduced row echelon form. pracma package. 
    nonzero <- rref_ei[,k+1]!=0 #Tell whether any of the comparison values are non-zero
    if(max(nonzero)>0){ #If there are any non-zero comparisons
    		row1 <- max(which(nonzero==T)) #which row in the adjusted matrix contains the value comparison
    		if(sum(abs(rref_ei[row1,1:k]))==0){ #if all the variable columns in the adjusted matrix are zero, stop function
    			stop("Default prior mean cannot be constructed from constraints.")
    			}
    	} #The above paragraph checks if a common boundary prior exists, e.g., in case b1 > 1, b1 < 2, we would now have two priors, but could be the average between these two
    
    if(comparisons == "only equality"){         
      #**3.1)only-equality----
      R_e <- matrices$equality$R_e
      r_e <- matrices$equality$r_e
      
      delta <- R_e %*% betahat #Posterior values we want to check
      delta_zero <- R_e %*% beta_zero
      
      #Scale matrix components, had to separate them to calculate RX
      X <- model.matrix(object) #X-values including intercept
      RX <- R_e %*% solve((t(X) %*% X)) %*% t(R_e) 
      s2 <- sum((model.frame(object)$y - X %*% betahat)^2)
      
      #Scale matrix for posterior t-distribution
      scale_post <- matrix(s2 * RX / (n - k), ncol = nrow(R_e)) #ncol = number of effects, needs to be in matrix for dmvt
      
      #Scale matrix for prior t-distribution
      scale_prior <- matrix(s2 * RX / (n*b - k), nrow(R_e)) #ncol = number of effects, needs to be in matrix for dmvt
      
      #Hypothesis test
      log_BF <- dmvt(x = r_e, delta = delta, sigma = scale_post, df = n - k, log = TRUE) - #using logs and backtransforming is more robust
        dmvt(x = r_e, delta = delta_zero, sigma = scale_prior, df = n*b - k, log = TRUE) 
      
      BF <- exp(log_BF)#end 'only equality' option
      
    } else if(comparisons == "only inequality"){
      #**3.2)only-inequality----
      
        R_i <- matrices$inequality$R_i
        r_i <- as.vector(matrices$inequality$r_i) #For pmvt must be a vector contrary to for dmvt.

          if(rankMatrix(R_i)[[1]] == nrow(R_i)){ #If matrix rank is equal to number of rows do exact test. NEW: line change of place
          	
          delta <- as.vector(R_i %*% betahat) #Posterior values we want to check
        	delta_zero <- as.vector(R_i %*% beta_zero) #Prior values.
        
         	#Scale matrix components
        	X <- model.matrix(object) #X-values including intercept
        	RX <- R_i %*% solve(t(X) %*% X) %*% t(R_i) #Linear transformation
        	s2 <- sum((model.frame(object)$y - X %*% betahat)^2) #sums of squares
        
        	#Scale matrix for posterior t-distribution
        	scale_post <- s2 * RX / (n - k) 
        
        	#Scale matrix for prior t-distribution
        	scale_prior <- s2 * RX / (n*b - k) 
          	
            if(nrow(scale_post) == 1){ #If univariate
              BF <- pt((r_i - delta) / sqrt(scale_post), df = n - k, lower.tail = FALSE)[1] / #posterior
                pt((r_i - delta_zero) / sqrt(scale_prior), df = n*b - k, lower.tail = FALSE)[1] #prior
            } else { #if multivariate
              BF <- pmvt(lower = r_i, upper = Inf, delta = delta, sigma = scale_post, df = n - k, type = "shifted")[1] / #posterior
                pmvt(lower = r_i, upper = Inf, delta = delta_zero, sigma = scale_prior, df = n*b - k, type = "shifted")[1] #prior
            }
          
          } else{#No transformation is possible. Alternative method using monte carlo draws if matrix rank not equal to numer of rows
          	
        	#Scale matrix for posterior t-distribution. NEW: NB! if RX not necessary this below is equal to vcov(object)
        	scale_post <- vcov(object)
        
        	#Scale matrix for prior t-distribution NEW: and this equal to vcov(object) * (n - k) / (n*b - k)
        	scale_prior <- vcov(object) * (n - k) / (n*b - k)
        	          	
          	draws_post <- rmvt(n = 1e6, delta = betahat, sigma = scale_post, df = n - k) #posterior draws NEW: deltas are changed         
          	draws_prior <- rmvt(n = 1e6, delta = beta_zero, sigma = scale_prior, df = n*b - k) #prior draws NEW: deltas are changed
          
          	BF <- mean(apply(draws_post%*%t(R_i) > rep(1, 1e6)%*%t(r_i), 1, prod)) / #NEW: changed computation. We now only have one row, prod superflous?
          		mean(apply(draws_prior%*%t(R_i) > rep(1, 1e6)%*%t(r_i), 1, prod)) #proportion posterior draws satisfying all constrains / prior draws satisfying all constraint
          }
      
    } else{ #If 'both comparisons'
      
      #**3.3)both-comparisons----
      
      #****Equality
      R_e <- matrices$equality$R_e
      r_e <- matrices$equality$r_e
      
      q_e <- nrow(R_e)

      #Scale matrix for posterior t-distribution
      scale_post <- vcov(object) #ncol = number of effects, needs to be in matrix for dmvt
      
      #Scale matrix for prior t-distribution
      scale_prior <- vcov(object) * (n - k) / (n*b - k)#ncol = number of effects, needs to be in matrix for dmvt
      
      #a)Transformation matrix
      D <- diag(k) - t(R_e) %*% solve(R_e %*% t(R_e)) %*% R_e 
      D2 <- unique(D) #Unique, must take unique first or else if only one row treats as vector
      D2 <- D2[as.logical(rowSums(D2 != 0)),] #Remove if only zeroes, this version keeps also rows where sum (+ -) ends up being zero
      Tm <- rbind(R_e, D2) #Transformation matrix, T is an object in base already (TRUE) so using Tm
      
      #b)
      w_post <- Tm %*% betahat #post. mean of xi in paper
      w_prior <- Tm %*% beta_zero 
      K_post <- Tm %*% scale_post %*% t(Tm) #post. scale of xi in paper
      K_prior <- Tm %*% scale_prior %*% t(Tm)
      
      #Equality BF 
      log_BF <- dmvt(x = r_e, delta = w_post[1:q_e], sigma = matrix(K_post[1:q_e, 1:q_e], ncol = q_e), df = n - k, log = TRUE) - #sigmas must be matrices due to code of dmvt
        dmvt(x = r_e, delta = w_prior[1:q_e], sigma = matrix(K_prior[1:q_e, 1:q_e], ncol = q_e), df = n*b - k, log = TRUE)  #using logs and backtransforming is more robust
      
      BFe <- exp(log_BF)
      
      #****Inequality
      R_i <- matrices$inequality$R_i
      r_i <- matrices$inequality$r_i       
      
      R_iv <- R_i %*% ginv(D2) #R_i tilde NEW: added ginv(D2) from package MASS instead of using solve(Tm). 
	    r_iv <- r_i - R_i %*% ginv(R_e) %*% r_e #r_i tilde #NEW: different position, not vector, different expression
      
      #Partitioning to make inequality computations more understandable
      #Posterior part
      w_1_post <- w_post[1:q_e]
      w_2_post <- w_post[(q_e + 1):k]
      
      K_11_post <- K_post[1:q_e, 1:q_e]
      K_12_post <- K_post[1:q_e, (q_e + 1):k]
      K_21_post <- K_post[(q_e + 1):k, 1:q_e]
      K_22_post <- K_post[(q_e + 1):k, (q_e + 1):k]
      
      #prior part
      w_1_prior <- w_prior[1:q_e]
      w_2_prior <- w_prior[(q_e + 1):k]
      
      K_11_prior <- K_prior[1:q_e, 1:q_e]
      K_12_prior <- K_prior[1:q_e, (q_e + 1):k]
      K_21_prior <- K_prior[(q_e + 1):k, 1:q_e]
      K_22_prior <- K_prior[(q_e + 1):k, (q_e + 1):k]
      
      #Conditional mean vectors and scale matrices
      #posterior
      w_2g1_post <- w_2_post + K_21_post %*% solve(K_11_post) %*% matrix(r_e - w_1_post) #w_2 given theta1, last part needs to be transposed to function as a vector in matrix calc.
      K_2g1_post <- as.vector((n - k + (t(matrix(r_e - w_1_post)) %*% solve(K_11_post) %*% matrix(r_e - w_1_post))) / #NEW
      		(n - k + q_e)) * (K_22_post - K_21_post %*% solve(K_11_post) %*% t(K_21_post)) #K_2 given theta1
      
      #prior
      w_2g1_prior <- w_2_prior + K_21_prior %*% solve(K_11_prior) %*% matrix(r_e - w_1_prior) #w_2 given theta1
      K_2g1_prior <- as.vector((n*b - k + (t(matrix(r_e - w_1_prior)) %*% solve(K_11_prior) %*% matrix(r_e - w_1_prior))) / #NEW: parenthesis added again
      		(n*b - k + q_e)) * (K_22_prior - K_21_prior %*% solve(K_11_prior) %*% t(K_21_prior)) #K_2 given theta1
      
      if(rankMatrix(R_iv)[[1]] == nrow(R_iv)){ #If matrix rank is equal to number of rows do exact test NEW: line moved. Now below all partitioning

        delta_post <- as.vector(R_iv %*% w_2g1_post) #Transformed posterior mean vector #NEW: w_2g1_post, NB! object name changed
        delta_prior <- as.vector(R_iv %*% w_2g1_prior) #Transformed prior mean vector NEW: w_2g1_prior, NB! object name changed
        
        scale_post_trans <- R_iv %*% K_2g1_post %*% t(R_iv) #Transformed posterior scale matrix NEW! computation and name, old name scale_m
        scale_prior_trans <- R_iv %*% K_2g1_prior %*% t(R_iv) #Transformed prior scale matrix NEW! computation and name, old name scale_star

        if(nrow(scale_post_trans) == 1){ #If univariate NEW: scale_post instead of scale_m
          BFi <- pt((r_iv - delta_post) / sqrt(scale_post_trans), df = n - k + q_e, lower.tail = FALSE)[1] / #posterior NEW: delta_post, scale_post_trans, added q_e to df
            pt((r_iv - delta_prior) / sqrt(scale_prior_trans), df = n*b - k + q_e, lower.tail = FALSE)[1] #prior NEW: delta_prior, scale_prior_trans, added q_e to df
          } else { #if multivariate
            BFi <- pmvt(lower = r_iv, upper = Inf, delta = delta_post, sigma = scale_post_trans, df = n - k + q_e, type = "shifted")[1] / #posterior NEW: same as for univariate
              pmvt(lower = r_iv, upper = Inf, delta = delta_prior, sigma = scale_prior_trans, df = n*b - k + q_e, type = "shifted")[1] #prior NEW: same as for univariate
            }
        
        } else{ #If rank smaller than number of rows, do monte carlo draws

          #Draw from prior and posterior
          draws_post <- rmvt(n = 1e6, delta = w_2g1_post, sigma = K_2g1_post, df = n - k + q_e) #posterior draws
		      draws_prior <- rmvt(n = 1e6, delta = w_2g1_prior, sigma = K_2g1_prior, df = n*b - k + q_e) #prior draws

          BFi <- mean(apply(draws_post%*%t(R_iv) > rep(1,1e6)%*%t(r_iv),1,prod)) / #NEW: as with only inequality. Same as before, differently expressed
          		mean(apply(draws_prior%*%t(R_iv) > rep(1,1e6)%*%t(r_iv),1,prod))
          
        } #End inequality part
      
      #Total BF
      BF <- BFe * BFi
      
      } #end 'both comparisons' option
    
    out[[h]] <- BF #Output for each specified hypothesis
    names(out)[[h]] <- paste0("Bayes Factor for '", hyp[[h]], "' vs. not '", hyp[[h]], "'") #name list with hypothesis as originally specified
    out[[length(out)]] <- out[[length(out)]] * BF #combined BF for all hypotheses
    
  }
  
  names(out)[[length(out)]] <- "Overall Bayes Factor for all hypotheses (if several)" #NEW: actually not, but does this make sense to have?
  out #final output is a list with all specified hypotheses
  
}
#End----


#***************************************************
#Testing the function----
#***************************************************
d <- sim_reg_data(c(0.2, 0.1, 0, 0.2, 1, 0.8))
q <- lm(y ~ X1 + X2 + X3 + X4 + X5 + X6, data = d)
# object <- q #for testing subsections of the function

hyp <- "X1 = X4 < X5; X3 = 0"

hyp_test(q, hyp)


