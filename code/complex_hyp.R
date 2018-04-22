
#Complex hypothesis testing with Bayes Factor
#Script purpose: Function to examine complex hypothesis for lm objects with a minimal prior and BF as output
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
#1) Initial setup and checks
#2) convert input into matrices, [gives warnings if hypothesis input improperly]
#3) Check which (if any) of the constraints matrices are NULL and choose computation option based on that
#**3.1) Only equality comparisons
#**3.2) Only inequality comparisons
#**3.3) both
#END

#Note: If several hypotheses specified in input (separated by semicolons) function loops over parts 2-3 for each hypothesis

hyp_test <- function(object, hyp, mcrep = 1e6){
  
  #1) initial setup and checks of input----
  varnames <- variable.names(object) #provides the variable names of the object, including intercept
  if(is.null(varnames)) stop("Please input proper linear model object")
  betahat <- object$coefficients # ML estimates for betas
  
  k <- length(varnames) #varnames length is the same as the number of parameters
  n <- length(object$fitted.values) # df posterior = n - k
  b <- (k + 1) / n #df prior = nb - k
  
  hyp2 <- gsub(" ", "", hyp) #removes all whitespace
  if(!grepl("^[0-9a-zA-Z><=,;().-]+$", hyp2)) stop("Impermissable characters in hypotheses.") #Self-explanatory
  if(grepl("[><=]{2,}", hyp2)) stop("Do not use combined comparison signs e.g., '>=' or '=='")
  
  step1 <- unlist(strsplit(hyp2, split = "[<>=,;()]")) #split by special characters and unlist
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
    
    if(grepl(",", framed$left) || grepl(",", framed$right)){ #Larger loop that deals with commas if the specified hypothesis contains any
      if(nrow(framed) > 1){
        for(r in 1:(nrow(framed)-1)){ #If a hypothesis has been specified with commas e.g., "X1 > 0, X2 > 0" or "(X1, X2) > X3"
          if(all.equal(framed$right[r], framed$left[r+1])){ #The right hand side of the hypothesis df will be equal to the next row left side
            if(substring(framed$right[r], 1, 1) == "(") { #If the first row begins with a ( as when "X1 > (X2, X3)" and opposed to "(X2, X3) > X1"
              framed$right[r] <- sub("),.+", ")", framed$right[r])#If so, remove everything to the right of the parenthesis on the right hand side
              framed$left[r+1] <- sub(".+),", "", framed$left[r +1])#and everything to the left of the parenthesis on the left hand side to correct the df
              } else{
                framed$right[r] <- sub(",.+", "", framed$right[r]) #else, remove everything to the right of the comma on the right hand side
                framed$left[r+1] <- sub("[^,]+,", "", framed$left[r+1]) #and everything to the left of the comma on the left hand side to correct the df
              }
            }
          }
        } 
    
      commas_left <- framed$left[grep(",", framed$left)] #At this point all remaining elements that contain commas should also have parentheses, check this
      commas_right <- framed$right[grep(",", framed$right)] #Necessary to use is isTRUE below in case one of these contains no commas, and 'any' for several rows
      if(isTRUE(any(!grepl("\\(.+)", commas_left))) || isTRUE(any(!grepl("\\(.+)", commas_right))) || #Check so rows contain parenthesis
         isTRUE(any(grepl(").+", commas_left))) || isTRUE(any(grepl(").+", commas_right))) || #Check so parentheses are not followed by anything
         isTRUE(any(grepl(".+\\(", commas_left))) || isTRUE(any(grepl(".+\\(", commas_right)))) { #chekc so parentheses are not preceded by anything
        stop("Incorrect hypothesis syntax or extra character, check specification")
        }
    
      framed$left <- gsub("[()]", "", framed$left) #drop remaining parentheses
      framed$right <- gsub("[()]", "", framed$right)
      commas <- unique(c(grep(",", framed$left), grep(",", framed$right))) #Gives us the unique rows that still contain commas (multiple comparisons) from left or right columns

      if(length(commas) > 0){ #If there are any multiple comparisons e.g., (X1, X2) below loop separates these in
        multiples <- vector("list", length = length(commas)) #Empty vector to store results for each row in loop below

        for(r in seq_along(commas)){ #for each row containing commas
          several <- framed[commas,][r, ] #select row r

          if(several$comp == "="){ #If e.g., (X1, X2) = X3, convert to X1 = X2 = X3
            
            several <- c(several$left, several$right)
            separate <- unlist(strsplit(several, split = ",")) #split by special characters and unlist
            hyp2 <- paste(separate, collapse = "=")
          
            pos_comparisons <- unlist(gregexpr("=", hyp2)) #Gives the positions of all comparison signs
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
            comp <- substring(hyp2, pos_comparisons, pos_comparisons) #Extract comparison signs
            multiples[[r]] <- data.frame(left = leftside, comp = comp, right = rightside, stringsAsFactors = FALSE) #hypotheses as a dataframe
    
            } else{ #If inequality comparison
            leftvars <- unlist(strsplit(several$left, split = ",")) #separate left hand var
            rightvars <- unlist(strsplit(several$right, split = ",")) #separate right hand vars

            left <- rep(leftvars, length.out = length(rightvars)*length(leftvars)) #repeat each leftvars the number of rightvars
            right <- rep(rightvars, each = length(leftvars)) #complement for rightvars
            comp <- rep(several$comp, length(left)) #repeat the comparison a corresponding number of times
    
            multiples[[r]] <- data.frame(left = left, comp = comp, right = right, stringsAsFactors = FALSE) #save as df to be able to combine with 'framed'
            }
          }

        framed <- framed[-commas,] #remove old unfixed rows with commas
        multiples <- do.call(rbind, multiples) #make list into dataframe
        framed <- rbind(multiples, framed) #recombine into one dataframe
      }
    } #end comma loop
    
    equality <- framed[framed$comp == "=",]
    inequality <- framed[!framed$comp == "=",]
    
    #****Equality part string-to-matrix
    if(nrow(equality) == 0) { #If there are no '=' comparisons set to NULL
      R_e <- r_e <- NULL
    } else{
    outcomes <- suppressWarnings(apply(equality[, -2], 2, as.numeric)) #Convert left/right to numeric, non-numeric values (variables) coerced to NA 
    outcomes <- matrix(outcomes, ncol = 2, byrow = FALSE) #Conversion to matrix in case there was only one row in outcomes
    if(any(rowSums(is.na(outcomes)) == 0)) stop("Value compared with value rather than variable, e.g., '2 = 2', check hypotheses")
    rows <- which(rowSums(is.na(outcomes)) < 2) #which rows contain a numeric value (comparing variable to value), that is not two NA-values
    specified <- t(outcomes[rows,]) #transpose so that specified comparison values are extracted in correct order below, e.g, in case when "X1 = 0, 2 = X2"
    specified <- specified[!is.na(specified)] #extract specified comparison values
    r_e <- ifelse(rowSums(is.na(outcomes)) == 2, 0, specified) #If variable = variable -> 0, if variable = value -> value
    r_e <- matrix(r_e) #convert to matrix
    
    var_locations <- apply(equality[, -2], 2, function(x) ifelse(x %in% varnames, match(x, varnames), 0)) #convert non-variables to 0 and others are given their locations in varnames
    var_locations <- matrix(var_locations, ncol = 2) #Necessary if only one comparison row
    
    R_e <- matrix(rep(0, nrow(equality)*length(varnames)), ncol = length(varnames)) #Create empty variable matrix
    
    for(i in seq_along(r_e)){ # for each row i in R_e, replace the columns specified in var_locations row i
      if(!all(var_locations[i, ] > 0)){ #If only one variable is specified (i.e., other one is set to zero)
        R_e[i, var_locations[i,]] <- 1 #Set this variable to 1 in R_e row i
      } else{ #If two variables specified
        R_e[i, var_locations[i,]] <- c(1, -1) #Set one column to 1 and the other to -1 in R_e row i
      }
    }
  }
    
    
    #****Inequality part string-to-matrix
  if(nrow(inequality) == 0) { #If there are no '>' or '<' comparisons set to NULL
    R_i <- r_i <- NULL 
  } else{
    outcomes <- suppressWarnings(apply(inequality[, -2], 2, as.numeric)) #Convert left/right to numeric, non-numeric values (variables) coerced to NA 
    outcomes <- matrix(outcomes, ncol = 2, byrow = FALSE) #Conversion to matrix in case there was only one row in outcomes
    if(any(rowSums(is.na(outcomes)) == 0)) stop("Value compared with value rather than variable, e.g., '2 > 2', check hypotheses")
    cols <- which(rowSums(is.na(outcomes)) < 2) #which columns contain a numeric value (comparing variable to value), that is not two NA-values
    specified <- t(outcomes[cols,]) #transpose so that specified comparison values are extracted in correct order below
    specified <- specified[!is.na(specified)] #extract specified comparison values
    r_i <- ifelse(rowSums(is.na(outcomes)) == 2, 0, specified) #If variable = variable -> 0, if variable = value -> value
    r_i <- matrix(r_i) #convert to matrix
    
    leq <- which(inequality$comp == "<") #gives the rows that contain '<' (lesser or equal) comparisons
    var_locations <- apply(inequality[, -2], 2, function(x) ifelse(x %in% varnames, match(x, varnames), 0)) #convert non-variables to 0 and others are given their locations
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
  }
    
    #3)check comparisons----------------
    if(is.null(R_i)){
      comparisons <- "only equality"
    } else if(is.null(R_e)){
      comparisons <- "only inequality"
    } else{
      comparisons <- "both comparisons"
    }
    
    #set prior mean
    R_ei <- rbind(R_e,R_i) #Sets prior mean around the specified boundary point instead of zero, eg., in case b1 = 2
    r_ei <- rbind(r_e,r_i) 
	  Rr_ei <- cbind(R_ei,r_ei) #Creates adjusted matrix
    beta_zero <- ginv(R_ei)%*%r_ei #NEW: ginv comes from MASS package
    
    if(nrow(Rr_ei) > 1){ #Only relevant if more than one row (and rref only works for >1 row)
      rref_ei <- rref(Rr_ei) #Checks if a common boundary prior exists, e.g., in case b1 > 1, b1 < 2, we would now have two priors, but could be the average between these two
      nonzero <- rref_ei[,k+1]!=0 #Tell whether any of the comparison values are non-zero
      if(max(nonzero)>0){ #If there are any non-zero comparisons
    	  	row1 <- max(which(nonzero==T)) #which row in the adjusted matrix contains the value comparison
    		  if(sum(abs(rref_ei[row1,1:k]))==0){ #if all the variable columns in the adjusted matrix are zero, stop function
    			  stop("Default prior mean cannot be constructed from constraints.")
    			}
      }
    }
    
    if(comparisons == "only equality"){         
      #**3.1)only-equality----

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
      log_BF <- dmvt(x = t(r_e), delta = delta, sigma = scale_post, df = n - k, log = TRUE) - #using logs and backtransforming is more robust
        dmvt(x = t(r_e), delta = delta_zero, sigma = scale_prior, df = n*b - k, log = TRUE) 
      
      BF <- exp(log_BF)#end 'only equality' option
      
    } else if(comparisons == "only inequality"){
      #**3.2)only-inequality----
      
        r_i <- as.vector(r_i) #For pmvt must be a vector contrary to for dmvt.

          if(rankMatrix(R_i)[[1]] == nrow(R_i)){ #If matrix rank is equal to number of rows do exact test.
          	
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
          	if(!is.numeric(mcrep) || !mcrep %% 1 == 0) stop("Input for mcrep should be an integer")
            reps <- mcrep
        	  #Scale matrix for posterior t-distribution.
        	scale_post <- vcov(object)
        
        	#Scale matrix for prior t-distribution
        	scale_prior <- vcov(object) * (n - k) / (n*b - k)
        	          	
          	draws_post <- rmvt(n = reps, delta = betahat, sigma = scale_post, df = n - k) #posterior draws NEW: deltas are changed         
          	draws_prior <- rmvt(n = reps, delta = beta_zero, sigma = scale_prior, df = n*b - k) #prior draws NEW: deltas are changed
          
          	BF <- mean(apply(draws_post%*%t(R_i) > rep(1, reps)%*%t(r_i), 1, prod)) / 
          		mean(apply(draws_prior%*%t(R_i) > rep(1, reps)%*%t(r_i), 1, prod)) #proportion posterior draws satisfying all constrains / prior draws satisfying all constraint
          }
      
    } else{ #If 'both comparisons'
      
      #**3.3)both-comparisons----
      
      #****Equality
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
      log_BF <- dmvt(x = t(r_e), delta = w_post[1:q_e], sigma = matrix(K_post[1:q_e, 1:q_e], ncol = q_e), df = n - k, log = TRUE) - #sigmas must be matrices due to code of dmvt
        dmvt(x = t(r_e), delta = w_prior[1:q_e], sigma = matrix(K_prior[1:q_e, 1:q_e], ncol = q_e), df = n*b - k, log = TRUE)  #using logs and backtransforming is more robust
      
      BFe <- exp(log_BF)
      
      #****Inequality
      R_iv <- R_i %*% ginv(D2) #R_i tilde 
	    r_iv <- r_i - R_i %*% ginv(R_e) %*% r_e #r_i tilde
      
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
      w_2g1_post <- w_2_post + K_21_post %*% solve(K_11_post) %*% matrix(r_e - w_1_post) #w_2 given theta1
      K_2g1_post <- as.vector((n - k + (t(matrix(r_e - w_1_post)) %*% solve(K_11_post) %*% matrix(r_e - w_1_post))) / 
      		(n - k + q_e)) * (K_22_post - K_21_post %*% solve(K_11_post) %*% t(K_21_post)) #K_2 given theta1
      
      #prior
      w_2g1_prior <- w_2_prior + K_21_prior %*% solve(K_11_prior) %*% matrix(r_e - w_1_prior) #w_2 given theta1
      K_2g1_prior <- as.vector((n*b - k + (t(matrix(r_e - w_1_prior)) %*% solve(K_11_prior) %*% matrix(r_e - w_1_prior))) / 
      		(n*b - k + q_e)) * (K_22_prior - K_21_prior %*% solve(K_11_prior) %*% t(K_21_prior)) #K_2 given theta1
      
      if(rankMatrix(R_iv)[[1]] == nrow(R_iv)){ #If matrix rank is equal to number of rows do exact test 
        r_iv <- as.vector(r_iv) #Necessary for pmvt
        
        delta_post <- as.vector(R_iv %*% w_2g1_post) #Transformed posterior mean vector 
        delta_prior <- as.vector(R_iv %*% w_2g1_prior) #Transformed prior mean vector N
        
        scale_post_trans <- R_iv %*% K_2g1_post %*% t(R_iv) #Transformed posterior scale matrix 
        scale_prior_trans <- R_iv %*% K_2g1_prior %*% t(R_iv) #Transformed prior scale matrix 

        if(nrow(scale_post_trans) == 1){ #If univariate 
          BFi <- pt((r_iv - delta_post) / sqrt(scale_post_trans), df = n - k + q_e, lower.tail = FALSE)[1] / #posterior 
            pt((r_iv - delta_prior) / sqrt(scale_prior_trans), df = n*b - k + q_e, lower.tail = FALSE)[1] #prior N
          } else { #if multivariate
            BFi <- pmvt(lower = r_iv, upper = Inf, delta = delta_post, sigma = scale_post_trans, df = n - k + q_e, type = "shifted")[1] / #posterior 
              pmvt(lower = r_iv, upper = Inf, delta = delta_prior, sigma = scale_prior_trans, df = n*b - k + q_e, type = "shifted")[1] #prior 
            }
        
        } else{ #If rank smaller than number of rows, do monte carlo draws
          if(!is.numeric(mcrep) || !mcrep %% 1 == 0) stop("Input for mcrep should be an integer")
          reps <- mcrep

          #Draw from prior and posterior
          draws_post <- rmvt(n = reps, delta = w_2g1_post, sigma = K_2g1_post, df = n - k + q_e) #posterior draws
		      draws_prior <- rmvt(n = reps, delta = w_2g1_prior, sigma = K_2g1_prior, df = n*b - k + q_e) #prior draws

          BFi <- mean(apply(draws_post%*%t(R_iv) > rep(1,reps)%*%t(r_iv),1,prod)) / 
          		mean(apply(draws_prior%*%t(R_iv) > rep(1,reps)%*%t(r_iv),1,prod))
          
        } #End inequality part
      
      #Total BF
      BF <- BFe * BFi
      
      } #end 'both comparisons' option
    
    out[[h]] <- BF #Output for each specified hypothesis
    names(out)[[h]] <- paste0("Bayes Factor for '", hyp[[h]], "' vs. not '", hyp[[h]], "'") #name list with hypothesis as originally specified
    out[[length(out)]] <- out[[length(out)]] * BF #combined BF for all hypotheses
    
  }
  
  names(out)[[length(out)]] <- "Overall Bayes Factor for all hypotheses (if several)" 
  out #final output is a list with all specified hypotheses
  
}
#End----


#***************************************************
#Testing the function----
#***************************************************
d <- sim_reg_data(c(0.2, 4, 0, 0.2, 0.2, 0.8))
q <- lm(y ~ X1 + X2 + X3 + X4 + X5 + X6, data = d)
# object <- q #for testing subsections of the function

hyp <- "X1 = X2, (X4, X5) > X6, X3 < -0.1; (X1, X4) = (0, X3)"
hyp <- "X2 > -10, X1 = X4" #When combined, inequlity comparison part is still bounded at BF = 2

hyp_test(q, hyp)


