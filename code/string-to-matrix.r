#Internship 3 - Joris Mulder, complex hypothesis testing with Bayes Factor
#Script purpose: Creation of function to convert input string of hypothesis to restriction matrices
#Code: Anton Ohlsson Collentine


#****************************************
#This function takes as input a linear model object, and hypothesis as a string e.g.,
hyp <- "X1 > X2 = X3 > X4 = 2"
#Where each variable must contain at least one letter
#and each variable is checked against those in the linear model object
#Permitted input is variables and hypothesized values, and any of the following symbols [><=]
#These symbols have to be single (i.e., not == or =>)
#comparisons are read from left to right in input string
#Output is a list with restriction matrices for equality and inequality (> or <) comparisons

#*************************************
#Testing purposes- function to simulate regression data----
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

#***************************************************
#For testing----
#***************************************************
d <- sim_reg_data(c(0.2, 0.1))
q <- lm(y ~ X1 + X2, data = d)
object <- q

hyp <- "X1 > X2 = 2"

#***************************************************
#Function string-to-matrices----
#***************************************************

create_matrices <- function(object, hyp){
  
  varnames <- variable.names(object) #provides the variable names of the linear model object, including intercept
  if(is.null(varnames)) stop("Please input proper linear model object")
 
  hyp2 <- gsub(" ", "", hyp) #removes all whitespace
  if(!grepl("^[0-9a-zA-Z><=]+$", hyp2)) stop("Impermissable characters in hypotheses. Only letters, numbers and '> < = ' permitted") #Self-explanatory
  if(grepl("[><=]{2,}", hyp2)) stop("Do not use combined comparison signs e.g., '>=' or '=='")
  
  step1 <- unlist(strsplit(hyp2, split = "[<>=]")) #split by comparison signs and unlist
  input_vars <- step1[grep("[a-zA-Z]+", step1)] #extract subunits that contain at least one letter
  if(!all(input_vars %in% varnames)) stop("Hypothesis variable(s) not in object, check spelling") #Checks if input variables exist in lm-object
  
  pos_comparisons <- unlist(gregexpr("[<>=]", hyp2)) #Gives the positions of all comparison signs
  left <- rep(NA, length(pos_comparisons) + 1) #empty vector for loop below
  right <- rep(NA, length(pos_comparisons) + 1) #empty vector for loop below
  pos1 <- c(-1, pos_comparisons) #positions to extract data to the left of comparisons
  pos2 <- c(pos_comparisons, nchar(hyp2) + 1) #positions to extract data to the right of comparisons
  for(i in seq_along(pos1)){
    left[i] <- substring(hyp2, pos1[i] + 1, pos1[i+1] - 1) #Extract all variables or outcomes to the left of a comparison sign
    right[i] <- substring(hyp2, pos2[i] + 1, pos2[i+1] - 1) #Extract all variables or outcomes to the right of a comparison sign
    }
  left <- left[-length(left)] #remove last element which is a NA due to loop formatting
  right <- right[-length(right)] #remove last element which is a NA due to loop formatting
  comparisons <- substring(hyp2, pos_comparisons, pos_comparisons) #Extract comparison signs
  framed <- data.frame(left = left, comp = comparisons, right = right, stringsAsFactors = FALSE) #hypotheses as a dataframe
  
  equality <- framed[framed$comp == "=",]
  inequality <- framed[!framed$comp == "=",]
  
  #********Equality
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
  
  
  #**************Inequality
  greater_than <- inequality[inequality$comp == ">",] #Separate between greater than 
  less_than <- inequality[inequality$comp == "<",] #and less than comparisons
  
  #For greater_than
  if(nrow(greater_than) == 0) { #If there are no '>' comparisons set to NULL
    geq <- NULL #outcome of loop
    } else{
      outcomes <- suppressWarnings(apply(greater_than[, -2], 2, as.numeric)) #Convert left/right to numeric, non-numeric values (variables) coerced to NA 
      outcomes <- matrix(outcomes, ncol = 2, byrow = TRUE) #Conversion to matrix in case there was only one row in outcomes
      if(any(rowSums(is.na(outcomes)) == 0)) stop("Value compared with value rather than variable, e.g., '2 > 2', check hypotheses")
      cols <- which(rowSums(is.na(outcomes)) < 2) #which columns contain a numeric value (comparing variable to value), that is not two NA-values
      specified <- t(outcomes[cols,]) #transpose so that specified comparison values are extracted in correct order below
      specified <- specified[!is.na(specified)] #extract specified comparison values
      r_i <- ifelse(rowSums(is.na(outcomes)) == 2, 0, specified) #If variable = variable -> 0, if variable = value -> value
      r_i <- matrix(r_i, ncol = length(r_i)) #convert to matrix
  
      var_locations <- t(apply(greater_than[, -2], 1, function(x) ifelse(x %in% varnames, which(varnames %in% x), 0))) #convert non-variables to NA, only worked with rows but gets transposed
  
      R_i <- matrix(rep(0, nrow(greater_than)*length(varnames)), ncol = length(varnames)) #Create empty variable matrix
  
      for(i in seq_along(r_i)){ # for each row i in R_i, replace the columns specified in var_locations row i
        if(!all(var_locations[i, ] > 0)){ #If only one variable is specified (i.e., other one is set to zero)
          R_i[i, var_locations[i,]] <- 1 #Set this variable to 1 in R_i row i
          } else{ #If two variables specified
            R_i[i, var_locations[i,]] <- c(1, -1) #Set one column to 1 and the other to -1 in R_i row i
          }
        }
  
      geq <- list(R_i = R_i, r_i = r_i) #list with greater or equal
      }
  
  #For less_than
  if(nrow(less_than) == 0) { #If there are no '<' comparisons set to NULL
    leq <- NULL #outcome of loop
    } else{
      outcomes <- suppressWarnings(apply(less_than[, -2], 2, as.numeric)) #Convert left/right to numeric, non-numeric values (variables) coerced to NA 
      outcomes <- matrix(outcomes, ncol = 2, byrow = TRUE) #Conversion to matrix in case there was only one row in outcomes
      if(any(rowSums(is.na(outcomes)) == 0)) stop("Value compared with value rather than variable, e.g., '2 < 2', check hypotheses")
      cols <- which(rowSums(is.na(outcomes)) < 2) #which columns contain a numeric value (comparing variable to value), that is not two NA-values
      specified <- t(outcomes[cols,]) #transpose so that specified comparison values are extracted in correct order below
      specified <- specified[!is.na(specified)] #extract specified comparison values
      r_i <- ifelse(rowSums(is.na(outcomes)) == 2, 0, specified) #If variable = variable -> 0, if variable = value -> value
      r_i <- matrix(r_i, ncol = length(r_i)) #convert to matrix
  
      var_locations <- t(apply(less_than[, -2], 1, function(x) ifelse(x %in% varnames, which(varnames %in% x), 0))) #convert non-variables to NA, only worked with rows but gets transposed
  
      R_i <- matrix(rep(0, nrow(less_than)*length(varnames)), ncol = length(varnames)) #Create empty variable matrix
  
      for(i in seq_along(r_i)){ # for each row i in R_i, replace the columns specified in var_locations row i
        if(!all(var_locations[i, ] > 0)){ #If only one variable is specified (i.e., other one is set to zero)
          R_i[i, var_locations[i,]] <- 1 #Set this variable to 1 in R_i row i
          } else{ #If two variables specified
            R_i[i, var_locations[i,]] <- c(1, -1) #Set one column to 1 and the other to -1 in R_i row i
          }
        }
  
      leq <- list(R_i = R_i, r_i = r_i) #list with less or equal to
      }
  
  if(is.null(geq) && is.null(leq)){
    list_inequality <- NULL #If no inequality comparisons, set to null
  } else{
    list_inequality <- list(R_i = list(geq = geq$R_i, leq = leq$R_i), r_i = list(geq = geq$r_i, leq = leq$r_i)) #Order list by R_i and r_i
  }
  
  matrices <- list(equality = list_equality, inequality = list_inequality) #final output
  
}

