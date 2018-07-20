
#Project: lmhyp - Informed hypothesis testing for regression
#Script purpose: Function to examine complex hypothesis for lm objects with a minimal prior and BF as output
#Code: Anton Ohlsson Collentine

#*************************************
#Hypothesis testing function----
#*************************************

#'Testing Informed Hypotheses
#'
#'Test hypotheses about continous predictors in an \code{\link{lm}}-object.
#'
#'This function is based on a method by Mulder (2014), a modification of the
#'fractional Bayes factor approach. In essence, it uses a number of observations
#'equal to the number of predictors in the \code{lm} model to construct a
#'minimally informative prior, and the remainder of the observations are then
#'used to test the hypotheses.
#'
#'It is advisable to standardize relevant variables before fitting the model
#'with \code{lm}. This is done simply by substracting the mean of a variable
#'from each observation and dividing by the standard deviation. A simple option
#'for achieving this is to use the \code{\link{scale}} function.
#'
#'Multiple hypotheses can be specified at the same time by separating them with
#'a semicolon. It is advisable to only specify competing hypotheses in this way,
#'that is, hypotheses regarding the same variables, e.g., \dQuote{X1 > 0; X1 <
#'0; X1 = 0}. If specifying multiple hypotheses and comparing against a value it
#'is currently only possible to compare against the same value, e.g., \dQuote{X1
#'= 0; X1 = 2} is not functional input. This is because the prior is centered
#'around the input value (or zero if no input value), which is not possible in
#'the case of several input values.
#'
#'Parentheses can be used to compare multiple variables with the same variable
#'or value. For example, \dQuote{(X1, X2) > X3} is read as \dQuote{X1 > X3 and
#'X2 > X3}. Each variable should only be specified once in a single hypothesis.
#'
#'For each specified hypothesis the posterior probability is output. If the
#'hypotheses are not exhaustive (i.e., do not cover the entire parameter space)
#'this includes the posterior probability of the complement to the input
#'hypotheses. The complement is the hypothesis that neither of the input
#'hypotheses is true. For example, inputting \dQuote{X1 > 0; X1 < 0} gives
#'posterior probabilities for only for these hypotheses, whereas inputting
#'\dQuote{(X1, X2) > 0} gives posterior probabilities for \dQuote{(X1, X2) > 0}
#'and \dQuote{not (X1, X2) > 0}.
#'
#'By saving the test as an object it is also possible to access the
#'\code{BF_matrix} which compares the hypotheses directly against each other
#'(see examples). This matrix divides the column hypothesis by each row
#'hypothesis and can be interpreted as \dQuote{given the data, [column
#'hypothesis] is [value] times as likely as [row hypothesis]}.
#'
#'@section References: Mulder, J. (2014). Prior adjusted default Bayes factors
#'  for testing (in) equality constrained hypotheses. Computational Statistics &
#'  Data Analysis, 71, 448-463.
#'
#'@examples
#'###Standardize variables and fit the linear model
#'dt <- as.data.frame(scale(mtcars[, c(1, 3:4, 6)]))
#'fit <- lm(mpg ~ disp + hp + wt, data = dt)
#'
#'###Define hypotheses based on theory and test them
#'hyp <- "(wt, hp) > disp > 0; (wt, hp) > disp = 0"
#'res <- test_hyp(fit, hyp)
#'res
#'
#'###Examine output that compares hypotheses directly with Bayes factors
#'res$BF_matrix
#'
#'@param object A regression model object fit using the \code{lm} function.
#'@param hyp A string specifying hypotheses to be tested using the variable
#'  names of the \code{lm} object.
#'@param mcrep Integer specifying the number of iterations if no analytical
#'  solutions is possible. This is rare and only the case if the rank of the
#'  constraint matrix is less than its number of rows.
#'
#'@export test_hyp

test_hyp <- function(object, hyp, mcrep = 1e6){

  #1) initial setup and checks of input----
  varnames <- variable.names(object)
  if(is.null(varnames)) stop("Please input proper linear model object")
  betahat <- object$coefficients

  k <- length(varnames)
  n <- length(object$fitted.values)
  b <- (k + 1) / n

  hyp2 <- gsub(" ", "", hyp)
  if(!grepl("^[0-9a-zA-Z><=,;().-]+$", hyp2)) stop("Impermissable characters in hypotheses.")
  if(grepl("[><=]{2,}", hyp2)) stop("Do not use combined comparison signs e.g., '>=' or '=='")

  step1 <- unlist(strsplit(hyp2, split = "[<>=,;()]"))
  input_vars <- step1[grep("[a-zA-Z]+", step1)]
  if(!all(input_vars %in% varnames)) stop("Hypothesis variable(s) not in object, check spelling")

  hyp <- unlist(strsplit(hyp2, split = ";"))
  for(no in seq_along(hyp)){names(hyp)[no] <- paste0("H", no)}
  hyps <- unlist(strsplit(hyp2, split = ";"))
  BFu <- rep(NA, length = length(hyps))

  BFip_posterior <- if(any(!grepl("=", hyps))) {rep(NA, sum(!grepl("=", hyps)))} else{NULL}
  if(!is.null(BFip_posterior)) {
    R_i_all <- vector("list", length =  length(BFip_posterior))
    r_i_all <- vector("list", length =  length(BFip_posterior))
    ineq_marker <- 0
    BFip_prior <- BFip_posterior
  }

  for(h in seq_along(hyps)){
    hyp2 <- hyps[[h]]
    step2 <- unlist(strsplit(hyp2, split = "[<>=,()]"))
    hyp_vars <- step2[grep("[a-zA-Z]+", step2)]
    if(any(duplicated(hyp_vars))) stop("Variables should occur only once in a hypothesis. Check semicolons.")

    #2)hyp-to-matrices----
    framer <- function(x){
      pos_comparisons <- unlist(gregexpr("[<>=]", x))
      leftside <- rep(NA, length(pos_comparisons) + 1)
      rightside <- rep(NA, length(pos_comparisons) + 1)
      pos1 <- c(-1, pos_comparisons)
      pos2 <- c(pos_comparisons, nchar(x) + 1)
      for(i in seq_along(pos1)){
        leftside[i] <- substring(x, pos1[i] + 1, pos1[i+1] - 1)
        rightside[i] <- substring(x, pos2[i] + 1, pos2[i+1] - 1)
      }
      leftside <- leftside[-length(leftside)]
      rightside <- rightside[-length(rightside)]
      comparisons <- substring(x, pos_comparisons, pos_comparisons)
      data.frame(left = leftside, comp = comparisons, right = rightside, stringsAsFactors = FALSE)
    }

    framed <- framer(hyp2)

    if(any(grepl(",", framed$left)) || any(grepl(",", framed$right))){
      if(nrow(framed) > 1){
        for(r in 1:(nrow(framed)-1)){
          if(all.equal(framed$right[r], framed$left[r+1])){
            if(substring(framed$right[r], 1, 1) == "(") {
              framed$right[r] <- sub("),.+", ")", framed$right[r])
              framed$left[r+1] <- sub(".+),", "", framed$left[r +1])
              } else{
                framed$right[r] <- sub(",.+", "", framed$right[r])
                framed$left[r+1] <- sub("[^,]+,", "", framed$left[r+1])
              }
            }
          }
        }

      commas_left <- framed$left[grep(",", framed$left)]
      commas_right <- framed$right[grep(",", framed$right)]
      if(isTRUE(any(!grepl("\\(.+)", commas_left))) || isTRUE(any(!grepl("\\(.+)", commas_right))) ||
         isTRUE(any(grepl(").+", commas_left))) || isTRUE(any(grepl(").+", commas_right))) ||
         isTRUE(any(grepl(".+\\(", commas_left))) || isTRUE(any(grepl(".+\\(", commas_right)))) {
        stop("Incorrect hypothesis syntax or extra character, check specification")
        }

      framed$left <- gsub("[()]", "", framed$left)
      framed$right <- gsub("[()]", "", framed$right)
      commas <- unique(c(grep(",", framed$left), grep(",", framed$right)))

      if(length(commas) > 0){
        multiples <- vector("list", length = length(commas))

        for(r in seq_along(commas)){
          several <- framed[commas,][r, ]

          if(several$comp == "="){

            several <- c(several$left, several$right)
            separate <- unlist(strsplit(several, split = ","))
            if(any(grepl("^$", several))) stop("Misplaced comma in hypothesis")
            converted_equality <- paste(separate, collapse = "=")
            multiples[[r]] <- framer(converted_equality)

            } else{
            leftvars <- unlist(strsplit(several$left, split = ","))
            rightvars <- unlist(strsplit(several$right, split = ","))
            if(any(grepl("^$", leftvars)) || any(grepl("^$", rightvars))) stop("Misplaced comma in hypothesis")

            left <- rep(leftvars, length.out = length(rightvars)*length(leftvars))
            right <- rep(rightvars, each = length(leftvars))
            comp <- rep(several$comp, length(left))

            multiples[[r]] <- data.frame(left = left, comp = comp, right = right, stringsAsFactors = FALSE)
            }
          }

        framed <- framed[-commas,]
        multiples <- do.call(rbind, multiples)
        framed <- rbind(multiples, framed)
      }
    }

    equality <- framed[framed$comp == "=",]
    inequality <- framed[!framed$comp == "=",]

    #****Equality part string-to-matrix
    if(nrow(equality) == 0) {
      R_e <- r_e <- NULL
    } else{
      outcomes <- suppressWarnings(apply(equality[, -2], 2, as.numeric))
      outcomes <- matrix(outcomes, ncol = 2, byrow = FALSE)
      if(any(rowSums(is.na(outcomes)) == 0)) stop("Value compared with value rather than variable, e.g., '2 = 2', check hypotheses")
      rows <- which(rowSums(is.na(outcomes)) < 2)
      specified <- t(outcomes[rows,])
      specified <- specified[!is.na(specified)]
      r_e <- ifelse(rowSums(is.na(outcomes)) == 2, 0, specified)
      r_e <- matrix(r_e)

      var_locations <- apply(equality[, -2], 2, function(x) ifelse(x %in% varnames, match(x, varnames), 0))
      var_locations <- matrix(var_locations, ncol = 2)

      R_e <- matrix(rep(0, nrow(equality)*length(varnames)), ncol = length(varnames))

      for(i in seq_along(r_e)){
        if(!all(var_locations[i, ] > 0)){
          R_e[i, var_locations[i,]] <- 1
        } else{
          R_e[i, var_locations[i,]] <- c(1, -1)
        }
      }
    }


    #****Inequality part string-to-matrix
  if(nrow(inequality) == 0) {
    R_i <- r_i <- NULL
  } else{
    outcomes <- suppressWarnings(apply(inequality[, -2], 2, as.numeric))
    outcomes <- matrix(outcomes, ncol = 2, byrow = FALSE)
    if(any(rowSums(is.na(outcomes)) == 0)) stop("Value compared with value rather than variable, e.g., '2 > 2', check hypotheses")
    cols <- which(rowSums(is.na(outcomes)) < 2)
    specified <- t(outcomes[cols,])
    specified <- specified[!is.na(specified)]
    r_i <- ifelse(rowSums(is.na(outcomes)) == 2, 0, specified)
    r_i <- matrix(r_i)

    leq <- which(inequality$comp == "<")
    var_locations <- apply(inequality[, -2], 2, function(x) ifelse(x %in% varnames, match(x, varnames), 0))
    var_locations <- matrix(var_locations, ncol = 2)

    R_i <- matrix(rep(0, nrow(inequality)*length(varnames)), ncol = length(varnames))

    for(i in seq_along(r_i)){
      if(!all(var_locations[i, ] > 0)){

        if(var_locations[i, 1] == 0){ 
          if(i %in% leq){
            value <-  1  
          } else{ 
            r_i[i] <- r_i[i]*-1 
            value <- -1 
          }
        } else{ 
          if(i %in% leq){ 
            r_i[i] <- r_i[i]*-1 
            value <-  -1  
          } else{
            value <- 1 
          }
        }
        
        R_i[i, var_locations[i,]] <- value 
        
      } else{
        value <- if(i %in% leq) c(-1, 1) else c(1, -1)
        R_i[i, var_locations[i,]] <- value
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
    R_ei <- rbind(R_e,R_i)
    r_ei <- rbind(r_e,r_i)
    Rr_ei <- cbind(R_ei,r_ei)
    beta_zero <- MASS::ginv(R_ei)%*%r_ei

    if(nrow(Rr_ei) > 1){
      rref_ei <- pracma::rref(Rr_ei)
      nonzero <- rref_ei[,k+1]!=0
      if(max(nonzero)>0){
        row1 <- max(which(nonzero==T))
        if(sum(abs(rref_ei[row1,1:k]))==0){
          stop("Default prior mean cannot be constructed from constraints.")
        }
      }
    }

    if(comparisons == "only equality"){
      #**3.1)only-equality----

      delta <- R_e %*% betahat
      delta_zero <- R_e %*% beta_zero

      #Scale matrix components
      X <- model.matrix(object)
      RX <- R_e %*% solve((t(X) %*% X)) %*% t(R_e)
      s2 <- sum((model.frame(object)[[1]] - X %*% betahat)^2)

      #Scale matrix for posterior t-distribution
      scale_post <- matrix(s2 * RX / (n - k), ncol = nrow(R_e))

      #Scale matrix for prior t-distribution
      scale_prior <- matrix(s2 * RX / (n*b - k), nrow(R_e))

      #Hypothesis test
      log_BF <- mvtnorm::dmvt(x = t(r_e), delta = delta, sigma = scale_post, df = n - k, log = TRUE) -
        mvtnorm::dmvt(x = t(r_e), delta = delta_zero, sigma = scale_prior, df = n*b - k, log = TRUE)

      BF <- exp(log_BF)

    } else if(comparisons == "only inequality"){
      #**3.2)only-inequality----

      ineq_marker <- ineq_marker + 1
      r_i <- as.vector(r_i)

          if(Matrix::rankMatrix(R_i)[[1]] == nrow(R_i)){

          delta <- as.vector(R_i %*% betahat)
          delta_zero <- as.vector(R_i %*% beta_zero)

          #Scale matrix components
          X <- model.matrix(object)
          RX <- R_i %*% solve(t(X) %*% X) %*% t(R_i)
          s2 <- sum((model.frame(object)[[1]] - X %*% betahat)^2)

          #Scale matrix for posterior t-distribution
          scale_post <- s2 * RX / (n - k)

          #Scale matrix for prior t-distribution
          scale_prior <- s2 * RX / (n*b - k)

            if(nrow(scale_post) == 1){ #If univariate
              prior_prob <- pt((r_i - delta_zero) / sqrt(scale_prior), df = n*b - k, lower.tail = FALSE)[1]
              posterior_prob <- pt((r_i - delta) / sqrt(scale_post), df = n - k, lower.tail = FALSE)[1]
              BF <- posterior_prob / prior_prob
            } else { #If multivariate
              prior_prob <- mvtnorm::pmvt(lower = r_i, upper = Inf, delta = delta_zero, sigma = scale_prior, df = n*b - k, type = "shifted")[1]
              posterior_prob <- mvtnorm::pmvt(lower = r_i, upper = Inf, delta = delta, sigma = scale_post, df = n - k, type = "shifted")[1]
              BF <- posterior_prob / prior_prob
            }

          } else{ #No transformation is possible.
            if(!is.numeric(mcrep) || !mcrep %% 1 == 0) stop("Input for mcrep should be an integer")

            #Scale matrix for posterior t-distribution.
            scale_post <- vcov(object)

            #Scale matrix for prior t-distribution
            scale_prior <- vcov(object) * (n - k) / (n*b - k)

            draws_post <- mvtnorm::rmvt(n = mcrep, delta = betahat, sigma = scale_post, df = n - k)
            draws_prior <- mvtnorm::rmvt(n = mcrep, delta = beta_zero, sigma = scale_prior, df = n*b - k)

            prior_prob <- mean(apply(draws_prior%*%t(R_i) > rep(1, mcrep)%*%t(r_i), 1, prod))
            posterior_prob <- mean(apply(draws_post%*%t(R_i) > rep(1, mcrep)%*%t(r_i), 1, prod))
            BF <- posterior_prob / prior_prob

          }

      BFip_prior[ineq_marker] <- prior_prob
      BFip_posterior[ineq_marker] <- posterior_prob
      R_i_all[[ineq_marker]] <- R_i
      r_i_all[[ineq_marker]] <- matrix(r_i)


    } else{

      #**3.3)both-comparisons----

      #****Equality
      q_e <- nrow(R_e)

      #Scale matrix for posterior t-distribution
      scale_post <- vcov(object)

      #Scale matrix for prior t-distribution
      scale_prior <- vcov(object) * (n - k) / (n*b - k)

      #a)Transformation matrix
      D <- diag(k) - t(R_e) %*% solve(R_e %*% t(R_e)) %*% R_e
      D2 <- unique(round(D, 5))
      D2 <- D2[as.logical(rowSums(D2 != 0)),]
      Tm <- rbind(R_e, D2)

      #b)
      w_post <- Tm %*% betahat
      w_prior <- Tm %*% beta_zero
      K_post <- Tm %*% scale_post %*% t(Tm)
      K_prior <- Tm %*% scale_prior %*% t(Tm)

      #Equality BF
      log_BF <- mvtnorm::dmvt(x = t(r_e), delta = w_post[1:q_e], sigma = matrix(K_post[1:q_e, 1:q_e], ncol = q_e), df = n - k, log = TRUE) -
        mvtnorm::dmvt(x = t(r_e), delta = w_prior[1:q_e], sigma = matrix(K_prior[1:q_e, 1:q_e], ncol = q_e), df = n*b - k, log = TRUE)

      BFe <- exp(log_BF)

      #****Inequality
      R_iv <- R_i %*% MASS::ginv(D2)
      r_iv <- r_i - R_i %*% MASS::ginv(R_e) %*% r_e

      ##Partitioning
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
      w_2g1_post <- w_2_post + K_21_post %*% solve(K_11_post) %*% matrix(r_e - w_1_post)
      K_2g1_post <- as.vector((n - k + (t(matrix(r_e - w_1_post)) %*% solve(K_11_post) %*% matrix(r_e - w_1_post))) /
                                (n - k + q_e)) * (K_22_post - K_21_post %*% solve(K_11_post) %*% t(K_21_post))

      #prior
      w_2g1_prior <- w_2_prior + K_21_prior %*% solve(K_11_prior) %*% matrix(r_e - w_1_prior)
      K_2g1_prior <- as.vector((n*b - k + (t(matrix(r_e - w_1_prior)) %*% solve(K_11_prior) %*% matrix(r_e - w_1_prior))) /
                                 (n*b - k + q_e)) * (K_22_prior - K_21_prior %*% solve(K_11_prior) %*% t(K_21_prior))

      if(Matrix::rankMatrix(R_iv)[[1]] == nrow(R_iv)){
        r_iv <- as.vector(r_iv)

        delta_post <- as.vector(R_iv %*% w_2g1_post)
        delta_prior <- as.vector(R_iv %*% w_2g1_prior)

        scale_post_trans <- R_iv %*% K_2g1_post %*% t(R_iv)
        scale_prior_trans <- R_iv %*% K_2g1_prior %*% t(R_iv)

        if(nrow(scale_post_trans) == 1){ #If univariate
          BFi <- pt((r_iv - delta_post) / sqrt(scale_post_trans), df = n - k + q_e, lower.tail = FALSE)[1] /
            pt((r_iv - delta_prior) / sqrt(scale_prior_trans), df = n*b - k + q_e, lower.tail = FALSE)[1]
          } else { #if multivariate
            BFi <- mvtnorm::pmvt(lower = r_iv, upper = Inf, delta = delta_post, sigma = scale_post_trans, df = n - k + q_e, type = "shifted")[1] /
              mvtnorm::pmvt(lower = r_iv, upper = Inf, delta = delta_prior, sigma = scale_prior_trans, df = n*b - k + q_e, type = "shifted")[1]
            }

        } else{
          if(!is.numeric(mcrep) || !mcrep %% 1 == 0) stop("Input for mcrep should be an integer")

          #Draw from prior and posterior
          draws_post <- mvtnorm::rmvt(n = mcrep, delta = w_2g1_post, sigma = K_2g1_post, df = n - k + q_e)
          draws_prior <- mvtnorm::rmvt(n = mcrep, delta = w_2g1_prior, sigma = K_2g1_prior, df = n*b - k + q_e)

          BFi <- mean(apply(draws_post%*%t(R_iv) > rep(1,mcrep)%*%t(r_iv),1,prod)) /
            mean(apply(draws_prior%*%t(R_iv) > rep(1,mcrep)%*%t(r_iv),1,prod))

        }

      #Total BF
      BF <- BFe * BFi

      } #end 'both comparisons' option

    BFu[h] <- BF #hypothesis vs. unconstrained
    names(BFu)[[h]] <- paste0("H", h)

  } #end loop over all hypotheses

  if(!is.null(BFip_posterior)){
    if(length(BFip_posterior) == 1){
     BFc <- (1 - BFip_posterior) / (1 - BFip_prior)
    } else{
      R_i_overlap <- do.call(rbind, R_i_all)
      r_i_overlap <- do.call(rbind, r_i_all)

      ineq_draws_prior <- mvtnorm::rmvt(n = 1e4, delta = beta_zero, sigma = vcov(object) * (n - k) / (n*b - k), df = (n*b - k))
      exhaustive <- mean(rowSums(ineq_draws_prior%*%t(R_i_overlap) > rep(1, 1e4)%*%t(r_i_overlap)) > 0)

      if(exhaustive == 1){
        BFc <- NULL
      } else{
        overlap <- mean(apply(ineq_draws_prior%*%t(R_i_overlap) > rep(1, 1e4)%*%t(r_i_overlap), 1, prod))

        if(overlap == 0){
          BFc <-  (1 - sum(BFip_posterior)) / (1 - sum(BFip_prior))
        } else{
          ineq_draws_posterior <- mvtnorm::rmvt(n = 1e4, delta = betahat, sigma = vcov(object), df = n - k)

          constraints_prior <- Map(function(Ri, ri){apply(ineq_draws_prior%*%t(Ri) > rep(1,1e4)%*%t(ri), 1, prod)}, R_i_all, r_i_all)
          constraints_posterior <- Map(function(Ri, ri){apply(ineq_draws_posterior%*%t(Ri) > rep(1,1e4)%*%t(ri), 1, prod)}, R_i_all, r_i_all)

          prob_prior <- mean(Reduce(`+`, constraints_prior) > 0)
          prob_posterior <- mean(Reduce(`+`, constraints_posterior) > 0)
          BFc <- (1 - prob_posterior) / (1 - prob_prior)
        }
      }
    }
  } else{
    BFc <- 1
  }

  if(!is.null(BFc)){names(BFc) <- "Hc"}
  BFu <- c(BFu, BFc)
  out_hyp_prob <- BFu / sum(BFu)

  BF_matrix <- matrix(rep(BFu, length(BFu)), ncol = length(BFu), byrow = TRUE)
  BF_matrix <- BF_matrix / BFu
  colnames(BF_matrix) <- rownames(BF_matrix) <- names(BFu)

  out <- list(BF_matrix_rounded = round(BF_matrix, digits = 3), BF_matrix_unrounded = BF_matrix,
              post_prob = out_hyp_prob, hypotheses = hyp)
  class(out) <- "hyp"
  out

}
#End----
