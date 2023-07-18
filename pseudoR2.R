# Script of the supportive functions for MRIC Pseudo-R2.



### MCFADDEN PSEUDO-R2 FUNCTION ###
# Functions: calculates the McFadden pseudo-R2 for the predictive outcome and
# response values.
# Parameters:
# 1. outcome_model: model combining outcome values from all outcome models
# 2. outcome_data: data used in the outcome_model
# 3. response_model: model combining response values from all response models
# 4. response_data: data used in the response_model

# Returns: A list than includes:
# 1. A vector that first has the pseudo-R2s for the C -1 predictive  outcome 
#   values and the pseudo-R2 for the predictive response value
# 2. Error count

R2_MF <- function(outcome_model, outcome_data, response_model, response_data){
  
  NA_weight_count <- 0 # Placeholder for error count
  
  ### STEP 1: Calculate McFadden pseudo R2 for C-1 predictive outcome scores 
  Y <- outcome_data$y # All Y values
  C <- ncol(outcome_model$fitted.values) # Number of categories of Y
  R <- response_data$R # Response indicator for each unit
  
  # Empty vector in which the pseudo-R2's for outcome values are saved
  R2_outcome <- rep(0, (C-1))
  probs_C <- fitted(outcome_model) # Prob. of belonging to categories of Y 
  
  # We loop through each C-1 category (c = 2...C)
  for(i in 1:(C -1)){

    # Indices of units belonging to category c or 1
    indices <- which(Y %in% c(1, i+1))
    
    # Probability of unit belonging to category c vs. 1
    p_outcome_c <- probs_C[ , i + 1]/(probs_C[, 1] + probs_C[ , i+1])

    # Chooses only the probabilities of units that belong to category c or 1
    p_outcome_c <- p_outcome_c[indices]

    # Chooses only Y values of units belonging to category c or 1
    y_c <- Y[indices]
    
    # If unit j belongs to category c, change Y value to 1, else to 0
    y_c <- ifelse(y_c == (i+1), 1, 0)
    
    # Calculates loglikelihood for estimated model (all predictors)
    loglik_alt <- sum(y_c * log(p_outcome_c) + (1 - y_c) * log(1 - p_outcome_c),
                      na.rm = T)
    
    # Calculates loglikelihood for the null model (only intercept)
    loglik_null <- sum(y_c * log(mean(y_c)) 
                       + (1 - y_c) * log(1 - mean(y_c)),
                       na.rm = T)
    
    # Calculates MCFadden's R2 for category c using the formula
    R2_c <- 1 - loglik_alt / loglik_null
    R2_outcome[i] <- R2_c
    
  }
  
  ### STEP 2: Calculate McFadden pseudo-R2 for the response model
  p_missing <- fitted(response_model) # Probability that Y is observed
  loglik_alt <- sum(R * log(p_missing) + (1 - R) * log(1 - p_missing),
                    na.rm = T)
  loglik_null <- sum(R * log(mean(R)) + (1 - R) * log(1 - mean(R)),
                     na.rm = T)
  
  
  R2_R <- 1 - loglik_alt / loglik_null 
  
  ### STEP 3: Return all the calculated pseudo-R2's in one vector
  R2_weights <- c(R2_outcome, R2_R)
  
  # error count: add 1 If one of the computed pseudo-R2s is 0
  if(any(is.na(R2_weights))){
    NA_weight_count <- NA_weight_count + 1
  }
  
  # Returns the pseudo-R2s and the error count
  return(list(R2_weights, NA_weight_count))
  
}



### COX AND SNELL PSEUDO-R2 FUNCTION ###
# Functions: calculates the Cox and Snell pseudo-R2 for the predictive outcome 
# and the response values.

# Parameters:
# 1. outcome_model: model combining outcome values from all outcome models
# 2. outcome_data: data used in the outcome_model
# 3. response_model: model combining response values from all response models
# 4. response_data: data used in the response_model

# Returns: A list than includes:
# 1. A vector that first has the pseudo-R2s for the C -1 predictive outcome 
#   values and the pseudo-R2 for the predictive response value
# 2. Error count

R2_CS <- function(outcome_model, outcome_data, response_model, response_data){
  
  NA_weight_count <- 0 # Placeholder for error count
  
  ### STEP 1: Calculate Cox and Snell pseudo-R2 for categories C-1 of Y 
  Y <- outcome_data$y # All Y values
  C <- ncol(outcome_model$fitted.values) # Number of categories of Y
  R <- response_data$R # Response indicator for each unit
  
  # Empty vector in which the pseudo R2's for C-1 outcome values are saved
  R2_outcome <- rep(0, (C-1))
  probs_C <- fitted(outcome_model) # Probabilities of belonging to each category
  
  # We loop through each C-1 category (c = 2 ... C)
  for(i in 1:(C -1)){
    # Indices of units belonging to category c or 1
    indices <- which(Y %in% c(1, i+1))
    
    # Probability of unit belonging to category c vs. 1
    p_outcome_c <- probs_C[ , i + 1]/(probs_C[, 1] + probs_C[ , i+1]) 
    
    # Chooses only the probabilities of units that belong to category c or 1
    p_outcome_c <- p_outcome_c[indices]
    
    # Chooses only Y values of units belonging to category c or 1
    y_c <- Y[indices]
    
    # If unit j belong to category c, change Y value to 1, else to 0
    y_c <- ifelse(y_c == (i+1), 1, 0)
    
    # Calculates the loglikelihood for estimated model (all predictors)
    loglik_alt <- sum(y_c * log(p_outcome_c) 
                      + (1 - y_c) * log(1 - p_outcome_c),
                      na.rm = T)
    
    # Calculates the loglikelihood for the null model (only intercept)
    loglik_null <- sum(y_c * log(mean(y_c)) 
                       + (1 - y_c) * log(1 - mean(y_c)),
                       na.rm = T)
    
    # Calculates Cox and Snell's R2 for category c using the formula
    R2_c <- 1 - exp(2 * (loglik_null - loglik_alt) / length(y_c))
    R2_outcome[i] <- R2_c
    
  }
  
  ### STEP 2: Calculate McFadden pseudo-R2 for the response values
  p_missing <- fitted(response_model) # Probability that Y is observed
  loglik_alt <- sum(R * log(p_missing) + (1 - R) * log(1 - p_missing),
                    na.rm = T)
  loglik_null <- sum(R * log(mean(R)) + (1 - R) * log(1 - mean(R)),
                     na.rm = T)
  
  R2_R <- 1 - exp(2 * (loglik_null - loglik_alt) / length(R))
  
  ### STEP 3: Return all the calculated pseudo-R2's in one vector
  R2_weights <- c(R2_outcome, R2_R)
  
  # Error ocunt: add 1 is one of the pseudo-R2s is 0
  if(any(is.na(R2_weights))){
    NA_weight_count <- NA_weight_count + 1
  }
  
  # Returns the oseudo-R2s and error count
  return(list(R2_weights, NA_weight_count))
}



### MCKELVEY AND ZAVOINA PSEUDO-R2 FUNCTION ###
# Functions: calculates the McKelvey and Zavoina pseudo-R2 for the predictive 
# outcome and the response values.

# Parameters:
# 1. outcome_model: model combining outcome values from all outcome models
# 2. outcome_data: data used in the outcome_model
# 3. response_model: model combining response values from all response models
# 4. response_data: data used in the response_model

# Returns: A list than includes:
# 1. A vector that first has the pseudo-R2s for the C -1 predictive outcome 
#   values and the pseudo-R2 for the predictive response value
# 2. Error count 1: How often one of the pseudo-R2s is 0
# 3. Error count 2: How often glm for pseudo-R2 did not converge

R2_MZ <- function(outcome_model, outcome_data, response_model, response_data){
  
  # Placeholder for error counts
  nonconv_count <- 0
  NA_weight_count <- 0
  
  ### STEP 1: Calculate McKelvey and Zavoina's R2 for each category C-1 of Y 
  C <- ncol(outcome_model$fitted.values) # Number of categories of Y
  
  # Empty vector in which the pseudo R2's for C-1 categories are saved
  R2_outcome <- rep(0, (C-1))
  probs_C <- fitted(outcome_model) # Probabilities of belonging to each category
  
  # We loop through each C-1 category (c = 2...C)
  for(i in 1:(C -1)){
    # Indices of units belonging to category c or to 1
    indices <- which(outcome_data[, 1] %in% c(1, i+1))
    
    # Only choose units that belong to category c or 1
    data_c <- outcome_data[indices, ]
    
    # If unit j belong to category c, change Y value to 1, else to 0
    data_c[, 1] <- ifelse(data_c[, 1] == (i+1), 1, 0)
    
    data_c$y <- as.factor(data_c$y)
    
    
    # Function that checks if glm converged when we regresses the category c on 
    # all predictive scores from all the outcome models.
    logreg_c <- fun_pseudoR2_convergence(data_c) 
    
    # If glm did not converge, pseudo-R2 would be replaced with fixed weights
    if(length(logreg_c) == 0) { 
      R2_c <- 1 / C
      R2_outcome[i] <- R2_c
      nonconv_count <- nonconv_count + 1 # adds 1 to the error count
      
    # If glm convergenced  
    } else {
      # Calculates McKelvey and Zavoina's R2 for category c using the formula
      R2_c <-  var(predict(logreg_c), na.rm = T) / 
        (var(predict(logreg_c), na.rm = T) + pi^2/3) 
      R2_outcome[i] <- R2_c
    }
  }
  
  ### STEP 2: Calculate McKelvey and Zavoina's R2 for the response value
  
  R2_R <- var(predict(response_model), na.rm = T) / 
    (var(predict(response_model), na.rm = T) + pi^2/3)
  
  ### STEP 3: Return all the calculated pseudo-R2's in one vector
  R2_weights <- c(R2_outcome, R2_R)
  
  # Error count: adds 1 if one of the pseudo-R2s is 0 
  if(any(is.na(R2_weights))){
    NA_weight_count <- NA_weight_count + 1
  }
  
  # Retuens a list with the pseudo-R2s and error counts
  return(list(R2_weights, NA_weight_count, nonconv_count))
}



### Nagelkerke PSEUDO-R2 FUNCTION ###
# Functions: calculates the Nagelkerke pseudo-R2 for the C-1 predictive 
# outcome and the response values.

# Parameters:
# 1. outcome_model: model combining predictive scores from all outcome models
# 2. outcome_data: data used in the outcome_model
# 3. response_model: model combining propensity scores from all response models
# 4. response_data: data used in the response_model

# Returns: A list than includes:
# 1. A vector that first has the pseudo-R2s for the C -1 predictive outcome 
#   values and the pseudo-R2 for the predictive response value
# 2. Error count: How often one of the pseudo-R2s is 0


R2_N <- function(outcome_model, outcome_data, response_model, response_data){
  
  ### STEP 1: Calculate Nagelkerke R2 for each category C-1 of Y 
  Y <- outcome_data$y
  C <- ncol(outcome_model$fitted.values) # Number of categories of Y
  R <- response_data$R
  
  NA_weight_count <- 0 # Placeholder for the error count
  
  # Empty vector in which the pseudo R2's for the outcome values are saved
  R2_outcome <- rep(0, (C-1))
  probs_C <- fitted(outcome_model) # Probabilities of belonging to each category
  
  
  # We loop through each C-1 category
  for(i in 1:(C -1)){
    # Indices of units belonging to category c or to 1
    indices <- which(Y %in% c(1, i+1))
    
    # Probability of unit belonging to category c vs. 1
    p_outcome_c <- probs_C[ , i + 1]/(probs_C[, 1] + probs_C[ , i+1]) 
    
    # Chooses only the probabilities of units that belong to category c or 1
    p_outcome_c <- p_outcome_c[indices]
    
    # Choose only Y values of units belonging to category c or 1
    y_c <- Y[indices]
    
    # If unit j belong to category c, change Y value to 1, else to 0
    y_c <- ifelse(y_c == (i+1), 1, 0)
    
    # Calculates loglikelihood for the estimated model (all predictors)
    loglik_alt <- sum(y_c * log(p_outcome_c) 
                      + (1 - y_c) * log(1 - p_outcome_c),
                      na.rm = T)
    
    # Calculates loglikelihood for the null model (only intercept)
    loglik_null <- sum(y_c * log(mean(y_c)) 
                       + (1 - y_c) * log(1 - mean(y_c)),
                       na.rm = T)
    
    # Calculates Nagelkerke's R2 for category c using the formula
    R2_c <-  (1 - exp(- (2 * (loglik_alt - loglik_null))/ length(y_c))) /
      (1 - exp(2 * loglik_null / length(y_c)))
    R2_outcome[i] <- R2_c
    
  }
  
  ### STEP 2: Calculate McFadden pseudo-R2 for the response values
  p_missing <- fitted(response_model) # Probability that Y is observed
  loglik_alt <- sum(R * log(p_missing) + (1 - R) * log(1 - p_missing),
                    na.rm = T)
  loglik_null <- sum(R * log(mean(R)) + (1 - R) * log(1 - mean(R)),
                     na.rm = T)
  
  R2_R <-  (1 - exp(- (2 * (loglik_alt - loglik_null))/ length(R))) /
    (1 - exp(2 * loglik_null / length(R)))
  
  ### STEP 3: Return all the calculated pseudo-R2's in one vector
  R2_weights <- c(R2_outcome, R2_R)
  
  # Error count: adds 1 if one of the pseudo-R2s is 0
  if(any(is.na(R2_weights))){
    NA_weight_count <- NA_weight_count + 1
  }
  
  # Returns a list including the pseudo-R2s and the error count
  return(list(R2_weights, NA_weight_count))
}



### FUNCTRION TO CHOOSE PSEUDO-R2 TYPE ###

# Functions: The function calls the corresponding pseudo-R2 type function.

# Parameters:
# 1. pseudoR2: Type of pseudo-R2 written in a string format ("McFadden", 
#             "CoxSnell","McKelveyZavoina", "Nagelkerke") 
# 2. outcome_model: model combining predictive scores from all outcome models
# 3. outcome_data: data used in the outcome_model
# 4. response_model: model combining propensity scores from all response models
# 5. response_data: data used in the response_model

# Returns: Returns the specific pseudo-R2s and error count(s)


R2_pseudo <- function(pseudoR2, outcome_model, outcome_data, response_model, 
                      response_data) {
  
  switch(pseudoR2,
         "McFadden" = return(R2_MF(outcome_model, outcome_data, 
                                   response_model, response_data)),
         
         "CoxSnell" = return(R2_CS(outcome_model, outcome_data, 
                                   response_model, response_data)),
         
         "McKelveyZavoina" = return(R2_MZ(outcome_model, outcome_data, 
                                          response_model, response_data)),
         
         "Nagelkerke" = return(R2_N(outcome_model, outcome_data, 
                                    response_model, response_data)))
}  



### PSEUDO-R2 ERROR TEST ###

# Functions: The function test that the glm can convergence when calculating 
# pseudo-R2 for the outcome values (only used with McKelvey and Zavoina)

# Parameters:
# 1. data: binary data used to calculate the pseudo-R2 for the outcome value

# Returns: If algorithm converges, it returns the binary regression. Otherwise
# it returns an empty list. 

fun_pseudoR2_convergence <- function(data) {
  tryCatch( expr = {
    logreg <- glm(y ~ ., 
                  data = data, 
                  family = binomial, 
                  control = list(maxit = 100))
    return(logreg)
  },
  
  error = function(e) { 
    logreg <- list()
    return(logreg )
  })
}
