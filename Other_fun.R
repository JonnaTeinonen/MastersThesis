# The Script includes supportive functions for the MI methods.

library(ResourceSelection) # For hoslem.test()


### AKAIKE-WEIGHTS FUNCTION ###

# Functions: The function calculates the Akaike weights of the models based
# on the given AIC scores.

# Parameters:
# 1. AIC_scores: A vector including the AIC scores of all the models

# Returns: A list that includes:
# 1. A vector including the Akaike weights for all the models
# 2. Error count: adds one if one of the Akaike weights is 0

Akaike_w <- function(AIC_scores){
  NA_weight_count <- 0 # Error vount
  AIC_diff <- rep(0, length(AIC_scores))
  A_weights <- rep(0, length(AIC_scores))
  
  # We calculate Akaike weight for all the outcome models: weight will be same 
  # for each category within the model
  
  # First the differences between each AIC score of the model and the smallest
  # AIC-score of all the models is computed
  for (i in seq_along(AIC_scores)){
    
    AIC_diff[i] <- AIC_scores[i] - min(AIC_scores)
  }
  
  # The differences are used to compute the Akaike weights
  for (i in seq_along(AIC_diff)){
    A_weights[i] <- exp((-1/2) * AIC_diff[i]) / (sum(exp((-1/2) * AIC_diff)))
  }

  # Error ocunt: adds 1 if any of the Akaike weights is 0
  if(any(is.na(A_weights))){
    NA_weight_count <- NA_weight_count + 1
  }
  
  # Returns Akaike weights and error count
  return(list(A_weights, NA_weight_count))
}



### HL -STATISTIC FUNCTION ###

# Functions: The function the HL chi-square statistic for each model

# Parameters:
# 1. C: number of categories of Y
# 2. outcome_data: a data frame including the y values of units j in the
#    bootstrap and C-1 predicted outcome values from each outcome model.
# 3.combined_response_model: a model in which response indicator is fitted 
#                           against all predictive response values obtained 
#                           from response models

# Returns: A list that includes:
# 1. A vector in which the first C-1 values are the HL-statistic of the
#    C-1 outcome values and the final value is the HL statistic of the 
#   response value.
#2. Error count

HL_fun <- function(C, outcome_data, combined_response_model){
  NA_weight_count <- 0
  
  ### STEP 1: Calculate HL statistic for C-1 predictive outcome scores
  
  # Empty vector in which the HL chi-square for C-1 categories are saved
  HL_outcome <- rep(0, (C-1))
  
  # We loop through each C-1 category
  for(i in 1:(C -1)){
    # Indices of units belonging to category c or to 1
    indices <- which(outcome_data[, 1] %in% c(1, i+1))
    
    # Only choose units that belong to category c or 1
    data_c <- outcome_data[indices, ]
    
    # If unit j belong to category c, change Y value to 1, else to 0
    data_c[, 1] <- ifelse(data_c[, 1] == (i+1), 1, 0)
    
    
    # Regress the category c on all explanatory variables (predictive
    # scores from all outcome models)
    logreg_c <- glm(y ~ ., data = data_c, 
                    family = binomial,
                    control = list(maxit = 100)) 
    
    # Calculates HL chi-square by dividing the units into 10 groups
    HL_c <- hoslem.test(logreg_c$y, fitted(logreg_c), g = 10)$statistic
    HL_outcome[i] <- HL_c
    
  }
  
  ### STEP 2: Calculate HL chi-square for R by dividing the units into 10 groups
  HL_R <-  hoslem.test(combined_response_model$y, 
                       fitted(combined_response_model), g = 10)$statistic
  
  ### STEP 3: Return all the calculated HL chi-square statistics in one vector
  HL_stat <- c(HL_outcome, HL_R)
  
  # Returns HL statistics and error ocunt
  return(list(HL_stat, NA_weight_count))
}



### EUCLIDEAN DISTANCE FUNCTION ###

# Functions: Calculates the Euclidean distance between unit i with missing Y 
# from the original data set and unit j with observed Y from bootstrap sample. 

# Parameters:
# 1. pred_i: a matrice including the final predicted outcome and response
#            values for units with missing Y from original sample.
# 2. pred_j: a matrice including the final predicted outcome and response
#            values for units with observed Y in the bootstrap sample.
# 3. weights: a vector including first the weights for the C-1 outcome values
#             and weight for the response value for each unit,

# Returns: an object in which each column represents one unit and the rows
#   represent the distances between unit i and each unit j.

fun_distance <- function(pred_i, pred_j, weights){
  all_distances <- sqrt(sum(weights * (pred_i - pred_j)^2))
  return(all_distances)
}



### POOLED PARAMETER FUNCTION ###

# Functions: Calculates the pooled average estimate/mean, standard error 
# of the estimate, 95% CI and if the coverage rate.

# Parameters:
# 1. MI_sets: Object including M imputed data sets.
# 2. population_proportions: A vector of the true population proportions of 
#    C categories of Y.

# Returns: a data frame in which each row represents category c of Y and 
# columns represent the pooled parameter estimates.


fun_pooled_param <- function(MI_sets, population_proportions){
  M <- length(MI_sets) # Number of imputed data sets
  n <- length(MI_sets[[1]]$y) # Number of units
  C <- length(levels(MI_sets[[1]]$y)) # Number of Y categories
  
  # Calculates the percentage of observations in each category of Y
  p_c <- sapply(MI_sets, function(x) table(x$y)/length(x$y))
  
  # Calculates the mean percentage of observations per category c over M
  # imputations
  mean_p_c <- apply(p_c, 1, mean)
  
  # Calculates the pooled standard error
  SE_pooled <- apply(p_c, 1, function(x){
    v_w <- (1/M) * sum(x * (1 - x) / n) # Within imputation variance
    v_b <-  sum((x- mean(x))^2) / (M - 1) # Between imputation variance
    v_t <- v_w + v_b + (v_b / M) # Total imputation variance
    SE <- sqrt(v_t) # Pooled SE
    return(SE)
  })
  
  
  # Calculates 95% CI
  CI_95 <- c(mean_p_c - qnorm(0.975) * SE_pooled, 
             mean_p_c + qnorm(0.975) * SE_pooled)
  
  # Checks if the true population proportions for C are within the 95% CI
  # (TRUE/FALSE)
  within_CI <- CI_95[1:C] <= population_proportions  & 
    population_proportions  <= CI_95[c(C + 1):length(CI_95)]
  
  
  # Creates the final output of the pooled estimates
  output <- data.frame(Est = mean_p_c,
                       SE = SE_pooled,
                       CI_lower = CI_95[1:C],
                       CI_upper = CI_95[c(C + 1):length(CI_95)],
                       within_CI = within_CI)
  return(output) 
}




