
library(nnet)

####################################################################
#################### MRNMMI WITH FIXED WEIGHTS ##################### 
####################################################################

# Functions: MRNMMI with fixed weights 
# Reference: Breemer (2022): https://github.com/LigayaBreemer/MRNNMI 
#            File: PILOT_STUDY.R L:73-217
#            File: MRNMMI_fun.R: L: 4-166

# Description: The function performs MRNMMI M times by using fixed weights.

# Parameters:
# 1. sample_data: data frame of the sample data: includes variables "y" 
#                 (categorical with more than 2 categories), auxiliary 
#                 variables and response indicator "response" 
#                 (binomial: y = is observed, 0 = y is missing)
# 2. outcome_models: a list of the outcome models written in a string format 
#                   (e.g. "y ~ x1 + x2 + x4")
# 3. response_models: a list of the response models written in a string 
#                     format (e.g. "y ~ x1 + x2 + x4")
# 4. M: number of times data is imputed. In default set to 5.

# Returns: a list of the imputed data sets.

MRNNMI <- function(sample_data, outcome_models, response_models, M = 5){
  
  source("Other_fun.R") # Script that the includes supportive functions
  n <- nrow(sample_data) # Number of units in the sample
  C <- length(levels(sample_data$y)) # Number of categories of Y
  
  imputed_data <- list() # List of the imputed data sets
  
  # We impute M = 5 times to create M data sets
  for(m in 1:M){
    
    
    ### STEP 1: We take a bootstrap of size n 
    boot_sample <- sample_data[sample(1:n, n, replace = TRUE), ]
    
    # We check that all the categories are represented in the bootstrap sample.
    # If not, we take a new sample.
    while(length(unique(boot_sample[boot_sample$response == 1, "y"])) != 3){
      boot_sample <- sample_data[sample(1:n, n, replace = TRUE), ]
    }
    
    
    ### STEP 2: We give the outcome models, save the predicted outcome values 
    ### from each model and use them in a combined regression model as 
    ### predictors
    
    y_pred_all <- list() # List of all predictions from all outcome models
    all_outcome_models <- list() # List of all outcome models
    
    # Fits the models, saves them and the predicted values for C-1 
    # categories of y
    for(i in 1:length(outcome_models)){
      outcome_mod <- multinom(outcome_models[[i]], 
                              data = boot_sample,
                              maxit = 1000) # maximum number of iterations
      all_outcome_models[[i]] <- outcome_mod # saves the model
      y_pred_all[[i]] <-  fitted(outcome_mod)[ , -1] # predicted values
    }
    
    # Combines all the predicted values into a dataframe
    y_pred_all <- as.data.frame(do.call(cbind, y_pred_all))
    
    # Data frame containing units with observed y and their predicted values
    # NOTE: We only choose units in the bootstrap sample with observed Y.
    outcome_data <- as.data.frame(cbind(
      boot_sample[boot_sample$response == 1, "y"],  y_pred_all))
    colnames(outcome_data) <- c("y", paste0("V", 1:ncol(y_pred_all)))
    
    # Next we regress the outcome variable on the saved predicted values and 
    # save the new predicted values.
    combined_outcome_model <- multinom(y ~ ., 
                                       data = outcome_data,
                                       maxit = 1000)
    y_pred <- fitted(combined_outcome_model)[, -1] # The final predicted values
    
    # Standardize the final predicted values.
    mean_y_pred <- colMeans(y_pred) 
    sd_y_pred <- apply(y_pred, 2, function(x) sd(x))
    
    for(i in 1:ncol(y_pred)){
      y_pred[,i] <- (y_pred[,i] - mean_y_pred[i]) / sd_y_pred[i]
    }
    
    
    ### STEP 3: We give the response models, save the predicted response values 
    ### from each model and use them in the combined regression model 
    ### as predictors.
    
    R_pred_all <- list() # List of all response values
    all_response_models <- list() # List of all response models
    
    # Fits each response model, saves the model and the response values
    for(i in 1:length(response_models)){
      response_mod <- glm(response_models[[i]], 
                          data = boot_sample, 
                          family = "binomial",
                          control = list(maxit = 100))
      
      all_response_models[[i]] <- response_mod # saves the model
      R_pred_all[[i]] <-  fitted(response_mod) # predicted response values
    }
    
    # Combines all response values into a data frame
    R_pred_all <- as.data.frame(do.call(cbind, R_pred_all))
    
    # Data frame containing all unit and their predicted response values.
    # NOTE: We choose all the units in the bootstrap sample regardless if Y 
    # is missing or not.
    response_data <- as.data.frame(cbind(boot_sample$response, R_pred_all))
    colnames(response_data)<- c("R", paste0("V", 1:ncol(R_pred_all)))
    
    # Next we regress the "response" variable on the saved predicted response 
    # values and save the new response values. 
    combined_response_model <- glm(R ~ ., 
                                   data = response_data ,
                                   family = "binomial",
                                   control = list(maxit = 100))
    
    R_pred <- fitted(combined_response_model) # The final response values
    
    
    # Standardizes the final response values
    mean_R_pred <- mean(R_pred)
    sd_R_pred <- sd(R_pred) 
    R_pred <- (R_pred - mean_R_pred) / sd_R_pred
    
    
    ### STEP 4: We calculate the predicted outcome values and response values
    ### for the nonrespondents in the original sample and standardize these.
    
    # Selects all cases in original data set with missing Y.
    to_be_imputed <- sample_data[sample_data$response == 0, -1]
    
    # First we calculate the predicted outcome values:
    
    # List of all predicted values for units with missing Y
    missing_pred_y_all <- list() 
    
    # Calculates the C - 1 predicted values by using the outcome models for
    # the unit(s) in to_be_imputed
    for(i in 1:length(all_outcome_models)){
      
      if (nrow(to_be_imputed) == 1){ # If to_be_imputed only has one unit
        missing_pred_y_all[[i]] <- predict(all_outcome_models[[i]], 
                                           newdata = to_be_imputed, 
                                           type = "probs")[-1]
        
        
      } else { # if to_be_imputed has more than one unit
        missing_pred_y_all[[i]] <- predict(all_outcome_models[[i]], 
                                           newdata = to_be_imputed, 
                                           type = "probs")[, -1]
        
      }
    }
    
    # If to_be_imputed has only one unit
    if(nrow(to_be_imputed) == 1){
      # Combines all the predicted values into a data frame
      missing_pred_y_all <- as.data.frame(t(unlist(missing_pred_y_all)))
      colnames(missing_pred_y_all)  <- c(paste0("V", 
                                                1:ncol(missing_pred_y_all)))
      
      # Calculates the final outcome values
      missing_pred_y <- data.frame(t(predict(combined_outcome_model, 
                                             newdata = missing_pred_y_all, 
                                             type = "probs")[-1]))
      
      # Standardizes the predicted values
      for(i in 1:ncol(missing_pred_y )){
        missing_pred_y[i] <- (missing_pred_y[i] - mean_y_pred[i]) / sd_y_pred[i]
      }
      
    
    # if to_be_imputed has more than one unit    
    } else {
      # Combines all the predicted values into a data frame
      missing_pred_y_all <- as.data.frame(do.call(cbind, missing_pred_y_all))
      colnames(missing_pred_y_all)  <- c(paste0("V", 
                                                1:ncol(missing_pred_y_all)))
      
      # Calculates the final predicted values
      missing_pred_y <- predict(combined_outcome_model, 
                                newdata = missing_pred_y_all, 
                                type = "probs")[, -1]
      
      # Standardizes the predicted values
      for(i in 1:ncol(missing_pred_y )){
        missing_pred_y[,i] <- (missing_pred_y[,i] - mean_y_pred[i]) / 
          sd_y_pred[i]
      }
      
      
    }
    
    # Next we calculate predicted response values:
    
    # List of all predicted response values for units with missing Y
    missing_pred_R_all <- list()  
    
    # Calculates the response values by using the saved response models 
    for(i in 1:length(all_response_models)){
      missing_pred_R_all[[i]] <- predict(all_response_models[[i]], 
                                         newdata = to_be_imputed, 
                                         type = "response")
    }
    
    # Combines all the response values into a data frame
    missing_pred_R_all <- as.data.frame(do.call(cbind, missing_pred_R_all))
    
    # Calculates the final response values
    missing_pred_R <- predict(combined_response_model, 
                              newdata = missing_pred_R_all, 
                              type = "response")
    
    
    # Standardizes the final response vlaues
    missing_pred_R <- (missing_pred_R - mean_R_pred) / sd_R_pred
    
    
    ### STEP 5: We use an Euclidean distance function  to calculate the 
    ### similarity  between unit i with missing Y in the original data set and 
    ### unit j with observed Y in the bootstrap sample.
    
    # Creates a matrix with the final predicted outcome and response values
    # for units with missing Y in the original sample
    pred_to_be_imputed <- cbind(missing_pred_y, missing_pred_R) 
    
    # Creates a matrix with the final predicted outcome and response values
    # for units with observed Y in the bootstrap sample
    pred_boot_obs <- cbind(y_pred, R_pred[boot_sample$response == 1])
    
    
    # Saves the indices of the bootstrap units with observed Y.
    ind <- which(boot_sample$response == 1)
    
    # Defines the fixed weights for the Euclidean distance function
    w <- 1/ C
    
    # We calculate Euclidean distances. Returns a matrix in which each column 
    # contains the distances  from a unit i with missing Y in the original data 
    # set to all units j with observed Y in the bootstrap sample.
    distances <- apply(pred_to_be_imputed, 1, function(x) 
      apply(pred_boot_obs, 1, fun_distance, x, weights = w))
    
    
    ### STEP 6: For each unit with missing Y in the original data set, the
    ### imputation set is defined as k-Nearest Neighbor (in default k = 5)
    ### and the missing values of Y are imputed by randomly drawing a value of
    ### one donor from the k-Nearest Neighborhood. 
    
    # Ranks the distances from smallest to largest
    distance_rank <- apply(distances, 2, order)
    
    # Units are represent by each column. We get for each unit with missing Y 
    # (in the original data set) the indices of five other units with 
    # observed Y (in the bootstrap sample) with the smallest calculated distance.
    indices_smallest_distance <- apply(distance_rank, 2, function(x) ind[x<=5])
    
    # Creates a list of possible donors for each unit with missing Y 
    donors <- apply(indices_smallest_distance, 2, function(x){
      boot_sample[x, "y"]
    })
    
    # Randomly draws one donor for each unit with missing Y in the original 
    # data set.
    new_Y_values <- apply(donors, 2, function(x) sample(x, 1))
    
    # Imputes
    imputed_set <- sample_data
    imputed_set$y[imputed_set$response == 0] <- new_Y_values
    imputed_data[[m]] <- imputed_set
    
    
  }
  
  
  return(imputed_data) # Returns a list of the M imputed data sets
}



####################################################################
############################# DRNNMI ############################### 
####################################################################

# Functions: Performs DRNNMI with fixed weights M times.
# Reference: Chen et al. (2021) amd Breemer (2022): 
#            https://github.com/LigayaBreemer/MRNNMI 
#            File: PILOT_STUDY.R L:225-300

# Parameters:
# 1. sample_data: data frame of the sample data: includes variables "y" 
#                 (categorical with more than 2 categories), auxiliary 
#                 variables and response indicator "response" 
#                 (binomial: y = is observed, 0 = y is missing)
# 2. outcome_model: One outcome model in a string format 
#                   (e.g. "y ~ x1 + x2 + x4")
# 3. response_model: One response model in a string format
#                    (e.g. "y ~ x1 + x2 + x4")
# 4. M: number of times data is imputed. In default set to 5.

# Returns: a list of the imputed data sets.

DRNNMI <- function(sample_data, outcome_model, response_model, M = 5){
  
  source("Other_fun.R") # Script that the includes supportive functions
  n <- nrow(sample_data) # Number of units in the sample
  C <- length(levels(sample_data$y)) # Number of categories of Y
  
  imputed_data <- list() # List of the imputed data sets
  
  # We impute M = 5 times to create M data sets
  for(m in 1:M){
    
    ### STEP 1: We take a bootstrap sample
    boot_sample <- sample_data[sample(1:n, n, replace = TRUE), ]
    
    # We check that all the categories are represented in the bootstrap sample.
    # If not, we take a new sample.
    while(length(unique(boot_sample[boot_sample$response == 1, "y"])) != 3){
      boot_sample <- sample_data[sample(1:n, n, replace = TRUE), ]
    }
    
    
    ### STEP 2: We give the outcome model, use the bootstrap sample and save the 
    ### predictive outcome values from the model and standardize them.
    ### NOTE: We only use units from the bootstrap sample with observed Y.
    outcome_mod <- multinom(outcome_model, 
                            data = boot_sample,
                            maxit = 1000)
    
    # Takes the predicted values of C-1 categories
    pred_y <- fitted(outcome_mod)[, -1] 
    
    # Standardizes the predicted values
    mean_y_pred <- colMeans(pred_y)
    sd_y_pred <- apply(pred_y, 2, function(x) sd(x)) 
    for(i in 1:ncol(pred_y)){
      pred_y[,i] <- (pred_y[,i] - mean_y_pred[i]) / sd_y_pred[i]
    }
    
    
    ### STEP 3: We give the response model, save the predictive response value
    ### from the the model and standardize them
    ### NOTE: We use all units form the bootstrap sample 
    response_mod <- glm(response_model, 
                        data = boot_sample, 
                        family = binomial,
                        control = list(maxit = 100))
    pred_R <- fitted(response_mod) # Predictive response value
  
    # Standardizes the predictive values
    mean_R_pred <- mean(pred_R)
    sd_R_pred <- sd(pred_R)
    pred_R <- (pred_R - mean_R_pred) / sd_R_pred
    
    
    ### STEP 4: We use a Euclidean distance function to calculate the similarity 
    ### between unit i with missing Y in the original data set and unit j with 
    ### observed Y in the bootstrap sample.
    
    # Selects all cases in original data set with missing Y
    to_be_imputed <- sample_data[sample_data$response == 0, -1]
    
    # We calculate the predictive outcome and response values for these units
    pred_y_missing <- predict(outcome_mod, # Saves C-1 predictive outcome values
                              newdata = to_be_imputed, 
                              type = "probs")[ , -1] 
    pred_R_missing <- predict(response_mod, # Saves the response value
                              newdata = to_be_imputed, 
                              type = "response") 
    
    # Standardizes all the predictive values
    for(i in 1:ncol(pred_y_missing)){
      pred_y_missing[, i] <- (pred_y_missing [, i] 
                              - mean_y_pred[i]) / sd_y_pred[i]
    }
    pred_R_missing <- (pred_R_missing - mean_R_pred) / sd_R_pred
    
    
    # Creates a matrix with the predictive outcome and response values for units 
    # with missing Y in the original sample
    pred_to_be_imputed <- cbind(pred_y_missing, pred_R_missing)
    
    # Creates a matrix with the predictive outcome and response values for units 
    # with observed Y in the bootstrap sample
    pred_boot_obs <- cbind(pred_y, pred_R[boot_sample$response == 1])
    
    # Saves the indices of the bootstrap units with observed Y
    ind <- which(boot_sample$response == 1)
    
    # Defines the fixed weights for the Euclidean distance function
    w <- 1/ C
    
    # We calculate Euclidean distances. Returns a matrix in which each column 
    # contains the distances from a unit i with missing Y in the original data 
    # set to all units j with observed Y in the bootstrap sample.
    distances <- apply(pred_to_be_imputed, 1, function(x) {
      apply(pred_boot_obs, 1, fun_distance, x, weights = w)
    })
    
    
    ### STEP 5: For each unit with missing Y in the original data set, the
    ### imputation set is defined as k-Nearest Neighbor (in default k = 5)
    ### and the missing values of Y are imputed by randomly drawing a value of
    ### one donor from the k-Nearest Neighborhood. 
    
    # Orders the distances from smallest to largest
    distance_ranks <- apply(distances, 2, order)
    
    # Units are represent by each column. We get for each unit with  missing Y 
    # (in the original data set) the indices of five other units with observed 
    # Y (in the bootstrap sample) with the smallest calculated distance.
    indices_smallest_distance <- apply(distance_ranks, 2, 
                                       function(x) ind[x<=5])
    
    
    # Creates a list of possible donors for each unit with missing Y 
    donors <- apply(indices_smallest_distance, 2, function(x){
      boot_sample[x, "y"]
    })
    
    # Randomly draws one donor for each unit with missing Y in the original 
    # data set
    new_Y_values <- apply(donors, 2, function(x) sample(x, 1))
    
    # Imputes
    imputed_set <- sample_data
    imputed_set$y[imputed_set$response == 0] <- new_Y_values
    imputed_data[[m]] <- imputed_set
    
  }
  
  
  return(imputed_data) # Returns a list of the M imputed data sets
}



####################################################################
############################## MRPMMI ############################## 
####################################################################

# Functions: Performs MRPMMI with fixed weights M times.

# Parameters:
# 1. sample_data: data frame of the sample data: includes variables "y" 
#                 (categorical with more than 2 categories), auxiliary 
#                 variables and response indicator "response" 
#                 (binomial: y = is observed, 0 = y is missing)
# 2. outcome_models: a list of the outcome models written in a string format 
#                   (e.g. "y ~ x1 + x2 + x4")
# 4. M: number of times data is imputed. In default set to 5.

# Returns: a list of the imputed data sets.

MRPMMI <- function(sample_data, outcome_models, M = 5){
  
  source("Other_fun.R") # Script that the includes supportive functions
  n <- nrow(sample_data) # Number of units in the sample
  C <- length(levels(sample_data$y)) # Number of categories of Y
  
  imputed_data <- list() # List of the imputed data sets
  
  # We impute M = 5 times to create M data sets
  for(m in 1:M){
    
    
    ### STEP 1: We take a bootstrap of size n 
    boot_sample <- sample_data[sample(1:n, n, replace = TRUE), ]
    
    # We check that all the categories are represented in the bootstrap sample.
    # If not, we take a new sample.
    while(length(unique(boot_sample[boot_sample$response == 1, "y"])) != 3){
      boot_sample <- sample_data[sample(1:n, n, replace = TRUE), ]
    }
    
    
    ### STEP 2: We give the outcome models, save the predicted outcome values 
    ### from each model and use them in a combined regression model as 
    ### predictors
    
    y_pred_all <- list() # List of all predictions from all outcome models
    all_outcome_models <- list() # List of all outcome models
    
    # Fits the models, saves then and the predicted values for C-1 
    # categories of y
    for(i in 1:length(outcome_models)){
      outcome_mod <- multinom(outcome_models[[i]], 
                              data = boot_sample,
                              maxit = 1000) 
      all_outcome_models[[i]] <- outcome_mod
      y_pred_all[[i]] <-  fitted(outcome_mod)[, -1] # predicted values
    }
    
    # Combines all the predicted values into a data frame
    y_pred_all <- as.data.frame(do.call(cbind, y_pred_all))
    
    # Data frame containing units with observed y and their predicted values
    # NOTE: We only choose units in the bootstrap sample with observed Y.
    outcome_data <- as.data.frame(cbind(
      boot_sample[boot_sample$response == 1, "y"], y_pred_all))
    colnames(outcome_data) <- c("y", paste0("V", 1:ncol(y_pred_all)))
    
    
    # Next we regress the outcome variable on the saved predicted values and 
    # save the new predicted values.
    combined_outcome_model <- multinom(y ~ ., 
                                       data = outcome_data,
                                       maxit = 1000)
    y_pred <- fitted(combined_outcome_model)[, -1] # The final predicted values
    
    # Standardizes the final predicted values.
    mean_y_pred <- colMeans(y_pred) 
    sd_y_pred <- apply(y_pred, 2, function(x) sd(x))
    
    for(i in 1:ncol(y_pred)){
      y_pred[, i] <- (y_pred[, i] - mean_y_pred[i]) / sd_y_pred[i]
    }
    
    
    ### STEP 3: We calculate the predicted outcome probabilities for the 
    ### nonrespondents in the original sample and standardize these.
    
    # Selects all cases in original data set with missing Y
    to_be_imputed <- sample_data[sample_data$response == 0, -1]
    
    # First we calculate the predicted outcome values:
    
    # List of all predicted values for units with missing Y
    missing_pred_y_all <- list() 
    
    # Calculates the C - 1 predicted values by using the outcome models
    for(i in 1:length(all_outcome_models)){
      
      # If to_be_imputed has one unit
      if (nrow(to_be_imputed) == 1){
        missing_pred_y_all[[i]] <- predict(all_outcome_models[[i]], 
                                           newdata = to_be_imputed, 
                                           type = "probs")[-1]
        
      # If to_be_imputed has more than one unit  
      } else {
        missing_pred_y_all[[i]] <- predict(all_outcome_models[[i]], 
                                           newdata = to_be_imputed, 
                                           type = "probs")[, -1]
        
      }
    }
    
    
    # If to_be_imputed has one unit
    if(nrow(to_be_imputed) == 1){
      # Combines all the predicted values into a data frame
      missing_pred_y_all <- as.data.frame(t(unlist(missing_pred_y_all)))
      colnames(missing_pred_y_all)  <- c(paste0("V", 1:ncol(missing_pred_y_all)))
      
      # Calculates the final predicted values
      missing_pred_y <- data.frame(t(predict(combined_outcome_model, 
                                             newdata = missing_pred_y_all, 
                                             type = "probs")[-1]))
      
      # Standardizes the predicted values
      for(i in 1:ncol(missing_pred_y )){
        missing_pred_y[i] <- (missing_pred_y[i] - mean_y_pred[i]) / sd_y_pred[i]
      }
      
    
      # If to_be_imputed has more than one unit  
    } else {
      # Combines all the predicted values into a data frame
      missing_pred_y_all <- as.data.frame(do.call(cbind, missing_pred_y_all))
      colnames(missing_pred_y_all)  <- c(paste0("V", 1:ncol(missing_pred_y_all)))
      
      # Calculates the final predicted values
      missing_pred_y <- predict(combined_outcome_model, 
                                newdata = missing_pred_y_all, 
                                type = "probs")[, -1]
      
      # Standardizes the predicted values
      for(i in 1:ncol(missing_pred_y )){
        missing_pred_y[,i] <- (missing_pred_y[,i] - mean_y_pred[i]) / sd_y_pred[i]
      }
      
      
    }
    
    
    ### STEP 4: We use a Euclidean distance function to calculate the similarity 
    ### between unit i with missing Y in the original data set and unit j with 
    ### observed Y in the bootstrap sample.
    
    # Creates a matrix with the final predicted outcome and response values
    # for units with missing Y in the original sample
    pred_to_be_imputed <- cbind(missing_pred_y) 
    
    # Creates a matrix with the final predicted outcome and response values
    # for units with observed Y in the bootstrap sample
    pred_boot_obs <- cbind(y_pred)
    
    
    # Saves the indices of the bootstrap units with observed Y
    ind <- which(boot_sample$response == 1)
    
    # Defines the fixed weights for the Euclidean distance function
    w <- 1/ (C - 1)
    
    # We calculate Euclidean distances. Returns a matrix in which each column 
    # contains the distances from a unit i with missing Y in the original data 
    # set to all units j with observed Y in the bootstrap sample.
    distances <- apply(pred_to_be_imputed, 1, function(x) 
      apply(pred_boot_obs, 1, fun_distance, x, weights = w))
    
    
    ### STEP 5: For each unit with missing Y in the original data set, the
    ### imputation set is defined as k-Nearest Neighbor (in default k = 5)
    ### and the missing values of Y are imputed by randomly drawing a value of
    ### one donor from the k-Nearest Neighborhood. 
    
    # Ranks the distances from smallest to largest
    distance_rank <- apply(distances, 2, order)
    
    # Units are represent by each column. We get for each unit with  missing Y 
    # (in the original data set) the indices of five other units with 
    # observed Y (in the bootstrap sample) with the smallest calculated distance.
    indices_smallest_distance <- apply(distance_rank, 2, function(x) ind[x<=5])
    
    # Creates a list of possible donors for each unit with missing Y 
    donors <- apply(indices_smallest_distance, 2, function(x){
      boot_sample[x, "y"]
    })
    
    # Randomly draws one donor for each unit with missing Y in the original 
    # data set.
    new_Y_values <- apply(donors, 2, function(x) sample(x, 1))
    
    # Imputes
    imputed_set <- sample_data
    imputed_set$y[imputed_set$response == 0] <- new_Y_values
    imputed_data[[m]] <- imputed_set
    
  }
  
  return(imputed_data) # Returns a list of the M imputed data sets
}



####################################################################
################ MRIC WEIGHTED BASED ON PSEUDO-R2 ##################
####################################################################

# Functions: Performs MRIC in which weights are based on pseudo-R2
# Reference: MI steps based on MRNMMI by Breemer (2022): 
#            https://github.com/LigayaBreemer/MRNNMI 
#            File: PILOT_STUDY.R L:73-217
#            File: MRNMMI_fun.R: L: 4-166


# Parameters:
# 1. sample_data: data frame of the sample data: includes variables "y" 
#                 (categorical with more than 2 categories), auxiliary 
#                 variables and response indicator "response" 
#                 (binomial: y = is observed, 0 = y is missing)
# 2. outcome_models: a list of the outcome models written in a string format 
#                   (e.g. "y ~ x1 + x2 + x4")
# 3. response_models: a list of the response models written in a string 
#                     format (e.g. "y ~ x1 + x2 + x4")
# 4. pseudoR2: The type of pseudoR2 which can be either "McFadden", "CoxSnell",
#              "McKelveyZavoina" or "Nagelkerke"
# 4. M: number of times data is imputed. In default set to 5.

# Returns: a list of the imputed data sets and error counts.


MRIC_PseudoR2 <- function(sample_data, 
                          outcome_models, 
                          response_models, 
                          pseudoR2, M = 5){
  
  source("pseudoR2.R") # Script that the includes supportive functions
  source("Other_fun.R")
  n <- nrow(sample_data) # Number of units in the sample
  C <- length(levels(sample_data$y)) # number of categories of Y
  
  # Placeholders for error that the function did not converge (only for
  # Mckelvey and Zavoina)
  nonconv_count <- 0 
  NA_weight_count <- 0 # Placeholder for the error that the weights are 0
  
  imputed_data <- list() # List of the imputed data sets
  
  # We impute M = 5 times to create M data sets
  for(m in 1:M){
    
    ### STEP 1: We take a bootstrap of size n 
    boot_sample <- sample_data[sample(1:n, n, replace = TRUE), ]
    
    # We check that all the categories are represented in the bootstrap sample.
    # If not, we take a new sample.
    while(length(unique(boot_sample[boot_sample$response == 1, "y"])) != 3){
      boot_sample <- sample_data[sample(1:n, n, replace = TRUE), ]
    }
    
    ### STEP 2: We give the outcome models, save the predicted outcome values
    ### from each model and use them in a combined regression model as 
    ### predictors
    
    y_pred_all <- list() # List of all predictions from all outcome models
    all_outcome_models <- list() # List of all outcome models
    
    # Fits the models, saves then and the predicted values for C-1 
    # categories of y
    for(i in 1:length(outcome_models)){
      outcome_mod <- multinom(outcome_models[[i]], 
                              data = boot_sample,
                              maxit = 1000) 
      all_outcome_models[[i]] <- outcome_mod # saves the model
      y_pred_all[[i]] <-  fitted(outcome_mod)[, -1] # predictive values
    }
    
    # Combines all the predictive values into a data frame
    y_pred_all <- as.data.frame(do.call(cbind, y_pred_all))
    
    # Data frame containing units with observed y and their predicted values
    # NOTE: We only choose units in the bootstrap sample with observed Y
    outcome_data <- as.data.frame(cbind(
      boot_sample[boot_sample$response == 1, "y"],  y_pred_all))
    colnames(outcome_data) <- c("y", paste0("V", 1:ncol(y_pred_all)))
    
    
    # Next we regress the outcome variable on the saved predicted values and 
    # save the new predicted values
    combined_outcome_model <- multinom(y ~ .,
                                       data = outcome_data,
                                       maxit = 1000)
    y_pred <- fitted(combined_outcome_model)[, -1] # The final predicted values
    
    # Standardizes the final predicted values
    mean_y_pred <- colMeans(y_pred) 
    sd_y_pred <- apply(y_pred, 2, function(x) sd(x))
    
    for(i in 1:ncol(y_pred)){
      y_pred[,i] <- (y_pred[,i] - mean_y_pred[i]) / sd_y_pred[i]
    }
    
    ### STEP 3: We give the response models, save the predictive response values 
    ### from each model and use them in the combined regression model 
    ### as predictors.
    
    R_pred_all <- list() # List of all predictive values
    all_response_models <- list() # List of all response models
    
    # Fits each response model, saves the model and predictive response values
    for(i in 1:length(response_models)){
      response_mod <- glm(response_models[[i]], 
                          data = boot_sample, 
                          family = "binomial",
                          control = list(maxit = 100))
      
      all_response_models[[i]] <- response_mod # saves the model
      R_pred_all[[i]] <-  fitted(response_mod) # predictive values
    }
    
    
    # Combines all predictive values into a data frame
    R_pred_all <- as.data.frame(do.call(cbind, R_pred_all))
    
    # Data frame containing all response indicators and their response values.
    # NOTE: We choose all the units in the bootstrap sample regardless if Y 
    # is missing or not.
    response_data <- as.data.frame(cbind(boot_sample$response, R_pred_all))
    colnames(response_data)<- c("R", paste0("V", 1:ncol(R_pred_all)))
    
    # Next we regress the "response" variable on the saved predictive values
    # and save the new predictive values. 
    combined_response_model <- glm(R ~ ., 
                                   data = response_data , 
                                   family = "binomial",
                                   control = list(maxit = 100))
    
    R_pred <- fitted(combined_response_model) # The final predictive values
    
    
    # Standardizes the final predictive response values
    mean_R_pred <- mean(R_pred )
    sd_R_pred <- sd(R_pred) 
    R_pred <- (R_pred - mean_R_pred) / sd_R_pred
    
    
    ### STEP 4: We calculate the predicted outcome and response values for the 
    ### nonrespondents in the original sample and standardize these.
    
    # Selects all cases in original data set with missing Y.
    to_be_imputed <- sample_data[sample_data$response == 0, -1]
    # 
    # First we calculate the predicted outcome values:
    
    # List of all predicted values for units with missing Y
    missing_pred_y_all <- list()
    
    # Calculates the C - 1 predicted values by using the outcome models
    for(i in 1:length(all_outcome_models)){
      
      # If to_be_imputed has one unit
      if (nrow(to_be_imputed) == 1){
        missing_pred_y_all[[i]] <- predict(all_outcome_models[[i]], 
                                           newdata = to_be_imputed, 
                                           type = "probs")[-1]
        
      } else {
        
        # If to_be_imputed has more than one unit
        missing_pred_y_all[[i]] <- predict(all_outcome_models[[i]], 
                                           newdata = to_be_imputed, 
                                           type = "probs")[, -1]
      }
    }
    
    
    # If to_be_imputed has one unit
    if(nrow(to_be_imputed) == 1){
      # Combines all the predicted values into a data frame
      missing_pred_y_all <- as.data.frame(t(unlist(missing_pred_y_all)))
      colnames(missing_pred_y_all)  <- c(paste0("V", 1:ncol(missing_pred_y_all)))
      
      # Calculates the final predictive values
      missing_pred_y <- data.frame(t(predict(combined_outcome_model, 
                                             newdata = missing_pred_y_all, 
                                             type = "probs")[-1]))
      
      # Standardizes the predictive values
      for(i in 1:ncol(missing_pred_y )){
        missing_pred_y[i] <- (missing_pred_y[i] - mean_y_pred[i]) / sd_y_pred[i]
      }
      
      
    } else {
      # Combines all the predictive values into a data frame
      missing_pred_y_all <- as.data.frame(do.call(cbind, missing_pred_y_all))
      colnames(missing_pred_y_all)  <- c(paste0("V", 1:ncol(missing_pred_y_all)))
      
      # Calculates the final predictive values
      missing_pred_y <- predict(combined_outcome_model, 
                                newdata = missing_pred_y_all, 
                                type = "probs")[, -1]
      
      # Standardizes the predicted values
      for(i in 1:ncol(missing_pred_y )){
        missing_pred_y[,i] <- (missing_pred_y[,i] - mean_y_pred[i]) / sd_y_pred[i]
      }
      
      
    }
    
    # Next we calculate predictive response values:
    
    # List of all predictive response values for units with missing Y
    missing_pred_R_all <- list()  
    
    # Calculates the predictive values by using response models 
    for(i in 1:length(all_response_models)){
      missing_pred_R_all[[i]] <- predict(all_response_models[[i]], 
                                         newdata = to_be_imputed, 
                                         type = "response")
    }
    
    
    # Combines all predictive response values into a data frame
    missing_pred_R_all <- as.data.frame(do.call(cbind, missing_pred_R_all))
    
    # Calculates the final predictive values
    missing_pred_R <- predict(combined_response_model, 
                              newdata = missing_pred_R_all, 
                              type = "response")
    
    
    # Standardizes the final predictive response values
    missing_pred_R <- (missing_pred_R - mean_R_pred) / sd_R_pred
    
    
    ### STEP 5: We use an Euclidean distance function to calculate the 
    ### similarity between unit i with missing Y in the original data set and 
    ### unit j with observed Y in the bootstrap sample.
    
    # Creates a matrix with the final predictive outcome and response values
    # for units with missing Y in the original sample
    pred_to_be_imputed <- cbind(missing_pred_y, missing_pred_R) 
    
    # Creates a matrix with the final predictive outcome and response values
    # for units with observed Y in the bootstrap sample
    pred_boot_obs <- cbind(y_pred, R_pred[boot_sample$response == 1])
    
    
    # Saves the indices of the bootstrap units with observed Y.
    ind <- which(boot_sample$response == 1)
    
    # Calculates the pseudo-R2's that will be used to calculate the weights
    # for the distance function. The first C-1 values in the vector are 
    # pseudo-R2s for the C-1 final predictive outcome values and the last one is 
    # the pseudo-R2 for the final response value
    
    # If pseudo-R2 is Mckelvey and Zavoina, we also count if the method
    # was able ot converge when calculating the pseudo-R2
    if(pseudoR2 == "McKelveyZavoina"){
      R2_output <- R2_pseudo(pseudoR2, combined_outcome_model, outcome_data,
                             combined_response_model, response_data)
      R2_values <- R2_output[[1]]
      
      NA_weight_count <- NA_weight_count + R2_output[[2]]
      nonconv_count <- nonconv_count + R2_output[[3]]
      

    } else {
      R2_output <- R2_pseudo(pseudoR2, combined_outcome_model, outcome_data,
                             combined_response_model, response_data)
      R2_values <- R2_output[[1]]
      NA_weight_count <- NA_weight_count + R2_output[[2]]
    }
    
    
    # Calculates 1-  response_model_pseudoR2: the percentage that is left
    # to weight the predictive outcome values
    weights_for_pred_scores <-  1 - tail(R2_values, n = 1)
    
    # Rescales the predictive outcome pseudo-R2's so that they will add up to 1
    R2_outcome_rescale <-  R2_values[-length(R2_values)] / 
      sum(R2_values[-length(R2_values)])
    
    # Next we calculate the final weights for the predictive outcome values 
    # by rescaling R2 again so that they will add up to 
    # 1- response_model_pseudoR2
    R2_weight_outcome <- R2_outcome_rescale * weights_for_pred_scores
    
    # Combines the final weights in one vectors: first C-1 are the rescaled
    # (pseudo-R2) weights for the final predictive outcome values and the last 
    # one is (pseudo-R2) weight for the final predictive response value. 
    R2_final_weights <- c(R2_weight_outcome, tail(R2_values, n = 1))
    
    
    # We calculate Euclidean distances. Returns a matrix in which each column 
    # contains the distances  from a unit i with missing Y in the original data 
    # set to all units j with observed Y in the bootstrap sample.
    distances <- apply(pred_to_be_imputed, 1, function(x) 
      apply(pred_boot_obs, 1, fun_distance, x, w = R2_final_weights))
    
    
    ### STEP 6: For each unit with missing Y in the original data set, the
    ### imputation set is defined as k-Nearest Neighbor (in default k = 5)
    ### and the missing values of Y are imputed by randomly drawing a value of
    ### one donor from the k-Nearest Neighborhood. 
    
    # Ranks the distances from smallest to largest
    distance_rank <- apply(distances, 2, order)
    
    
    # Units are represent by each column. We get for each unit with  missing Y 
    # (in the original data set) the indices of five other units with 
    # observed Y (in the bootstrap sample) with the smallest calculated distance.
    indices_smallest_distance <- apply(distance_rank, 2, function(x) ind[x<=5])
    
    # Creates a list of possible donors for each unit with missing Y 
    donors <- apply(indices_smallest_distance, 2, function(x){
      boot_sample[x, "y"] 
    })
    
    # Randomly draws one donor for each unit with missing Y in the original 
    # data set.
    new_Y_values <- apply(donors, 2, function(x) sample(x, 1))
    
    # Imputes
    imputed_set <- sample_data
    imputed_set$y[imputed_set$response == 0] <- new_Y_values
    imputed_data[[m]] <- imputed_set
    
    
  }
  
  # Returns a list of the M imputed data sets and the error counts
  return(list(imputed_data, NA_weight_count, nonconv_count))
}



####################################################################
################### MRIC WEIGHTED BASED ON AIC ##################### 
####################################################################

# Functions: Performs MRIC in which weights are based on AIC/Akaike weights
# Reference: MI partly based on MRNMMI by Breemer (2022): 
#            https://github.com/LigayaBreemer/MRNNMI 
#            File: PILOT_STUDY.R L:73-217
#            File: MRNMMI_fun.R: L: 4-166


# Parameters:
# 1. sample_data: data frame of the sample data: includes variables "y" 
#                 (categorical with more than 2 categories), auxiliary 
#                 variables and response indicator "response" 
#                 (binomial: y = is observed, 0 = y is missing)
# 2. outcome_models: a list of the outcome models written in a string format 
#                   (e.g. "y ~ x1 + x2 + x4")
# 3. response_models: a list of the response models written in a string 
#                     format (e.g. "y ~ x1 + x2 + x4")
# 4. M: number of times data is imputed. In default set to 5.

# Returns: a list of the imputed data sets and error count.

MRIC_AIC <- function(sample_data, outcome_models, response_models, M = 5){
  
  source("Other_fun.R") # Script that includes the supportive functions
  n <- nrow(sample_data) # Number of units in the sample
  
  NA_weight_count <- 0 # error count
  
  imputed_data <- list() # List of the imputed data sets
  
  # We impute M = 5 times to create M data sets
  for(m in 1:M){
    
    
    ### STEP 1: We take a bootstrap of size n 
    boot_sample <- sample_data[sample(1:n, n, replace = TRUE), ]
    
    # We check that all the categories are represented in the bootstrap sample.
    # If not, we take a new sample.
    while(length(unique(boot_sample[boot_sample$response == 1, "y"])) != 3){
      boot_sample <- sample_data[sample(1:n, n, replace = TRUE), ]
    }
    
    
    ### STEP 2: We give the outcome models, save the predicted outcome values
    ### from each model and use them in a combined regression model as 
    ### predictors.
    
    y_pred_all <- list() # List of all predictions from all outcome models
    all_outcome_models <- list() # List of all outcome models
    AIC_scores_outcome <- list() # List of AIC scores of the outcome models
    
    # Fits the models, saves then and the predicted values for C-1 
    # categories of Y, the model and the AIC-score of the model
    for(i in 1:length(outcome_models)){
      outcome_mod <- multinom(outcome_models[[i]], 
                              data = boot_sample,
                              maxit = 1000)
      AIC_scores_outcome[[i]] <- outcome_mod$AIC
      all_outcome_models[[i]] <- outcome_mod
      y_pred_all[[i]] <-  fitted(outcome_mod)[, -1]
    }
    
    # Combines all the predicted values into a data frame
    y_pred_all <- as.data.frame(do.call(cbind, y_pred_all))
    
    # Standardizes the final predicted values
    mean_y_pred <- colMeans(y_pred_all) 
    sd_y_pred <- apply(y_pred_all, 2, function(x) sd(x))
    
    for(i in 1:ncol(y_pred_all)){
      y_pred_all[,i] <- (y_pred_all[,i] - mean_y_pred[i]) / sd_y_pred[i]
    }
    
    
    ### STEP 3: We give the response models, save the predictive response values
    ### from each model and use them in the combined regression model 
    ### as predictors.
    
    R_pred_all <- list() # List of all predictive response values
    all_response_models <- list() # List of all response models
    AIC_scores_response <- list() # List of AIC scores of response models 
    
    # Fits each response model, saves the model the predictive values and
    # AIC-score of each model
    for(i in 1:length(response_models)){
      response_mod <- glm(response_models[[i]], 
                          data = boot_sample, 
                          family = "binomial",
                          control = list(maxit = 100))
      AIC_scores_response[[i]] <- response_mod$aic
      
      all_response_models[[i]] <- response_mod
      R_pred_all[[i]] <-  fitted(response_mod)
    }
    
    # Combines all response values into a data frame
    R_pred_all <- as.data.frame(do.call(cbind, R_pred_all))
    R_pred_all <- cbind(boot_sample$response, R_pred_all)
    
    # Standardizes the final response values
    mean_R_pred <- colMeans(R_pred_all)
    sd_R_pred <- apply(R_pred_all, 2, function(x) sd(x))
    
    for(i in 1:ncol(R_pred_all)){
      R_pred_all[,i] <- (R_pred_all[,i] - mean_R_pred[i]) / sd_R_pred[i]
    }
    
    
    ### STEP 4: We calculate the predicted outcome and response values for the 
    ### nonrespondents in the original sample and standardize these.
    
    # Selects all cases in original data set with missing Y.
    to_be_imputed <- sample_data[sample_data$response == 0, -1]
    
    # First we calculate the predictive outcome values:
    
    # List of all predictive values for units with missing Y
    missing_pred_y_all <- list() 
    
    # Calculates the C - 1 predictive values by using the outcome models
    for(i in 1:length(all_outcome_models)){
      
      # If to_be_imputed includes one unit
      if (nrow(to_be_imputed) == 1){
        missing_pred_y_all[[i]] <- predict(all_outcome_models[[i]], 
                                           newdata = to_be_imputed, 
                                           type = "probs")[-1]
        
      # If to_be_imputed includes more than one unit  
      } else {
        missing_pred_y_all[[i]] <- predict(all_outcome_models[[i]], 
                                           newdata = to_be_imputed, 
                                           type = "probs")[, -1]
        
      }
    }
    
    # If to_be_imputed includes one unit
    if(nrow(to_be_imputed) == 1){
      
      # Combines all predictive values into a data frame
      missing_pred_y_all <- as.data.frame(t(unlist(missing_pred_y_all)))
      colnames(missing_pred_y_all)  <- c(paste0("V", 
                                                1:ncol(missing_pred_y_all)))
      
      # Standardizes the predictive values
      for(i in 1:ncol(missing_pred_y )){
        missing_pred_y_all[i] <- (missing_pred_y_all[i] - mean_y_pred[i]) / 
          sd_y_pred[i]
      }
      
    # If to_be_imputed includes more than one unit  
    } else {
      # Combines all the predictive values into a data frame
      missing_pred_y_all <- as.data.frame(do.call(cbind, missing_pred_y_all))
      colnames(missing_pred_y_all)  <- c(paste0("V", 1:ncol(missing_pred_y_all)))
      
      # Standardizes the predictive values
      for(i in 1:ncol(missing_pred_y_all)){
        missing_pred_y_all[,i] <- (missing_pred_y_all[,i] - mean_y_pred[i])  / 
          sd_y_pred[i]
      }
      
      
    }
    
    
    # Next we calculate predicted response values:
    
    # List of all predicted response values for units with missing Y
    missing_pred_R_all <- list()  
    
    # Calculates the repsonse values by using response models 
    for(i in 1:length(all_response_models)){
      missing_pred_R_all[[i]] <- predict(all_response_models[[i]], 
                                         newdata = to_be_imputed, 
                                         type = "response")
    }
    
    # Combines all predicted response values into a data frame
    missing_pred_R_all <- as.data.frame(do.call(cbind, missing_pred_R_all))
    
    
    # Standardizes the final response values
    for(i in 1:ncol(missing_pred_R_all)){
      missing_pred_R_all[,i] <- (missing_pred_R_all[,i] 
                                 - mean_R_pred[i]) / sd_R_pred[i]
    }
    
    
    ### STEP 5: We use a Euclidean distance function to calculate the similarity 
    ### between unit i with missing Y in the original data set and unit j with 
    ### observed Y in the bootstrap sample.
    
    # Creates a matrix with the final predicted outcome and response values
    # for units with missing Y in the original sample
    pred_to_be_imputed <- cbind(missing_pred_y_all, missing_pred_R_all) 
    
    # Creates a matrix with the final predicted outcome and response values
    # for units with observed Y in the bootstrap sample
    pred_boot_obs <- cbind(y_pred_all, 
                           R_pred_all[boot_sample$response == 1, -1])
    
    
    # Saves the indices of the bootstrap units with observed Y
    ind <- which(boot_sample$response == 1)
    
    A_weights_outcome_output <- Akaike_w(unlist(AIC_scores_outcome))
    A_weights_response_output <- Akaike_w(unlist(AIC_scores_response))
    
    NA_weight_count <- NA_weight_count + A_weights_outcome_output[[2]] +
      A_weights_response_output[[2]]
    
    # Calculates the weights that will be used in the distance function
    A_weights <- c(A_weights_outcome_output[[1]], 
                   A_weights_response_output[[1]]) 
  
    
    n_model = length(AIC_scores_outcome) # Number of outcome models
    
    # Copies Akaike weight for the outcome values C-1 times for each model 
    # (so that all C-1  categories per model have a weight) and saves these in 
    # "A_weights" together with response Akaike weights.
    A_weights <- c(rep(A_weights[1:n_model], 
                       each = length(levels(sample_data$y)) - 1), 
                   A_weights[(n_model + 1):length(A_weights)])
    
    
    
    # We calculate Euclidean distances. Returns a matrix in which each column 
    # contains the distances from a unit i with missing Y in the original data 
    # set to all units j with observed Y in the bootstrap sample.
    distances <- apply(pred_to_be_imputed, 1, function(x) apply(pred_boot_obs, 1, 
                                                                fun_distance, 
                                                                x, w = A_weights))
    
    ### STEP 6: For each unit with missing Y in the original data set, the
    ### imputation set is defined as k-Nearest Neighbor (in default k = 5)
    ### and the missing values of Y are imputed by randomly drawing a value of
    ### one donor from the k-Nearest Neighborhood. 
    
    # Orders the distances from smallest to largest
    distance_rank <- apply(distances, 2, order)
    
    # Units are represent by each column. We get for each unit with  missing Y 
    # (in the original data set) the indices of five other units with 
    # observed Y (in the bootstrap sample) with the smallest calculated distance.
    indices_smallest_distance <- apply(distance_rank, 2, function(x) ind[x<=5])
    
    # Creates a list of possible donors for each unit with missing Y 
    donors <- apply(indices_smallest_distance, 2, function(x){
      boot_sample[x, "y"]
    })
    
    # Randomly draws one donor for each unit with missing Y in the original 
    # data set.
    new_Y_values <- apply(donors, 2, function(x) sample(x, 1))
    
    # Imputes
    imputed_set <- sample_data
    imputed_set$y[imputed_set$response == 0] <- new_Y_values
    imputed_data[[m]] <- imputed_set
    
  }
  
  # Returns a list of the M imputed data sets and the error count
  return(list(imputed_data, NA_weight_count))
}



####################################################################
####### MRIC WEIGHTED BASED ON HOSMER-LEMESHOW CHI-SQUARE ##########
####################################################################

# Functions: Performs MRIC in which weights are based on Hosmer-Lemeshow (HL)
#            chi-square statistic.
# Reference: MI based on MRNMMI by Breemer (2022): 
#            https://github.com/LigayaBreemer/MRNNMI 
#            File: PILOT_STUDY.R L:73-217
#            File: MRNMMI_fun.R: L: 4-166


# Parameters:
# 1. sample_data: data frame of the sample data: includes variables "y" 
#                 (categorical with more than 2 categories), auxiliary 
#                 variables and response indicator "response" 
#                 (binomial: y = is observed, 0 = y is missing)
# 2. outcome_models: a list of the outcome models written in a string format 
#                   (e.g. "y ~ x1 + x2 + x4")
# 3. response_models: a list of the response models written in a string 
#                     format (e.g. "y ~ x1 + x2 + x4")
# 4. M: number of times data is imputed. In default set to 5.

# Returns: a list of the imputed data sets and error count.

MRIC_HL <- function(sample_data, outcome_models, response_models, M = 5){
  
  source("Other_fun.R") # Script that the includes supportive functions
  n <- nrow(sample_data) # Number of units in the sample
  C <- length(levels(sample_data$y)) # Number of categories of Y
  
  NA_weight_count <- 0 # error count
  imputed_data <- list() # List of the imputed data sets
  
  # We impute M = 5 times to create M data sets
  for(m in 1:M){
    
    
    ### STEP 1: We take a bootstrap of size n 
    boot_sample <- sample_data[sample(1:n, n, replace = TRUE), ]
    
    # We check that all the categories are represented in the bootstrap sample.
    # If not, we take a new sample.
    while(length(unique(boot_sample[boot_sample$response == 1, "y"])) != 3){
      boot_sample <- sample_data[sample(1:n, n, replace = TRUE), ]
    }
    
    ### STEP 2: We give the outcome models, save the predicted outcome values 
    ### from each model and use them in a combined regression model as 
    ### predictors
    
    y_pred_all <- list() # List of all predictions from all outcome models
    all_outcome_models <- list() # List of all outcome models
    
    # Fits the models, saves then and the predicted outcome values for C-1 
    # categories of Y
    for(i in 1:length(outcome_models)){
      outcome_mod <- multinom(outcome_models[[i]],
                              data = boot_sample,
                              maxit = 1000) 
      all_outcome_models[[i]] <- outcome_mod
      y_pred_all[[i]] <-  fitted(outcome_mod)[ , -1] # predicted values
    }
    
    # Combines all the predicted values into a data frame
    y_pred_all <- as.data.frame(do.call(cbind, y_pred_all))
    
    # Data frame containing units with observed y and their predicted values
    # NOTE: We only choose units in the bootstrap sample with observed Y.
    outcome_data <- as.data.frame(cbind(
      boot_sample[boot_sample$response == 1, "y"], y_pred_all))
    colnames(outcome_data) <- c("y", paste0("V", 1:ncol(y_pred_all)))
    
    # Next we regress the outcome variable on the saved predicted values and 
    # save the new predicted values.
    combined_outcome_model <- multinom(y ~ .,
                                       data = outcome_data,
                                       maxit = 1000)
    y_pred <- fitted(combined_outcome_model)[, -1] # The final predicted values
    
    
    # Standardizes the final predicted values
    mean_y_pred <- colMeans(y_pred) 
    sd_y_pred <- apply(y_pred, 2, function(x) sd(x))
    
    for(i in 1:ncol(y_pred)){
      y_pred[,i] <- (y_pred[,i] - mean_y_pred[i]) / sd_y_pred[i]
    }
    
    
    ### STEP 3: We give the response models, save the predicted response values 
    ### from each model and use them in the combined regression model as 
    ### predictors.
    
    R_pred_all <- list() # List of all predicted response values
    all_response_models <- list() # List of all response models
    
    # Fits each response model, saves the model and predicted values
    for(i in 1:length(response_models)){
      response_mod <- glm(response_models[[i]], 
                          data = boot_sample, 
                          family = "binomial",
                          control = list(maxit = 100))
      
      all_response_models[[i]] <- response_mod
      R_pred_all[[i]] <-  fitted(response_mod)
    }
    
    # Combines all response values into a data frame
    R_pred_all <- as.data.frame(do.call(cbind, R_pred_all))
    
    # Data frame containing all response indicators and their response values
    # NOTE: We choose all the units in the bootstrap sample regardless if Y 
    # is missing or not.
    response_data <- as.data.frame(cbind(boot_sample$response, R_pred_all))
    colnames(response_data)<- c("R", paste0("V", 1:ncol(R_pred_all)))
    
    # Next we regress the "response" variable on the saved response values 
    # and save the new predicted response values 
    combined_response_model <- glm(R ~ ., 
                                   data = response_data , 
                                   family = "binomial",
                                   control = list(maxit = 100))
    
    R_pred <- fitted(combined_response_model) # The final response values
    
    
    # Standardizes the final predicted response values
    mean_R_pred <- mean(R_pred)
    sd_R_pred <- sd(R_pred) 
    R_pred <- (R_pred - mean_R_pred) / sd_R_pred
    
    ### STEP 4: We calculate the predicted outcome and response values for the 
    ### nonrespondents in the original sample and standardize these.
    
    # Selects all cases in original data set with missing Y.
    to_be_imputed <- sample_data[sample_data$response == 0, -1]
    
    # First we calculate the predicted outcome values:
    
    # List of all predicted values for units with missing Y
    missing_pred_y_all <- list() 
    
    # Calculates the C - 1 predicted values by using the outcome models
    for(i in 1:length(all_outcome_models)){
      missing_pred_y_all[[i]] <- predict(all_outcome_models[[i]], 
                                         newdata = to_be_imputed, 
                                         type = "probs")[, -1]
    }
    
    # Combines all the predicted values into a data frame
    missing_pred_y_all <- as.data.frame(do.call(cbind, missing_pred_y_all))
    colnames(missing_pred_y_all)  <- c(paste0("V", 1:ncol(missing_pred_y_all)))
    
    # Calculates the final predicted values
    missing_pred_y <- predict(combined_outcome_model, 
                              newdata = missing_pred_y_all, 
                              type = "probs")[, -1]
    
    # Standardizes the predicted values
    for(i in 1:ncol(missing_pred_y)){
      missing_pred_y[,i] <- (missing_pred_y[,i] 
                             - mean_y_pred[i]) / sd_y_pred[i]
    }
    
    # Next we calculate predicted response values:
    
    # List of all predicted response values for units with missing Y
    missing_pred_R_all <- list()  
    
    # Calculates the response values by using response models 
    for(i in 1:length(all_response_models)){
      missing_pred_R_all[[i]] <- predict(all_response_models[[i]], 
                                         newdata = to_be_imputed, 
                                         type = "response")
    }
    
    # Combines all response values into a data frame
    missing_pred_R_all <- as.data.frame(do.call(cbind, missing_pred_R_all))
    
    # Calculates the final predicted response values
    missing_pred_R <- predict(combined_response_model, 
                              newdata = missing_pred_R_all, 
                              type = "response")
    
    
    # Standardizes the final response values
    missing_pred_R <- (missing_pred_R - mean_R_pred) / sd_R_pred
    
    ### STEP 5: We use a Euclidean distance function to calculate the similarity 
    ### between unit i with missing Y in the original data set and unit j with 
    ### observed Y in the bootstrap sample.
    
    # Creates a matrix with the final predicted outcome and response values
    # for units with missing Y in the original sample
    pred_to_be_imputed <- cbind(missing_pred_y, missing_pred_R) 
    
    # Creates a matrix with the final predicted outcome and response values
    # for units with observed Y in the bootstrap sample
    pred_boot_obs <- cbind(y_pred, R_pred[boot_sample$response == 1])
    
    # Saves the indices of the bootstrap units with observed Y.
    ind <- which(boot_sample$response == 1)
    
    # Calculates the HL chi-square statistics that will be used to calculate 
    # the weights for the distance function. The first C-1 values in the vector 
    # are HL statistics for the C-1 final predictive outcome values and the last one is 
    # the pseudo-R2 for the final response value.
    HL_stat_output <- HL_fun(C = C, 
                      outcome_data = outcome_data,
                      combined_response_model = combined_response_model)
    
    NA_weight_count <- NA_weight_count + HL_stat_output[[2]]
    
    # Takes the inverse of HL statistics
    HL_stat <- 1 / HL_stat_output[[1]]
    
    # Divides each HL stat by the sum to that they add up to 1
    HL_stat <- HL_stat / (sum(HL_stat))
    
    
    # We calculate Euclidean distances. Returns a matrix in which each column 
    # contains the distances from a unit i with missing Y in the original data 
    # set to all units j with observed Y in the bootstrap sample.
    distances <- apply(pred_to_be_imputed, 1, function(x) 
      apply(pred_boot_obs, 1, fun_distance, x, w = HL_stat))
    
    
    ### STEP 6: For each unit with missing Y in the original data set, the
    ### imputation set is defined as k-Nearest Neighbor (in default k = 5)
    ### and the missing values of Y are imputed by randomly drawing a value of
    ### one donor from the k-Nearest Neighborhood. 
    
    # Ranks the distances from smallest to largest
    distance_rank <- apply(distances, 2, order)
    
    # Units are represent by each column. We get for each unit with  missing Y 
    # (in the original data set) the indices of five other units with 
    # observed Y (in the bootstrap sample) with the smallest calculated distance.
    indices_smallest_distance <- apply(distance_rank, 2, function(x) ind[x<=5])
    
    # Creates a list of possible donors for each unit with missing Y 
    donors <- apply(indices_smallest_distance, 2, function(x){
      boot_sample[x, "y"]
    })
    
    # Randomly draws one donor for each unit with missing Y in the original 
    # data set.
    new_Y_values <- apply(donors, 2, function(x) sample(x, 1))
    
    # Imputes
    imputed_set <- sample_data
    imputed_set$y[imputed_set$response == 0] <- new_Y_values
    imputed_data[[m]] <- imputed_set
    
  }
  
  # Returns a list of the M imputed data sets and the error count
  return(list(imputed_data, NA_weight_count)) 
}



