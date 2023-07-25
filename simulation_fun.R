library(mice)
source("MI_methods.R")
source("Other_fun.R")


##############################################################
#################### SUPPORTIVE FUNCTIONS ####################
##############################################################

# Function that calculates pooled parameters for MRIC R2 (McFadden)
MRIC_pseudoR2_MF_simfun <- function(sample_data, 
                                    outcome_models, 
                                    response_models, 
                                    population_proportions, # true proportions
                                    it){
  
  # Calling the MRIC method
  MRIC_PseudoR2_MF_output <- MRIC_PseudoR2(sample_data, 
                                           outcome_models,
                                           response_models, 
                                           "McFadden") 
  
  MRIC_PseudoR2_MF_sets <- MRIC_PseudoR2_MF_output[[1]] # Imputed datasets
  
  # Calls function that calculates the pooled parameters 
  MCIR_PseudoR2_MF_pooled_param <- 
    fun_pooled_param(MRIC_PseudoR2_MF_sets, population_proportions)
  MCIR_PseudoR2_MF_pooled_param$var <- c(1, 2, 3) # Category of Y
  MCIR_PseudoR2_MF_pooled_param$it <- it # iteration number
  
  # Computes if any of the iteration resulted in 0 weights
  NA_weight_count <- sum(MRIC_PseudoR2_MF_output[[2]])
  
  # Returns list including the pooled parameters and the number of times
  # weights were 0 
  return(list(MCIR_PseudoR2_MF_pooled_param, NA_weight_count))
} 



# Function that calculates pooled parameters for MRIC R2 (Cox and Snell)
MRIC_pseudoR2_CS_simfun <-function(sample_data, 
                                   outcome_models, 
                                   response_models,
                                   population_proportions, 
                                   it){
  
  # Calling the MRIC method
  MRIC_PseudoR2_CS_output <- MRIC_PseudoR2(sample_data, 
                                           outcome_models, 
                                           response_models,
                                           "CoxSnell") 
  
  MRIC_PseudoR2_CS_sets <- MRIC_PseudoR2_CS_output[[1]] # Imputed datasets
  
  # Calls function that calculates the pooled parameters 
  MCIR_PseudoR2_CS_pooled_param <- 
    fun_pooled_param(MRIC_PseudoR2_CS_sets, population_proportions)
  
  MCIR_PseudoR2_CS_pooled_param$var <- c(1, 2, 3)
  MCIR_PseudoR2_CS_pooled_param$it <- it
  
  # Computes if any of the iteration resulted in 0 weights
  NA_weight_count <- sum(MRIC_PseudoR2_CS_output[[2]])
  
  # Returns list including the pooled parameters and the number of times
  # weights were 0 
  return(list(MCIR_PseudoR2_CS_pooled_param, NA_weight_count))
} 



# Function that calculates pooled parameters for MRIC R2 (McKelvey and Zavoina)
MRIC_pseudoR2_MZ_simfun <-function(sample_data, 
                                   outcome_models, 
                                   response_models,
                                   population_proportions, 
                                   it){
  
  # Calling the MRIC method
  MRIC_PseudoR2_MZ_output <- MRIC_PseudoR2(sample_data, 
                                           outcome_models, 
                                           response_models,
                                           "McKelveyZavoina")
  
  MRIC_PseudoR2_MZ_sets <- MRIC_PseudoR2_MZ_output[[1]] # Imputed datasets
  
  # Calls function that calculates the pooled parameters 
  MCIR_PseudoR2_MZ_pooled_param <- 
    fun_pooled_param(MRIC_PseudoR2_MZ_sets, population_proportions)
  
  MCIR_PseudoR2_MZ_pooled_param$var <- c(1, 2, 3)
  MCIR_PseudoR2_MZ_pooled_param$it <- it
  
  # Computes if any of the iteration resulted in 0 weights
  NA_weight_count <- sum(MRIC_PseudoR2_MZ_output[[2]]) 
  
  # Computes the number of times whether calculating R2 failed
  nonconv_count <- sum(MRIC_PseudoR2_MZ_output[[3]])
  
  # Returns list including the pooled parametersthem  number of times
  # weights were 0 and number of times the function was unable to compute R2
  return(list(MCIR_PseudoR2_MZ_pooled_param, NA_weight_count, nonconv_count))
} 



# Function that calculates pooled parameters for MRIC R2 (Nagelkerke)
MRIC_pseudoR2_N_simfun <-function(sample_data, 
                                  outcome_models, 
                                  response_models,
                                  population_proportions, 
                                  it){
  
  # Calling the MRIC method
  MRIC_PseudoR2_N_output <- MRIC_PseudoR2(sample_data, 
                                          outcome_models, 
                                          response_models,
                                          "Nagelkerke") 
  
  
  MRIC_PseudoR2_N_sets <- MRIC_PseudoR2_N_output[[1]] # Imputed datasets
  
  # Calls function that calculates the pooled parameters 
  MCIR_PseudoR2_N_pooled_param <- 
    fun_pooled_param(MRIC_PseudoR2_N_sets, population_proportions)
  
  MCIR_PseudoR2_N_pooled_param$var <- c(1, 2, 3)
  MCIR_PseudoR2_N_pooled_param$it <- it
  
  # Computes if any of the iteration resulted in 0 weights
  NA_weight_count <- sum(MRIC_PseudoR2_N_output[[2]])
  
  
  # Returns list including the pooled parameters and the number of times
  # weights were 0 
  return(list(MCIR_PseudoR2_N_pooled_param, NA_weight_count))
}



# Function that calculates pooled parameters for MRIC with Akaike weights
MRIC_AIC_simfun <-function(sample_data, 
                           outcome_models, 
                           response_models,
                           population_proportions, 
                           it){
  
  # Calling the MRIC method
  MRIC_AIC_output <- MRIC_AIC(sample_data, 
                              outcome_models, 
                              response_models)
  
  MRIC_AIC_sets <-  MRIC_AIC_output[[1]] # Imputed datasets
  
  # Calls function that calculates the pooled parameters 
  MCIR_AIC_pooled_param <- 
    fun_pooled_param(MRIC_AIC_sets, population_proportions)
  
  MCIR_AIC_pooled_param$var <- c(1, 2, 3)
  MCIR_AIC_pooled_param$it <- it
  
  # Computes if any of the iteration resulted in 0 weights 
  NA_weight_count <- sum(MRIC_AIC_output[[2]])
  
  # Returns list including the pooled parameters and the number of times
  # weights were 0 
  return(list(MCIR_AIC_pooled_param, NA_weight_count))
} 


# Function that calculates pooled parameters for MRIC Hosmer-Lemeshow test 
# statistic
MRIC_HL_simfun <-function(MRIC_HL_sets, # Imputed datasets 
                          population_proportions,
                          it){
  
  # Calls function that calculates the pooled parameters 
  MCIR_HL_pooled_param <- 
    fun_pooled_param(MRIC_HL_sets[[1]], population_proportions)
  
  MCIR_HL_pooled_param$var <- c(1, 2, 3)
  MCIR_HL_pooled_param$it <- it
  
  # Computes if any of the iteration resulted in 0 weights 
  NA_weight_count <- sum(MRIC_HL_sets[[2]])
  
  # Returns list including the pooled parameters and the number of times
  # weights were 0 
  return(list(MCIR_HL_pooled_param, NA_weight_count))
} 


# Function that calculated pooled parameters for MRNNMI
MRNNMI_simfun <-function(sample_data, 
                         outcome_models, 
                         response_models,
                         population_proportions,
                         it){
  
  # Imputed datasets
  MRNNMI_sets <- MRNNMI(sample_data, outcome_models, response_models) 
  
  # Pooled parameters
  MRNNMI_pooled_param <- fun_pooled_param(MRNNMI_sets, 
                                          population_proportions)
  MRNNMI_pooled_param$var <- c(1, 2, 3)
  MRNNMI_pooled_param$it <- it
  
  # Returns pooled parameters
  return(MRNNMI_pooled_param)
} 



# Function that calculated pooled parameters for DRNNMI
DRNNMI_simfun <-function(sample_data, 
                         outcome_models, 
                         response_models,
                         population_proportions, 
                         it){
  # Imputed datasets
  DRNNMI_sets <- DRNNMI(sample_data, outcome_models[[1]],  response_models[[1]])
  
  # Pooled parameters
  DRNNMI_pooled_param <- fun_pooled_param(DRNNMI_sets, population_proportions)
  DRNNMI_pooled_param$var <- c(1, 2, 3)
  DRNNMI_pooled_param$it <- it
  
  # Returns pooled parameters
  return(DRNNMI_pooled_param)
} 



# Function that calculated pooled parameters for MRPMMI
MRPMMI_simfun <-function(sample_data, 
                         outcome_models,
                         population_proportions, 
                         it){
  
  MRPMMI_sets <- MRPMMI(sample_data, outcome_models) # Imputed datasets
  
  # Pooled parameters
  MRPMMI_pooled_param <- fun_pooled_param(MRPMMI_sets, 
                                          population_proportions)
  MRPMMI_pooled_param$var <- c(1, 2, 3)
  MRPMMI_pooled_param$it <- it
  
  # Returns pooled parameters
  return(MRPMMI_pooled_param)
} 



# Function that calculates pooled parameters for MICE with default setting
MICE_default_simfun <- function(sample_data,
                                population_proportions,
                                it){
  
  MICE_sets <- mice(data = sample_data, m = 5, maxit = 20) # MICE
  
  # Extracts imputed datasets
  MICE_complete <- complete(MICE_sets, action = "long")
  MICE_sets <- list() 
  
  # Creating a list including all the imputed datasets
  for (i in 1:length(unique(MICE_complete$.imp))){
    MICE_sets[[i]] <- MICE_complete[MICE_complete$.imp == i, ]
  }
  
  # Pooled parameters
  MICE_pooled_param <- fun_pooled_param(MICE_sets,
                                        population_proportions)
  MICE_pooled_param$var <- c(1, 2, 3)
  MICE_pooled_param$it <- it
  
  # Returns pooled parameters
  return(MICE_pooled_param)
} 



# Function that calculates pooled parameters for MICE with outcome model
MICE_model_simfun <- function(sample_data,
                              population_proportions,
                              outcome_model, # outcome model used by MICE
                              it){
  
  MICE_sets <- mice(data = sample_data, # MICE 
                    m = 5, 
                    maxit = 20,
                    formula = outcome_model)
  
  # Extracts imputed datasets
  MICE_complete <- complete(MICE_sets, action = "long")
  MICE_sets <- list()
  
  # Creating a list including all the imputed datasets
  for (i in 1:length(unique(MICE_complete$.imp))){
    MICE_sets[[i]] <- MICE_complete[MICE_complete$.imp == i, ]
  }
  
  # Pooled parameters
  MICE_pooled_param <- fun_pooled_param(MICE_sets,
                                        population_proportions)
  MICE_pooled_param$var <- c(1, 2, 3)
  MICE_pooled_param$it <- it
  
  # Returns pooled parameters
  return(MICE_pooled_param)
} 



# Function that checks that the Hosmer-Lemeshow test can be done with the
# current sample
fun_HLError_check <- function(sample_data, outcome_models, response_models) {
  tryCatch( expr = {
    
    # Callms MRIC with Akaike weights
    MRIC_HL_output <- MRIC_HL(sample_data,  outcome_models,  response_models)
    
    # If there is no error, returns the imputed datasets
    return(MRIC_HL_output)
  },
  
  error = function(e) { 
    
    # If wrror occurs, returns an empty list
    MRIC_HL_output <- list()
    return(MRIC_HL_output)
  })
}



##############################################################
################## MAIN SIMULATION FUNCTIONS #################
##############################################################


# Simulation function for robust methods
# Parameters:
# 1. x:list of parameters that include the iteration number (i), population,
# sample size (n), outcome models, response models and wanted response rate. 
fun_parallel_robust <- function(x){
  source("simulation_fun.R")
  
  # Renaming parameters
  i               <- x[[1]]
  population      <- x[[2]]
  n               <- x[[3]]
  outcome_models  <- x[[4]]
  response_models <- x[[5]]
  response_rate   <- x[[6]]
  
  # Placeholders for the iterations errors and ARR
  it_stat_info <- c(0, 0, 0, 0, 0, 0, 0, 0, 0)
  names(it_stat_info) <- c("HL errors",
                           "MRIC MF NA weights",
                           "MRIC CS NA weights",
                           "MRIC MZ NA weights",
                           "MRIC MZ could not converge", 
                           "MRIC N NA weights",
                           "MRIC AIC NA weights",
                           "MRIC HL NA weights",
                           "Average response rate")
  
  N <- nrow(population) # Number of observations in the population
  possible_HL_error <- TRUE # HL error is true
  
  # Calculates the true population proportions
  population_proportions <- table(population$y)/nrow(population)
  
  set.seed(i) # Setting the seed
  
  # While HL error is true we draw a sample of size n and and compute the
  # the response rate
  while(possible_HL_error == TRUE){ 
    sample_data <-  population[sample(1:N, n, replace = FALSE), ]
    
    if(response_rate == 60){
      response_prob <- (exp(-19.3 
                            + 1.5 * sample_data$x1 
                            + 2.5 * sample_data$x2 
                            + 1.2 * sample_data$x3) /
                          (1 + exp(-19.3 
                                   + 1.5 * sample_data$x1 
                                   + 2.5 * sample_data$x2 
                                   + 1.2 * sample_data$x3)))
      
    } else if(response_rate == 70){
      
      response_prob <- (exp(-18.2 
                            + 1.5 * sample_data$x1 
                            + 2.5 * sample_data$x2 
                            + 1.2 * sample_data$x3) /
                          (1 + exp(-18.2 
                                   + 1.5 * sample_data$x1 
                                   + 2.5 * sample_data$x2 
                                   + 1.2 * sample_data$x3)))
      
    } else if(response_rate == 80) {
      
      response_prob <- (exp(-17
                            + 1.5 * sample_data$x1 
                            + 2.5 * sample_data$x2 
                            + 1.2 * sample_data$x3) /
                          (1 + exp(-17 
                                   + 1.5 * sample_data$x1 
                                   + 2.5 * sample_data$x2 
                                   + 1.2 * sample_data$x3)))
    } else{
      print("Correct response probability was not given. 
              The response probabiltiy can be either 60, 70 or 80.")
      break
    }
    
    # Calculates response indicator based on response probability
    sample_data$response <- rbinom(n, 1, response_prob)
    
    # Removes some Y values from the sample based on the response indicator
    sample_data$y[sample_data$response == 0] <- NA
    
    # We draw a new sample if the sample does not contain observations from all 
    # y categories (= {1, 2, 3, NA})
    while(length(unique(sample_data$y))!=4){ 
      sample_data <- population[sample(1:N, n, replace = FALSE), ]
      
      if(response_rate == 60){
        response_prob <- (exp(-19.3 
                              + 1.5 * sample_data$x1 
                              + 2.5 * sample_data$x2 
                              + 1.2 * sample_data$x3) /
                            (1 + exp(-19.3 
                                     + 1.5 * sample_data$x1 
                                     + 2.5 * sample_data$x2 
                                     + 1.2 * sample_data$x3)))
        
      } else if (response_rate == 70){
        
        response_prob <- (exp(-18.2 
                              + 1.5 * sample_data$x1 
                              + 2.5 * sample_data$x2 
                              + 1.2 * sample_data$x3) /
                            (1 + exp(-18.2 
                                     + 1.5 * sample_data$x1 
                                     + 2.5 * sample_data$x2 
                                     + 1.2 * sample_data$x3)))
        
      } else if (response_rate == 80) {
        
        response_prob <- (exp(-17
                              + 1.5 * sample_data$x1 
                              + 2.5 * sample_data$x2 
                              + 1.2 * sample_data$x3) /
                            (1 + exp(-17 
                                     + 1.5 * sample_data$x1 
                                     + 2.5 * sample_data$x2 
                                     + 1.2 * sample_data$x3)))
      } else{
        print("Correct response probability was not given. 
              The response probabiltiy can be either 60, 70 or 80.")
        break
      }
      
      
      
      sample_data$response <- rbinom(n, 1, response_prob)
      
      sample_data$y[sample_data$response == 0] <- NA
    }
    
    # Checks what the HL error functions returns
    MRIC_HL_sets_test <- fun_HLError_check(sample_data, 
                                           outcome_models, 
                                           response_models) 
    
    # If the function returns empty list, error remains TRUE
    if(length(MRIC_HL_sets_test) == 0) { 
      possible_HL_error <- TRUE
      
      # Computes how many times HL error occurs
      it_stat_info[1] <- it_stat_info[1] + 1
      
    } else {
      
      # If the function returns the imputed datasets, error changes to FALSE
      possible_HL_error <- FALSE
      
    }
  }
  
  
  # Calls the function that computes the pooled parameter estimates for MRIC
  # R2 (McFadden)
  MRIC_PseudoR2_MF_vec <- MRIC_pseudoR2_MF_simfun(sample_data,
                                                  outcome_models,
                                                  response_models,
                                                  population_proportions,
                                                  i)
  MRIC_PseudoR2_MF_result <- MRIC_PseudoR2_MF_vec[[1]] # pooled parameters
  it_stat_info[2] <- it_stat_info[2] + MRIC_PseudoR2_MF_vec[[2]] # error count
  
  # Calls the function that computes the pooled parameter estimates for MRIC
  # R2 (Cox and Snell)
  MRIC_PseudoR2_CS_vec <- MRIC_pseudoR2_CS_simfun(sample_data,
                                                  outcome_models,
                                                  response_models,
                                                  population_proportions,
                                                  i)
  MRIC_PseudoR2_CS_result <- MRIC_PseudoR2_CS_vec[[1]] # pooled parameters
  it_stat_info[3] <- it_stat_info[3] + MRIC_PseudoR2_CS_vec[[2]] # error count
  
  # Calls the function that computes the pooled parameter estimates for MRIC
  # R2 (McKelvey and Zavoina)
  MRIC_PseudoR2_MZ_vec <- MRIC_pseudoR2_MZ_simfun(sample_data,
                                                  outcome_models,
                                                  response_models,
                                                  population_proportions,
                                                  i)
  
  MRIC_PseudoR2_MZ_result <- MRIC_PseudoR2_MZ_vec[[1]] # pooled parameters
  
  it_stat_info[4] <- it_stat_info[4] + MRIC_PseudoR2_MZ_vec[[2]] # error counts
  it_stat_info[5] <- it_stat_info[5] + MRIC_PseudoR2_MZ_vec[[3]]
  
  # Calls the function that computes the pooled parameter estimates for MRIC
  # R2 (Nagelkerke)
  MRIC_PseudoR2_N_vec <- MRIC_pseudoR2_N_simfun(sample_data,
                                                outcome_models,
                                                response_models,
                                                population_proportions,
                                                i)
  MRIC_PseudoR2_N_result <- MRIC_PseudoR2_N_vec[[1]] # pooled parameters
  it_stat_info[6] <- it_stat_info[6] + MRIC_PseudoR2_N_vec[[2]] # error count
  
  
  # Calls the function that computes the pooled parameter estimates for MRIC
  # with Akaike weights
  MRIC_AIC_vec <- MRIC_AIC_simfun(sample_data,
                                  outcome_models,
                                  response_models,
                                  population_proportions,
                                  i)
  MRIC_AIC_result <- MRIC_AIC_vec[[1]] # pooled parameters
  it_stat_info[7] <- it_stat_info[7] + MRIC_AIC_vec[[2]] # error count
  
  # Calls the function that computes the pooled parameter estimates for MRIC
  # with HL statistic 
  MRIC_HL_vec <- MRIC_HL_simfun(MRIC_HL_sets_test,
                                population_proportions,
                                i)
  MRIC_HL_result <- MRIC_HL_vec[[1]] # pooled parameters
  it_stat_info[8] <- it_stat_info[8] + MRIC_HL_vec[[2]] # error count
  
  # Calls the function that computes the pooled parameter estimates for MRNNMI
  MRNNMI_result <- MRNNMI_simfun(sample_data,
                                 outcome_models,
                                 response_models,
                                 population_proportions,
                                 i)
  
  # Calls the function that computes the pooled parameter estimates for DRNNMI
  DRNNMI_result <- DRNNMI_simfun(sample_data,
                                 outcome_models,
                                 response_models,
                                 population_proportions,
                                 i)
  
  # Calls the function that computes the pooled parameter estimates for MRPMMI
  MRPMMI_result <- MRPMMI_simfun(sample_data,
                                 outcome_models,
                                 population_proportions,
                                 i)
  
  # Computed the average response rate in the sample
  it_stat_info[9] <- sum(sample_data$response, na.rm = T) / n
  
  # Returns the pooled estimates and the error counts
  all_pooled_results <- list(
    MRIC_PseudoR2_MF_result,
    MRIC_PseudoR2_CS_result,
    MRIC_PseudoR2_MZ_result,
    MRIC_PseudoR2_N_result,
    MRIC_AIC_result,
    MRIC_HL_result,
    MRNNMI_result,
    DRNNMI_result,
    MRPMMI_result,
    it_stat_info)
  
  # Renaming parameters
  i               <- x[[1]]
  population      <- x[[2]]
  n               <- x[[3]]
  outcome_models  <- x[[4]]
  response_models <- x[[5]]
  response_rate   <- x[[6]]
  file_name <- paste0("C:\\...\\robust\\n", n, "\\parallel_robust_results_i_", i, "_n_", n, "_RR_", response_rate, ".rds")
  
  # Saves all the results
  saveRDS(all_pooled_results, file = file_name)
  
}



# Simulation function for MICE variants
# Parameters:
# 1. x:list of parameters that include the iteration number (i), population,
# sample size (n), outcome models, response models and wanted response rate. 
fun_parallel_mice <- function(x){
  source("simulation_fun.R")
  
  # Renaming parameters
  i               <- x[[1]]
  population      <- x[[2]]
  n               <- x[[3]]
  outcome_models  <- x[[4]]
  response_models <- x[[5]]
  response_rate   <- x[[6]]
  
  it_stat_info <- c(0, 0)   # Placeholders for the iterations error and ARR
  names(it_stat_info) <- c("HL errors", "Average response rate")
  
  N <- nrow(population)  # Number of observations in the population
  possible_HL_error <- TRUE # HL error is true
  
  # Calculates the true population proportions
  population_proportions <- table(population$y)/nrow(population)
  
  set.seed(i) # Setting the seed
  
  # While HL error is true we draw a sample of size n and and compute the
  # the response rate
  while(possible_HL_error == TRUE){ 
    sample_data <-  population[sample(1:N, n, replace = FALSE), ]
    
    if(response_rate == 60){
      response_prob <- (exp(-19.3 
                            + 1.5 * sample_data$x1 
                            + 2.5 * sample_data$x2 
                            + 1.2 * sample_data$x3) /
                          (1 + exp(-19.3 
                                   + 1.5 * sample_data$x1 
                                   + 2.5 * sample_data$x2 
                                   + 1.2 * sample_data$x3)))
      
    } else if (response_rate == 70){
      
      response_prob <- (exp(-18.2 
                            + 1.5 * sample_data$x1 
                            + 2.5 * sample_data$x2 
                            + 1.2 * sample_data$x3) /
                          (1 + exp(-18.2 
                                   + 1.5 * sample_data$x1 
                                   + 2.5 * sample_data$x2 
                                   + 1.2 * sample_data$x3)))
      
    } else if (response_rate == 80) {
      
      response_prob <- (exp(-17
                            + 1.5 * sample_data$x1 
                            + 2.5 * sample_data$x2 
                            + 1.2 * sample_data$x3) /
                          (1 + exp(-17 
                                   + 1.5 * sample_data$x1 
                                   + 2.5 * sample_data$x2 
                                   + 1.2 * sample_data$x3)))
    } else{
      print("Correct response probability was not given. 
              The response probabiltiy can be either 60, 70 or 80.")
      break
    }
    
    # Calculates response indicator based on response probability
    sample_data$response <- rbinom(n, 1, response_prob)
    
    # Removes some Y values from the sample based on the response indicator
    sample_data$y[sample_data$response == 0] <- NA
    
    # We draw a new sample if the sample does not contain observations from all 
    # y categories (= {1, 2, 3, NA})
    while(length(unique(sample_data$y))!=4){ 
      sample_data <- population[sample(1:N, n, replace = FALSE), ]
      
      if(response_rate == 60){
        response_prob <- (exp(-19.3 
                              + 1.5 * sample_data$x1 
                              + 2.5 * sample_data$x2 
                              + 1.2 * sample_data$x3) /
                            (1 + exp(-19.3 
                                     + 1.5 * sample_data$x1 
                                     + 2.5 * sample_data$x2 
                                     + 1.2 * sample_data$x3)))
        
      } else if (response_rate == 70){
        
        response_prob <- (exp(-18.2 
                              + 1.5 * sample_data$x1 
                              + 2.5 * sample_data$x2 
                              + 1.2 * sample_data$x3) /
                            (1 + exp(-18.2 
                                     + 1.5 * sample_data$x1 
                                     + 2.5 * sample_data$x2 
                                     + 1.2 * sample_data$x3)))
        
      } else if (response_rate == 80) {
        
        response_prob <- (exp(-17
                              + 1.5 * sample_data$x1 
                              + 2.5 * sample_data$x2 
                              + 1.2 * sample_data$x3) /
                            (1 + exp(-17 
                                     + 1.5 * sample_data$x1 
                                     + 2.5 * sample_data$x2 
                                     + 1.2 * sample_data$x3)))
      } else{
        print("Correct response probability was not given. 
              The response probabiltiy can be either 60, 70 or 80.")
        break
      }
      
      
      sample_data$response <- rbinom(n, 1, response_prob)
      
      # We change some values to missing based on response 
      sample_data$y[sample_data$response == 0] <- NA
    }
    
    
    # Checks what the HL error functions returns
    MRIC_HL_sets_test <- fun_HLError_check(sample_data, 
                                           outcome_models, 
                                           response_models) 
    
    # If the function returns empty list, error remains TRUE
    if(length(MRIC_HL_sets_test) == 0) { 
      possible_HL_error <- TRUE
      
      # Computes how many times HL error occurs
      it_stat_info[1] <- it_stat_info[1] + 1
      
    } else {
      # If the function returns the imputed datasets, error changes to FALSE
      possible_HL_error <- FALSE
      
    }
  }
  
  # Outcome models used in different MICE variants
  MICE_models <- list(as.formula("y ~ x1 + x2 + x3"), 
                      as.formula("y ~ x1 + x3"), 
                      as.formula("y ~ x1 + x2 + x4 + x3:x4"))
  
  
  # Calculating the average response rate of the sample
  it_stat_info[2] <- sum(sample_data$response, na.rm = T) / n
  
  # Calls the function that computes the pooled parameter estimates for MICE
  # with the default settings
  MICE_default_result <- MICE_default_simfun(sample_data, 
                                             population_proportions, 
                                             i)
  
  # Calls the function that computes the pooled parameter estimates for MICE
  # with the correct outcome model
  MICE_CM_result <- MICE_model_simfun(sample_data, 
                                      population_proportions,
                                      MICE_models[1],
                                      i)
  
  # Calls the function that computes the pooled parameter estimates for MICE
  # with the incorrect nested outcome model
  MICE_INM_result <-  MICE_model_simfun(sample_data, 
                                        population_proportions,
                                        MICE_models[2],
                                        i)
  
  # Calls the function that computes the pooled parameter estimates for MICE
  # with the incorrect nnonnested outcome model
  MICE_INNM_result <-  MICE_model_simfun(sample_data, 
                                         population_proportions,
                                         MICE_models[3],
                                         i)
  
  # Returns the pooled parameter estimates
  all_pooled_results <- list(
    MICE_default_result,
    MICE_CM_result,
    MICE_INNM_result,
    MICE_INM_result,
    it_stat_info)
  
  file_name <- paste0("C:\\...\\mice\\n", n, "\\RR", response_rate, "\\parallel_mice_results__i_", i, "_n_", n, "_RR_", response_rate, ".rds")
  
  # Saves the results
  saveRDS(all_pooled_results, file = file_name)
  
  
}



# Simulation function to compute the complete-cases analysis
# Parameters:
# 1. x:list of parameters that include the iteration number (i), population,
# sample size (n), outcome models, response models and wanted response rate. 

fun_sample_proportions <- function(x){
  source("simulation_fun.R")
  # Renaming parameters
  i               <- x[[1]]
  population      <- x[[2]]
  n               <- x[[3]]
  outcome_models  <- x[[4]]
  response_models <- x[[5]]
  response_rate   <- x[[6]]
  
  
  N <- nrow(population)
  possible_HL_error <- TRUE
  
  # Calculates the true population proportions
  population_proportions <- table(population$y)/nrow(population)
  
  set.seed(i) # Setting the seed
  
  # Calculates the true population proportions
  population_proportions <- table(population$y)/nrow(population)
  
  
  # We draw a new sample if the sample does not contain observations from all 
  # y categories (= {1, 2, 3, NA})
  while(possible_HL_error == TRUE){ 
    sample_data <-  population[sample(1:N, n, replace = FALSE), ]
    
    if(response_rate == 60){
      response_prob <- (exp(-19.3 + 1.5 * sample_data$x1 + 2.5* sample_data$x2 + 1.2 * sample_data$x3) /
                          (1 + exp(-19.3 + 1.5 * sample_data$x1 + 2.5* sample_data$x2 
                                   + 1.2 * sample_data$x3)))
      
    } else if (response_rate == 70){
      
      response_prob <- (exp(-18.2 + 1.5 * sample_data$x1 + 2.5* sample_data$x2 + 1.2 * sample_data$x3) /
                          (1 + exp(-18.2 + 1.5 * sample_data$x1 + 2.5* sample_data$x2 
                                   + 1.2 * sample_data$x3)))
      
    } else if (response_rate == 80) {
      
      response_prob <- (exp(-17+ 1.5 * sample_data$x1 + 2.5* sample_data$x2 + 1.2 * sample_data$x3) /
                          (1 + exp(-17 + 1.5 * sample_data$x1 + 2.5* sample_data$x2 
                                   + 1.2 * sample_data$x3)))
    } else{
      print("Correct response probability was not given. 
              The response probabiltiy can be either 60, 70 or 80.")
      break
    }
    
    sample_data$response <- rbinom(n, 1, response_prob)
    
    # We change some values to missing based on response 
    sample_data$y[sample_data$response == 0] <- NA
    
    # We draw a new sample if the sample does not contain observations from all 
    # y categories (= {1, 2, 3, NA})
    while(length(unique(sample_data$y))!=4){ 
      sample_data <- population[sample(1:N, n, replace = FALSE), ]
      
      if(response_rate == 60){
        response_prob <- (exp(-19.3 + 1.5 * sample_data$x1 + 2.5* sample_data$x2 + 1.2 * sample_data$x3) /
                            (1 + exp(-19.3 + 1.5 * sample_data$x1 + 2.5* sample_data$x2 
                                     + 1.2 * sample_data$x3)))
        
      } else if (response_rate == 70){
        
        response_prob <- (exp(-18.2 + 1.5 * sample_data$x1 + 2.5* sample_data$x2 + 1.2 * sample_data$x3) /
                            (1 + exp(-18.2 + 1.5 * sample_data$x1 + 2.5* sample_data$x2 
                                     + 1.2 * sample_data$x3)))
        
      } else if (response_rate == 80) {
        
        response_prob <- (exp(-17+ 1.5 * sample_data$x1 + 2.5* sample_data$x2 + 1.2 * sample_data$x3) /
                            (1 + exp(-17 + 1.5 * sample_data$x1 + 2.5* sample_data$x2 
                                     + 1.2 * sample_data$x3)))
      } else{
        print("Correct response probability was not given. 
              The response probabiltiy can be either 60, 70 or 80.")
        break
      }
      
      
      
      sample_data$response <- rbinom(n, 1, response_prob)
      
      # We change some values to missing based on response 
      sample_data$y[sample_data$response == 0] <- NA
    }
    
    
    
    MRIC_HL_sets_test <- fun_HLError_check(sample_data, 
                                           outcome_models, 
                                           response_models) 
    
    
    
    if(length(MRIC_HL_sets_test) == 0) { 
      possible_HL_error <- TRUE
      
    } else {
      possible_HL_error <- FALSE
      
    }
    
  }
  
  sample_data_n <- na.omit(sample_data)
  
  sample_proportions <- table(sample_data_n$y)/nrow(sample_data_n)
  
  
  file_name <- paste0("C:\\...\\Imputed data\\sample proportions\\s5\\n", n, "\\RR", response_rate, "\\CCA_i_", i, "_n_", n, "_RR_", response_rate, ".rds")
  
  saveRDS(sample_proportions, file = file_name)
  
}


