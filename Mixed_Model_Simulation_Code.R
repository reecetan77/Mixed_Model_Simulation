#Loading relevant libraries
library(tidyverse)
library(lme4)
library(nlme)
library(MASS)
library(VGAM)
library(qqplotr)
library(furrr)
library(tictoc)
library(AICcmodavg)
library(parallel)

#Clearing environment
rm(list=ls())

#Control statements for model fitting
cl_1 <- lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000,
                 msMaxEval = 1000, opt = c("optim"), optimMethod = "BFGS",
                 returnObject = FALSE, apVar = FALSE)

cl_2 <- lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, 
  msMaxEval = 1000, opt = "nlminb",
  returnObject = FALSE, apVar = FALSE)


cl_a <- lmerControl(optimizer = "bobyqa",
                    optCtrl=list(maxfun=1e9),
                    check.conv.grad = .makeCC("warning", tol = 1e-3, relTol = NULL),
                    check.conv.singular = .makeCC(action = "message", tol = 1e-9),
                    check.conv.hess = .makeCC(action = "warning", tol = 1e-6))

#Function to simulate data with random intercepts, linear, and quadratic time effects
sim_function <- 
  function(
           #Setting random seed for simulation 
           seed = NULL,
           #Number of Subjects and Time Points
           n_subject = 100,
           n_time = 10,
           #Polynomial order
           poly_order = c("Linear", "Quadratic", "Cubic"),
           #Fixed effects parameters
           fixed_intercept = 4,
           fixed_time = 2,
           fixed_time_sq = -0.5,
           fixed_time_cb = 0.05,
           #Random effect variances
           rand_intercept_sd = 2,
           rand_time_sd = 2,
           rand_time_sq_sd = 0.5,
           rand_time_cb_sd = 0.1,
           #Random effect correlations
           intercept_time_corr = 0.1,
           intercept_time_sq_corr = 0.1,
           intercept_time_cb_corr = 0.1,
           time_time_sq_corr = -0.1,
           time_time_cb_corr = 0.1,
           time_sq_time_cb_corr = 0.1,
           #Residual/Innovation sd
           resid_sd = 1,
           #Within-Cluster correlation structure
           #p: autoregressive order
           #d: order of differencing
           #q: moving average order
           corr_mod_order = c(0,0,0),
           #AR Coefficient
           ar_coef = NULL,
           #MA Coefficient
           ma_coef = NULL
  ){
      #Setting random seed
      set.seed(seed)
    
      # #Ensuring null values don't get passed to arima.sim function
      # ar_coef <- if (is.na(ar_coef)) numeric(0) else ar_coef
      # ma_coef <- if (is.na(ma_coef)) numeric(0) else ma_coef
    
      #Setting highest polynomial order
      poly_order <- match.arg(poly_order)
          
      #Simulating data
      fixed_random_effects <- switch(poly_order,
                              
        "Linear" = expand_grid(Subject_ID = 1:n_subject, Time = 1:n_time) %>% 
          #Adding fixed effects
          mutate(Sim_Fixed_Intercept = fixed_intercept,
                 Sim_Fixed_Time = fixed_time,
                 Sim_Fixed_Time_Sq = 0,
                 Sim_Fixed_Time_Cb = 0) %>% 
          #Adding random effects
          left_join(tibble(Subject_ID = 1:n_subject, 
                           faux::rnorm_multi(
                             #number of observations
                             n = n_subject,
                             #all of mean = 0
                             mu = c(Sim_Rand_Intercept = 0, Sim_Rand_Time = 0),
                             #sd for each random effect
                             sd = c(rand_intercept_sd, rand_time_sd),
                             #correlation matrix for random effects
                             r = c(1, intercept_time_corr,
                                   intercept_time_corr, 1)
                           )), by = "Subject_ID"),
        
        "Quadratic" = expand_grid(Subject_ID = 1:n_subject, Time = 1:n_time) %>% 
          #Adding fixed effects
          mutate(Sim_Fixed_Intercept = fixed_intercept,
                 Sim_Fixed_Time = fixed_time,
                 Sim_Fixed_Time_Sq = fixed_time_sq,
                 Sim_Fixed_Time_Cb = 0) %>% 
          #Adding random effects
          left_join(tibble(Subject_ID = 1:n_subject, 
                            faux::rnorm_multi(
                              #number of observations
                              n = n_subject,
                              #all of mean = 0
                              mu = c(Sim_Rand_Intercept = 0, Sim_Rand_Time = 0, Sim_Rand_Time_Sq = 0),
                              #sd for each random effect
                              sd = c(rand_intercept_sd, rand_time_sd, rand_time_sq_sd),
                              #correlation matrix for random effects
                              r = c(1, intercept_time_corr, intercept_time_sq_corr,
                                    intercept_time_corr, 1, time_time_sq_corr,
                                    intercept_time_sq_corr, time_time_sq_corr, 1)
                            )), by = "Subject_ID"),
        
        "Cubic" = expand_grid(Subject_ID = 1:n_subject, Time = 1:n_time) %>% 
          #Adding fixed effects
          mutate(Sim_Fixed_Intercept = fixed_intercept,
                 Sim_Fixed_Time = fixed_time,
                 Sim_Fixed_Time_Sq = fixed_time_sq,
                 Sim_Fixed_Time_Cb = fixed_time_cb) %>% 
          #Adding random effects
          left_join(tibble(Subject_ID = 1:n_subject, 
                           faux::rnorm_multi(
                             #number of observations
                             n = n_subject,
                             #all of mean = 0
                             mu = c(Sim_Rand_Intercept = 0, Sim_Rand_Time = 0, Sim_Rand_Time_Sq = 0, Sim_Rand_Time_Cb = 0),
                             #sd for each random effect
                             sd = c(rand_intercept_sd, rand_time_sd, rand_time_sq_sd, rand_time_cb_sd),
                             #correlation matrix for random effects
                             r = c(1, intercept_time_corr, intercept_time_sq_corr,intercept_time_cb_corr,
                                   intercept_time_corr, 1, time_time_sq_corr, time_time_cb_corr,
                                   intercept_time_sq_corr, time_time_sq_corr, 1, time_sq_time_cb_corr,
                                   intercept_time_cb_corr, time_time_cb_corr, time_sq_time_cb_corr, 1)
                           )), by = "Subject_ID")
        )
      
      #Error Terms
      #Model order c(p,d,q)
      #p: autoregressive order
      #d: order of differencing
      #q: moving average order
      
      #Combining fixed and random effects, adding error term outcome column
      sim_data <- fixed_random_effects %>% 
        group_by(Subject_ID) %>% 
        mutate(
        bind_cols(switch(poly_order,
                         
                 "Linear" = tibble(Sim_Rand_Intercept_SD = rand_intercept_sd,
                                   Sim_Rand_Time_SD = rand_time_sd,
                                   Sim_Rand_Time_Sq_SD = 0,
                                   Sim_Rand_Time_Cb_SD = 0,
                                   Sim_Rand_Intercept_Time_Corr = intercept_time_corr,
                                   Sim_Rand_Intercept_Time_Sq_Corr = NA,
                                   Sim_Rand_Intercept_Time_Cb_Corr = NA,
                                   Sim_Rand_Time_Time_Sq_Corr = NA,
                                   Sim_Rand_Time_Time_Cb_Corr = NA,
                                   Sim_Rand_Time_Sq_Time_Cb_Corr = NA,
                                   Sim_Resid_SD = resid_sd,
                                   Sim_Residual = as.numeric(arima.sim(model = list(order = corr_mod_order, ar = ar_coef, ma = ma_coef), sd = resid_sd, n = n_time)),
                                   Sim_Outcome = (Sim_Fixed_Intercept + Sim_Rand_Intercept) + 
                                     (Sim_Fixed_Time + Sim_Rand_Time)*Time + 
                                     Sim_Residual
                                  ),
                 
                 "Quadratic" = tibble(Sim_Rand_Intercept_SD = rand_intercept_sd,
                                      Sim_Rand_Time_SD = rand_time_sd,
                                      Sim_Rand_Time_Sq_SD = rand_time_sq_sd,
                                      Sim_Rand_Time_Cb_SD = 0,
                                      Sim_Rand_Intercept_Time_Corr = intercept_time_corr,
                                      Sim_Rand_Intercept_Time_Sq_Corr = intercept_time_sq_corr,
                                      Sim_Rand_Intercept_Time_Cb_Corr = NA,
                                      Sim_Rand_Time_Time_Sq_Corr = time_time_sq_corr,
                                      Sim_Rand_Time_Time_Cb_Corr = NA,
                                      Sim_Rand_Time_Sq_Time_Cb_Corr = NA,
                                      Sim_Resid_SD = resid_sd,
                                      Sim_Residual = as.numeric(arima.sim(model = list(order = corr_mod_order, ar = ar_coef, ma = ma_coef), sd = resid_sd, n = n_time)),
                                      Sim_Outcome = (Sim_Fixed_Intercept + Sim_Rand_Intercept) + 
                                        (Sim_Fixed_Time + Sim_Rand_Time)*Time + 
                                        (Sim_Fixed_Time_Sq + Sim_Rand_Time_Sq)*Time^2 + 
                                        Sim_Residual
                                      ),
                 
                 "Cubic" = tibble(Sim_Rand_Intercept_SD = rand_intercept_sd,
                                  Sim_Rand_Time_SD = rand_time_sd,
                                  Sim_Rand_Time_Sq_SD = rand_time_sq_sd,
                                  Sim_Rand_Time_Cb_SD = rand_time_cb_sd,
                                  Sim_Rand_Intercept_Time_Corr = intercept_time_corr,
                                  Sim_Rand_Intercept_Time_Sq_Corr = intercept_time_sq_corr,
                                  Sim_Rand_Intercept_Time_Cb_Corr = intercept_time_cb_corr,
                                  Sim_Rand_Time_Time_Sq_Corr = time_time_sq_corr,
                                  Sim_Rand_Time_Time_Cb_Corr = time_time_cb_corr,
                                  Sim_Rand_Time_Sq_Time_Cb_Corr = time_sq_time_cb_corr,
                                  Sim_Resid_SD = resid_sd,
                                  Sim_Residual = as.numeric(arima.sim(model = list(order = corr_mod_order, ar = ar_coef, ma = ma_coef), sd = resid_sd, n = n_time)),
                                  Sim_Outcome = (Sim_Fixed_Intercept + Sim_Rand_Intercept) + (Sim_Fixed_Time + Sim_Rand_Time)*Time + 
                                    (Sim_Fixed_Time_Sq + Sim_Rand_Time_Sq)*Time^2 + (Sim_Fixed_Time_Cb + Sim_Rand_Time_Cb)*Time^3 +
                                    Sim_Residual
                                  )
        ))) %>% 
          ungroup() %>% 
          #Appending the AR and MA coefficients to the dataframe
          mutate(Sim_AR_Coef = list(ar_coef),
                 Sim_MA_Coef = list(ma_coef),
                 Sim_Corr_Order = list(corr_mod_order)) %>% 
          #Rearranging columns to put generative model parameters 
          #(fixed effects, random effect sds/correlations and ar/ma coefs) before the generated values
          relocate(matches("_SD$|_Corr$|Resid_SD$|_Coef$|_Order$"), .after = "Sim_Fixed_Time_Cb")        
        
      return(sim_data)
  }

#Function to run mixed-effects model using nlme
#Using tryCatch to catch errors/warnings/convergence issues
#1/24/2025 - Fitting models using ARMA(1,1) often resulted in a "Coefficient matrix not invertible" error
#when using cl_1. Sometimes trying cl_2 would fix the issue, but sometimes an error would still present. 
#This function returns the model when no errors are encountered. If it runs into the above error message,
#cl_2 is used. If different error is encountered at any stage, the message is returned.

model_function <- function(poly_order = c("Linear", "Quadratic", "Cubic"),
                           corr_struct = c("ID", "AR(1)","ARMA(1,1)"),
                           data_set){
  
  poly_order <- match.arg(poly_order)
  
  #Predefined formulas
  fixef_formula <- list(
    Linear = Sim_Outcome ~ 1 + Time,
    Quadratic = Sim_Outcome ~ 1 + Time + I(Time^2),
    Cubic = Sim_Outcome ~ 1 + Time + I(Time^2) + I(Time^3)
  )
  
  ranef_formula <- list(
    Linear = ~ 1 + Time | Subject_ID,
    Quadratic = ~ 1 + Time + I(Time^2) | Subject_ID,
    Cubic = ~ 1 + Time + I(Time^2) + I(Time^3) | Subject_ID
  )
  
  #Correlation structure to be included in the model
  corr_struct <- match.arg(corr_struct)
  
  modeled_corr <- switch(corr_struct,
                        "ID" = NULL,
                        "AR(1)" = corAR1(form = ~ Time | Subject_ID),
                        "ARMA(1,1)" = corARMA(form = ~ Time | Subject_ID, p = 1, q = 1)
                        )
  
  
  sim_mod <- tryCatch(
            #Attempt to fit model with the control statement
    suppressMessages(
            lme(fixed = fixef_formula[[poly_order]], 
                 random = ranef_formula[[poly_order]],
                 correlation = modeled_corr, 
                 data = data_set,
                 method = "ML",
                 control = cl_1)),
            
            error = function(e){
              #Check if error message matches the target
              if (str_detect(e$message, pattern = "Coefficient matrix not invertible")) {
                suppressMessages(message("Target error detected: ", e$message, " - Retrying with second control."))
            
            #Retrying without control
            retry_mod <- tryCatch(
              suppressMessages(
                lme(fixed = fixef_formula[[poly_order]],
                    random = ranef_formula[[poly_order]],
                    correlation = modeled_corr,
                    data = data_set,
                    method = "ML",
                    control = cl_2)),
                
                error = function(e2){
                  suppressMessages(message("Retry with second control failed: ", e2$message, " - Returning NULL."))
                  return(e2$message) # Return error message from retry
                }
              )
              } else {
                #Return the original error message for non-target errors
                suppressMessages(message("Non-target error detected: ", e$message))
                return(e$message)
              }
            }
            
  )
          
  return(sim_mod)
}



#Plan for parallel processing
plan(multisession, workers = detectCores() - 2)
options(future.rng.onMisuse = "ignore")

#Number of simulated data sets for each condition
n_data_sims <- 1

#Setting generative fixed effects
sim_fixed_intercept <- 5
sim_fixed_time <- 3
sim_fixed_time_sq <- -1
sim_fixed_time_cb <- -0.5

#Setting generative AR and MA coefficients
sim_ar_coef <- c(0.3, 0.7)
sim_ma_coef <- c(0.3, 0.7)

#Setting generative correlation structures
sim_corr_struct <- c("ID", "AR(1)","ARMA(1,1)")

#Setting generative random effect sd's
#Vector contains a "smaller" and "larger" value for each polynomial order
sim_rand_intercept_sd <- c(1,3)
sim_rand_time_sd <- c(1,3)
sim_rand_time_sq_sd <- c(1,3)
sim_rand_time_cb_sd <- c(1,3)

#Setting generative residual sd
#Vector contains a "smaller" and "larger" value
sim_resid_sd <- c(1,3)

#Creating smaller subset of simulation conditions and models to use for adding model metrics
#function to sample n simulative conditions
sample_n_of <- function(data, size, ...) {
  dots <- quos(...)
  
  group_ids <- data %>% 
    group_by(!!! dots) %>% 
    group_indices()
  
  sampled_groups <- sample(unique(group_ids), size)
  
  data %>% 
    filter(group_ids %in% sampled_groups)
}


#Parameter grid for simulations
sim_parameter_grid <- 
  #Obtaining combinations of each condition: Polynomial Order, Subject_n, Time_n,
  expand_grid(#Generative Polynomial Order
              Sim_Poly_Order = c("Linear", "Quadratic", "Cubic"),
              #Number of Subjects
              Subject_n = c(25,50),
              #Number of Time points
              Time_n = c(6,8)) %>% 
  #Adding fixed effects that are constant across all conditions
  mutate(Sim_Fixed_Intercept = sim_fixed_intercept,
         Sim_Fixed_Time = sim_fixed_time,
         Sim_Fixed_Time_Sq = ifelse(Sim_Poly_Order != "Linear", sim_fixed_time_sq, 0),
         Sim_Fixed_Time_Cb = ifelse(Sim_Poly_Order == "Cubic", sim_fixed_time_cb, 0)) %>% 
  #Adding random effect sds, residual sd and level-1 correlation structures
  mutate(Sim_Rand_Intercept_SD = list(sim_rand_intercept_sd),
         Sim_Rand_Time_SD = list(sim_rand_time_sd),
         Sim_Rand_Time_Sq_SD = ifelse(Sim_Poly_Order == "Linear", 0, list(sim_rand_time_sq_sd)),
         Sim_Rand_Time_Cb_SD = ifelse(Sim_Poly_Order == "Cubic", list(sim_rand_time_cb_sd), 0),
         Sim_Resid_SD = list(sim_resid_sd),
         Sim_Corr_Struct = list(sim_corr_struct)) %>% 
  #Unnesting list columns to obtain all combinations of generative parameters
  unnest(Sim_Rand_Intercept_SD) %>% 
  unnest(Sim_Rand_Time_SD) %>% 
  unnest(Sim_Rand_Time_Sq_SD) %>% 
  unnest(Sim_Rand_Time_Cb_SD) %>% 
  unnest(Sim_Resid_SD) %>% 
  unnest(Sim_Corr_Struct) %>% 
  #Adding generative correlation order and ar/ma coefficients 
  mutate(Sim_Corr_Order = case_match(Sim_Corr_Struct,
                                     "ID" ~ list(c(0,0,0)),
                                    "AR(1)" ~ list(c(1,0,0)),
                                    "ARMA(1,1)" ~ list(c(1,0,1))),
         Sim_AR_Coef = case_match(Sim_Corr_Struct,
                                  "ID" ~ list(NULL),
                                  c("AR(1)", "ARMA(1,1)") ~ list(sim_ar_coef)),
         Sim_MA_Coef = case_match(Sim_Corr_Struct,
                                  "ID" ~ list(NULL),
                                  "ARMA(1,1)" ~ list(sim_ma_coef))) %>% 
  #Unnesting AR and MA coefficients
  unnest(Sim_AR_Coef, keep_empty = TRUE) %>% 
  unnest(Sim_MA_Coef, keep_empty = TRUE) %>% 
  #Replacing NAs with NULL in AR/MA Coef
  mutate(Sim_AR_Coef = ifelse(is.na(Sim_AR_Coef), list(NULL), Sim_AR_Coef),
         Sim_MA_Coef = ifelse(is.na(Sim_MA_Coef), list(NULL), Sim_MA_Coef)) %>% 
  #Adding IDs to each generative condition
  mutate(Sim_Condition_ID = dplyr::row_number(), .before = Sim_Poly_Order) %>% 
  #Repeating each condition however many times we want
  uncount(n_data_sims) %>% 
  mutate(Sim_Data_ID = row_number(), .before = Sim_Condition_ID) 

#Simulating data based on grid parameters
tic()
sim_data <-
  sim_parameter_grid %>% 
  #This line samples first half of data sets 
  slice_head(prop = 0.1) %>% 
  #This line samples the second half
  #slice_tail(prop = 0.50) %>% 
  #Simulating data based on conditions in each row
  mutate(Sim_Data = pmap(list(
                              #Setting seed to the Sim_Data_ID
                              Sim_Data_ID,
                              #Poly order
                              Sim_Poly_Order,
                              #Number of subjects
                              Subject_n,
                              #Number of time points
                              Time_n,
                              #Fixed effects
                              Sim_Fixed_Intercept,
                              Sim_Fixed_Time,
                              Sim_Fixed_Time_Sq,
                              Sim_Fixed_Time_Cb,
                              #Random effects
                              Sim_Rand_Intercept_SD,
                              Sim_Rand_Time_SD,
                              Sim_Rand_Time_Sq_SD,
                              Sim_Rand_Time_Cb_SD,
                              #Residual sd
                              Sim_Resid_SD,
                              #Level-1 correlation structure
                              Sim_Corr_Order,
                              Sim_AR_Coef,
                              Sim_MA_Coef),
                         ~sim_function(seed = ..1,
                                      poly_order = ..2,
                                      n_subject = ..3,
                                      n_time = ..4,
                                      fixed_intercept = ..5,
                                      fixed_time = ..6,
                                      fixed_time_sq = ..7,
                                      fixed_time_cb = ..8,
                                      rand_intercept_sd = ..9,
                                      rand_time_sd = ..10,
                                      rand_time_sq_sd = ..11,
                                      rand_time_cb_sd = ..12,
                                      resid_sd = ..13, 
                                      corr_mod_order = ..14,
                                      ar_coef = ..15,
                                      ma_coef = ..16)))
toc()

sim_data_dir <- "/Users/Tanner/Desktop/Keith_and_Tanner_Correlation_Specification_Paper/Mixed_Model_Data_Sims"

#saveRDS(sim_data, file = paste0(sim_data_dir,"/Mixed_Model_Data_Sims_First_Half_", format(Sys.Date(), "%Y_%m_%d"), ".rds"))
saveRDS(sim_data, file = paste0(sim_data_dir,"/Mixed_Model_Data_Sims_Second_Half_", format(Sys.Date(), "%Y_%m_%d"), ".rds"))


tic()
#Modeling data using different combinations of polynomial orders and level-1 correlation structures
sim_data_models <- sim_data %>% 
  #Adding on the modeling conditions
  #Modeled correlation structure
  mutate(Model_Corr_Struct = list(c("ID", "AR(1)","ARMA(1,1)"))) %>% 
  unnest(Model_Corr_Struct) %>% 
  #Modeled polynomial structure
  mutate(Model_Poly_Order = list(c("Linear", "Quadratic", "Cubic"))) %>% 
  #Unnesting Model_Poly_Order to obtain make sure lower order polynomial terms are modeled for each simulated dataset
  unnest(Model_Poly_Order) %>% 
  #Adding ID to each model condition
  mutate(group_key = paste(Model_Corr_Struct, Model_Poly_Order, sep = "_")) %>% 
  mutate(Model_ID = dense_rank(factor(group_key, levels = unique(group_key))), .before = Model_Corr_Struct) %>% 
  dplyr::select(-group_key) %>% 
  #Running models on each data set according to correlation structure specified in Model_Corr_Struct column
  #Also extracting the apVar object from the model fit. This represents the approximate Variance-Covariance matrix
  #of the random effect estimates and can produce NAs or singularities without throwing a warning/error in model fit
  mutate(Model_Result = 
           future_pmap(list(Model_Poly_Order, Model_Corr_Struct, Sim_Data), 
                       ~ model_function(poly_order = ..1, corr_struct = ..2, data_set = ..3)
         #Model_apVar = map(Model_Result, "apVar"
                           )) 
toc()

sim_models_dir <- "/Users/Tanner/Desktop/Keith_and_Tanner_Correlation_Specification_Paper/Mixed_Model_Model_Sims"

saveRDS(sim_data_models, file = paste0(sim_models_dir,"/Mixed_Model_Model_Sims_", format(Sys.Date(), "%Y_%m_%d"), ".rds"))



###Adding Model Performance Metrics
#Tabulating Proportions of Models that Ran with Errors/Convergence Issues
model_convergence_issues <-
  sim_data_models %>% group_by(Subject_n, Time_n, Sim_Poly_Order, Model_Poly_Order, Sim_Corr_Struct, Model_Corr_Struct) %>%
  summarize(Convergence_Issue_Count = sum(map_lgl(Model_Result, ~ class(.x) == "character")),
            Convergence_Issue_Prop = round(Convergence_Issue_Count/n(), digits = 3), .groups = "drop")

#Tabulating Proportion of models that ran with issues in apVar
model_apvar_issues <-
  sim_data_models %>% group_by(Subject_n, Time_n, Sim_Poly_Order, Model_Poly_Order, Sim_Corr_Struct, Model_Corr_Struct) %>%
  summarize(apVar_Issue_Count = sum(map_lgl(Model_apVar, ~ length(class(.x)) == 1 & "character" %in% class(.x))),
            apVar_Issue_Prop = round(apVar_Issue_Count/n(), digits = 3), .groups = "drop")


#Function that creates named vector of random effect variances
ranef_extractor <- function(model){
  
  ranef_estimates <- VarCorr(model) %>% as.data.frame.matrix()
  
  ranef_sd <- ranef_estimates %>% pull(StdDev) %>% as.numeric()
  
  names(ranef_sd) <- rownames(ranef_estimates)
  
  return(ranef_sd)
  
}


##Removing models that ran with issues for this data frame
sim_model_metrics <- sim_data_models %>% 
  filter(map_lgl(Model_Result, ~ class(.x) != "character")) %>% 
  #Dropping irrelevant columns 
  dplyr::select(-Sim_Corr_Order, -Model_apVar) %>% 
  mutate(
         #Adding Information Criteria
         Model_AIC = map_dbl(Model_Result, ~ AIC(.x)),
         Model_AICc = map_dbl(Model_Result, ~ AICc(.x)),
         Model_BIC = map_dbl(Model_Result, ~ BIC(.x))
         ) %>% 
  mutate(
         #Adding Fixed Effects Estimates
         Model_Fixed = map(Model_Result, ~ fixef(.x)),
         #Adding Random Effect SDs
         Model_Rand = map(Model_Result, ~ ranef_extractor(.x))
         ) %>% 
  #Unpacking estimates into their own columns
  unnest_wider(Model_Fixed, names_sep = "_") %>% 
  unnest_wider(Model_Rand, names_sep = "_") %>% 
  #Tidying column names
  rename_with(~ str_remove_all(., "\\(|\\)|I\\(")) %>% 
  rename_with(~ str_replace_all(., c("\\^2" = "_Sq","\\^3" = "_Cb", "Rand_Residual" = "Resid_SD"))) %>% 
  rename_with(~ paste0(., "_SD"), .cols = contains("Model_Rand"))


#Tabulating which models are the best fit to each dataset using AIC, AICc, and BIC
sim_model_IC_comparisons <- 
  sim_model_metrics %>% 
  dplyr::select(-Sim_Data, -Model_Result, - Sim_AR_Coef, - Sim_MA_Coef, -c(Model_Fixed_Intercept:Model_Rand_Time_Cb_SD)) %>% 
  relocate(Sim_Poly_Order, .after = Sim_Corr_Struct) %>% 
  group_by(across(Sim_Data_ID:Sim_Poly_Order)) %>% 
  summarize(
    AIC_min = min(Model_AIC),
    AIC_Model_Corr_Struct = Model_Corr_Struct[which.min(Model_AIC)],
    AIC_Model_Poly_Order = Model_Poly_Order[which.min(Model_AIC)],
    
    AICc_min = min(Model_AICc),
    AICc_Model_Corr_Struct = Model_Corr_Struct[which.min(Model_AICc)],
    AICc_Model_Poly_Order = Model_Poly_Order[which.min(Model_AICc)],
    
    BIC_min = min(Model_BIC),
    BIC_Model_Corr_Struct = Model_Corr_Struct[which.min(Model_BIC)],
    BIC_Model_Poly_Order = Model_Poly_Order[which.min(Model_BIC)],
    .groups = "drop"
  ) %>% 
  unite(Sim_Model, Sim_Corr_Struct, Sim_Poly_Order, sep = ", ", remove = FALSE) %>% 
  unite(AIC_Best_Model, starts_with("AIC_Model"), sep = ", ", remove = FALSE) %>% 
  unite(AICc_Best_Model, starts_with("AICc_Model"), sep = ", ", remove = FALSE) %>% 
  unite(BIC_Best_Model, starts_with("BIC_Model"), sep = ", ", remove = FALSE) %>% 
  relocate(Sim_Model, AIC_Best_Model, AICc_Best_Model, BIC_Best_Model, .after = Sim_Poly_Order)


#Adding Bias and Precision Metrics for Bias and Precision
col_names <- names(sim_model_metrics %>% 
                     dplyr::select(matches("Fixed|Rand|Resid")))
#Pulls suffixes out
suffixes <- unique(str_extract(col_names, "_.*$"))

sim_model_bias_precision <- sim_model_metrics %>% 
  #Mutating across all columns that start with "Sim_" and have a suffix of "Fixed", "Rand" or "Resid
  mutate(across(
    all_of(names(.)[str_detect(names(.), paste0("^Sim*(", paste(suffixes, collapse = "|"),")$"))]),
    #Calculating difference and squared difference each column from its counterpart that starts with "Model"
    .fns = list(
      Diff = ~ sim_model_metrics[[sub("^Sim_", "Model_", cur_column())]] - .,
      Diff_Sq = ~ (sim_model_metrics[[sub("^Sim_", "Model_", cur_column())]] - .)^2
    ))) %>% 
  rename_with(~ str_remove(.x, "^Sim_"), matches("Diff|Diff_Sq")) %>% 
  #Computing Bias and Precision for each Sim_ID and Model_ID combination
  group_by(Sim_Condition_ID, Model_ID) %>% 
  #Averaging differences between simulation parameter and model estimate for bias
  dplyr::mutate(across(ends_with("_Diff"), ~ mean(.x, na.rm = TRUE), .names = "{.col}_Bias")) %>% 
  #Taking sqrt of average squared differences between simulation parameter and model estimate for precision
  dplyr::mutate(across(ends_with("_Diff_Sq"), ~ sqrt(mean(.x, na.rm = TRUE)), .names = "{.col}_Precision")) %>% 
  #Modifying column names
  rename_with(~ str_remove(.x, "_Diff"), contains("Diff_Bias")) %>% 
  rename_with(~ str_remove(.x, "_Diff_Sq"), contains("Diff_Sq_Precision")) %>% 
  #Selecting relevant columns
  dplyr::select(c(Sim_Condition_ID:Time_n, Sim_Resid_SD, Sim_Corr_Struct, Model_ID:Model_Poly_Order), 
                matches("Bias|Precision")) %>% 
  #Rearranging to get bias and precision next to each other for each term
  dplyr::select(c(Sim_Condition_ID:Model_Poly_Order), 
                matches(paste0(rep(str_remove(suffixes, "^_"), each = 2), c("_Bias","_Precision")))) %>% 
  distinct() %>% 
  ungroup()
  
  
#Adding Model p-values for fixed effects t-tests
#Function to extract p-values and put them in table that follows naming convention
p_value_function <- function(model_){
  
p_vals <- summary(model_)$tTable %>% 
  data.frame() %>% 
  rownames_to_column(var = "Term") %>% 
  mutate(Term = case_match(Term,
                           "(Intercept)" ~ "Model_Fixed_Intercept",
                           "Time" ~ "Model_Fixed_Time",
                           "I(Time^2)" ~ "Model_Fixed_Time_Sq",
                           "I(Time^3)" ~ "Model_Fixed_Time_Cb")) %>% 
  dplyr::select(Term, p.value) %>% 
  pivot_wider(names_from = Term, values_from = p.value) %>% 
  rename_with(~ paste0(.x,"_p_value"), everything())
  
return(p_vals)
}

#Data frame to hold p-values
sim_model_p_values <-
  sim_model_metrics %>% 
  dplyr::select(Sim_Data_ID:Model_Result) %>% 
  mutate(p_values = map(Model_Result, ~p_value_function(.x))) %>% 
  unnest_wider(p_values) %>% 
  dplyr::select(-Sim_Data,-Model_Result)

  

