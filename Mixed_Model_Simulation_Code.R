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

#Clearing environment
rm(list=ls())

#Control statements for model fitting
cl_1 <- lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000,
                 msMaxEval = 1000, opt = c("optim"), optimMethod = "BFGS",
                 returnObject = FALSE)

cl_2 <- lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, 
  msMaxEval = 1000, opt = "nlminb",
  returnObject = FALSE)


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
           sigma_sq = 1,
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
                                   Sim_Sigma_Sq = sigma_sq,
                                   Sim_Residual = as.numeric(arima.sim(model = list(order = corr_mod_order, ar = ar_coef, ma = ma_coef), sd = sigma_sq, n = n_time)),
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
                                      Sim_Residual = as.numeric(arima.sim(model = list(order = corr_mod_order, ar = ar_coef, ma = ma_coef), sd = sigma_sq, n = n_time)),
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
                                  Sim_Sigma_Sq = sigma_sq,
                                  Sim_Residual = as.numeric(arima.sim(model = list(order = corr_mod_order, ar = ar_coef, ma = ma_coef), sd = sigma_sq, n = n_time)),
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
          relocate(matches("_SD$|_Corr$|Sigma_Sq$|_Coef$|_Order$"), .after = "Sim_Fixed_Time_Cb")        
        
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
                           data_set
                           ){
  
  poly_order <- match.arg(poly_order)
  
  #Fixed effects and random effects formulas based on poly_order
  model_formulas <- switch(poly_order,
                           "Linear" = list(fixef_formula = "Sim_Outcome ~ 1 + Time",
                                           ranef_formula = "~ 1 + Time | Subject_ID"),
                           
                           "Quadratic" = list(fixef_formula = "Sim_Outcome ~ 1 + Time + I(Time^2)",
                                              ranef_formula = "~ 1 + Time + I(Time^2) | Subject_ID"),
                           
                           "Cubic" = list(fixef_formula = "Sim_Outcome ~ 1 + Time + I(Time^2) + I(Time^3)",
                                          ranef_formula = "~ 1 + Time + I(Time^2) + I(Time^3) | Subject_ID")
                           
                           )
  
  #Extracting fixed and random effects formulas as variables
  fixef_formula <- model_formulas$fixef_formula
  ranef_formula <- model_formulas$ranef_formula
  
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
            lme(fixed = as.formula(fixef_formula), 
                 random = as.formula(ranef_formula),
                 correlation = modeled_corr, 
                 data = data_set,
                 #Using ML for LRT tests
                 method = "REML",
                 control = cl_1)),
            
            error = function(e){
              #Check if error message matches the target
              if (str_detect(e$message, pattern = "Coefficient matrix not invertible")) {
                suppressMessages(message("Target error detected: ", e$message, " - Retrying with second control."))
            
            #Retrying without control
            retry_mod <- tryCatch(
              suppressMessages(
                lme(fixed = as.formula(fixef_formula),
                    random = as.formula(ranef_formula),
                    correlation = modeled_corr,
                    data = data_set,
                    method = "REML",
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


#Number of simulated data sets for each condition
n_data_sims <- 3

#Plan for parallel processing
plan(multisession, workers = 6)
options(future.rng.onMisuse = "ignore")

tic()
simulations <- 
  #Obtaining combinations of each condition: Subject_n, Time_n, Sim_Ranef_Dist, Sim_Corr_Struct,
  expand_grid(Sim_Poly_Order = c("Linear", "Quadratic", "Cubic"),
              Subject_n = c(25,50), 
              Time_n = c(6,8),
              Sim_Corr_Struct = c("ID", "AR(1)","ARMA(1,1)")) %>% 
  mutate(Sim_Condition_ID = row_number(), .before = Sim_Poly_Order) %>% 
  mutate(Sim_Corr_Order = case_match(Sim_Corr_Struct,
                                     "ID" ~ list(c(0,0,0)),
                                    "AR(1)" ~ list(c(1,0,0)),
                                    "ARMA(1,1)" ~ list(c(1,0,1))),
         Sim_AR_Coef = case_match(Sim_Corr_Struct,
                                  "ID" ~ list(NULL),
                                  c("AR(1)", "ARMA(1,1)") ~ list(0.5)),
         Sim_MA_Coef = case_match(Sim_Corr_Struct,
                                  "ID" ~ list(NULL),
                                  "ARMA(1,1)" ~ list(0.5))) %>% 
  #Repeating each condition however many times we want
  uncount(n_data_sims) %>% 
  mutate(Sim_Data_ID = row_number(), .before = Sim_Condition_ID) %>% 
  #Simulating data based on conditions in each row
  mutate(Sim_Data = future_pmap(list(Sim_Poly_Order,
                              Subject_n,
                              Time_n,
                              Sim_Ranef_Dist,
                              Sim_Corr_Order,
                              Sim_AR_Coef,
                              Sim_MA_Coef),
                         ~sim_function(poly_order = ..1,
                                      n_subject = ..2,
                                      n_time = ..3,
                                      ranef_dist = ..4,
                                      corr_mod_order = ..5,
                                      ar_coef = ..6,
                                      ma_coef = ..7))) %>% 
  #Adding on the modeling conditions
  #Modeled correlation structure
  mutate(Model_Corr_Struct = list(c("ID", "AR(1)","ARMA(1,1)"))) %>% 
  unnest(Model_Corr_Struct) %>% 
  #Modeled polynomial structure
  mutate(Model_Poly_Order = case_when(
                                  Sim_Poly_Order == "Linear" ~ list("Linear"),
                                  Sim_Poly_Order == "Quadratic" ~ list(c("Linear", "Quadratic")),
                                  Sim_Poly_Order == "Cubic" ~ list(c("Linear", "Quadratic", "Cubic"))
                                  )
         ) %>% 
  #Unnesting Model_Poly_Order to obtain make sure lower order polynomial terms are modeled for each simulated dataset
  unnest(Model_Poly_Order) %>% 
  #Running models on each data set according to correlation structure specified in Model_Corr_Struct column
  #Also extracting the apVar object from the model fit. This represents the approximate Variance-Covariance matrix
  #of the random effect estimates and can produce NAs or singularities without throwing a warning/error in model fit
  mutate(Model_Result = 
           future_pmap(list(Model_Poly_Order, Model_Corr_Struct, Sim_Data), 
                       ~ model_function(poly_order = ..1, corr_struct = ..2, data_set = ..3)),
         Model_apVar = map(Model_Result, "apVar")) 

toc()

dir <- "/Users/Tanner/Desktop/Keith and Tanner Correlation Specification Paper/Mixed_Model_Sims_"

#saveRDS(simulations, file = paste0(dir,format(Sys.Date(), "%Y_%m_%d")))



###Adding Model Performance Metrics
#Tabulating Proportions of Models that Ran with Errors/Convergence Issues
model_convergence_issues <-
  simulations %>% group_by(Subject_n, Time_n, Sim_Poly_Order, Model_Poly_Order, Sim_Corr_Struct, Model_Corr_Struct) %>%
  summarize(Convergence_Issue_Count = sum(map_lgl(Model_Result, ~ class(.x) == "character")),
            Convergence_Issue_Prop = round(Convergence_Issue_Count/n(), digits = 3), .groups = "drop")

#Tabulating Proportion of models that ran with issues in apVar
model_apvar_issues <-
  simulations %>% group_by(Subject_n, Time_n, Sim_Poly_Order, Model_Poly_Order, Sim_Corr_Struct, Model_Corr_Struct) %>%
  summarize(apVar_Issue_Count = sum(map_lgl(Model_apVar, ~ length(class(.x)) == 1 & "character" %in% class(.x))),
            apVar_Issue_Prop = round(apVar_Issue_Count/n(), digits = 3), .groups = "drop")


#Wald Test Statistic Function (used to calculate empirical Type I error rate for fixed effects)
#wald_function <- function(data_col_, model_col_, fixef_name_){

  wald_stat <- map2_dbl(data_col_, model_col_, function(data_,model_){

    #Fixed effect population parameter
    fixef_param <- data_ %>% dplyr::select(matches(paste0("Fixed_",fixef_name_))) %>%
      unique() %>% as.numeric()

    #Fixed effect model summary
    fixef_mod_summary <- as.data.frame(summary(model_)$tTable) %>% rownames_to_column(var = "Fixef")

    #Fixed effect model estimate
    fixef_est <- fixef_mod_summary %>% filter(str_detect(Fixef, pattern = fixef_name_)) %>% pull(Value)

    #Fixed effect model standard error
    fixef_se <- fixef_mod_summary %>% filter(str_detect(Fixef, pattern = fixef_name_)) %>% pull(Std.Error)

    stat <- (fixef_est - fixef_param)/fixef_se

    return(stat)

  }
   )

  return(wald_stat)
}


#Function that creates named vector of random effect variances
ranef_extractor <- function(model){
  
  ranef_estimates <- VarCorr(model) %>% as.data.frame.matrix()
  
  ranef_sd <- ranef_estimates %>% pull(StdDev) %>% as.numeric()
  
  names(ranef_sd) <- rownames(ranef_estimates)
  
  return(ranef_sd)
  
}



##Removing models that ran with issues for this data frame
sim_model_metrics <- simulations %>% 
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
  rename_with(~ str_replace_all(., c("\\^2" = "_Sq","\\^3" = "_Cb", "Rand_Residual" = "Sigma_Sq"))) %>% 
  rename_with(~ paste0(., "_SD"), .cols = contains("Rand"))



#Function to calculate relative bias for model estimates
rel_bias <- function(data_col, param_name_, model_col_){
  
  if(param_name_ %in% names(data_col)){
    
    param <- unique(data_col[[param_name_]])
    
    rel_bias <- (model_col_ - param)/param %>% as.double()
    
  } else {rel_bias = NA_real_}
  
  
  return(rel_bias)
}

##Calculating Relative Biases for fixed effects and random effect estimates
sim_rel_bias <- 
  sim_model_metrics %>% 
  dplyr::select(-c(Model_AIC:Model_BIC)) %>% 
  mutate(across(matches("Fixed|Rand|Sigma"),
                ~ map2_dbl(Sim_Data, .x, 
                           ~ rel_bias(data_col = .x, 
                                      param_name_ = paste0("Sim_", str_remove(cur_column(), "^Model_")), 
                                      model_col_ = .y)),
                .names = "Rel_Bias_{.col}")) %>% 
  rename_with(~str_remove(., pattern = "_Model"), starts_with("Rel_Bias"))
  





# ##Model Comparisons LRT for all models
# sim_model_lrt <-
#   sim_model_metrics %>% 
#   group_by(Sim_Data_ID) %>% 
#   mutate(
#     LRT = list(do.call(anova, Model_Result)),
#     LRT = map(LRT, ~ .x %>% 
#                 dplyr::select(-call) %>% 
#                 mutate(Sim_Corr_Struct = Sim_Corr_Struct, .before = Model) %>% 
#                 mutate(Sim_Data_ID = Sim_Data_ID, .before = 1) %>% 
#                 tibble()),
#     Correlations = list(unique(Model_Corr_Struct)),
#     LRT = map2(LRT, Correlations, ~ .x %>% mutate(Model_Corr_Struct = .y, .after = Sim_Corr_Struct))
#   ) %>% dplyr::select(-Correlations)
  
#Tabulating which models are the best fit to each dataset using AIC, AICc, and BIC
sim_model_comparisons <- 
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

sim_model_comparisons %>% 
  dplyr::summarize(AIC_Prop = mean(Sim_Model == AIC_Best_Model),
                   AICc_Prop = mean(Sim_Model == AICc_Best_Model),
                   BIC_Prop = mean(Sim_Model == BIC_Best_Model))


