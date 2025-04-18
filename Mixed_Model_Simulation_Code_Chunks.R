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
library(qs)

cl_1 <- lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000,
                   msMaxEval = 1000, opt = c("optim"), optimMethod = "BFGS",
                   returnObject = FALSE, apVar = FALSE)

cl_2 <- lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, 
                   msMaxEval = 1000, opt = "nlminb",
                   returnObject = FALSE, apVar = FALSE)



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




data_dir <- "/Users/Tanner/Desktop/Keith_and_Tanner_Correlation_Specification_Paper/Mixed_Model_Data_Sims"

model_dir <- "/Users/Tanner/Desktop/Keith_and_Tanner_Correlation_Specification_Paper/Mixed_Model_Model_Sims"

model_grid_dir <- "/Users/Tanner/Desktop/Keith_and_Tanner_Correlation_Specification_Paper/Mixed_Model_Model_Grids"

#Files to be completed
sim_data_files <- list.files(model_grid_dir)


#Plan for parallel processing
plan(multicore, workers = detectCores() - 1)
options(future.rng.onMisuse = "ignore")


setwd(model_grid_dir)


tic()
imap(sim_data_files, function(sim_data_chunk_, sim_data_chunk_id_){
  
  #Reading in chunk of model grid
  sim_model_grid <-qread(sim_data_chunk_) 
  
  #Computing Models
  sim_model_result <- future_pmap(list(sim_model_grid$Model_Poly_Order, sim_model_grid$Model_Corr_Struct, 
                                         sim_model_grid$Sim_Data),
                                    ~ model_function(poly_order = ..1, corr_struct = ..2, data_set = ..3))

  
  
  #saving each result
  qsave(sim_model_result, 
        file = paste0(model_dir,"/Sim_Model_Chunk_Sim_Condition_",sim_data_chunk_id_,".qs"),
        preset = "fast",
        nthreads = 4)
  
  rm(sim_model_result)
  gc()
  return(TRUE)
  
})
toc()


                   