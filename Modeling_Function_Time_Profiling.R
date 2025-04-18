library(tidyverse)
library(qs)
library(nlme)
library(tictoc)
library(glmmTMB)
library(parallel)
library(furrr)

cl_1 <- lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000,
                   msMaxEval = 1000, opt = c("optim"), optimMethod = "BFGS",
                   returnObject = FALSE, apVar = FALSE)

cl_2 <- lmeControl(maxIter = 1000, msMaxIter = 1000, niterEM = 1000, 
                   msMaxEval = 1000, opt = "nlminb",
                   returnObject = FALSE, apVar = FALSE)

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

#Loading sim param grid for simulation conditions
sim_grid <- qread("/Users/Tanner/Desktop/Keith_and_Tanner_Correlation_Specification_Paper/Mixed_Model_Sim_Grid/Sim_Parameter_Grid.qs")

#To test the extreme case, let's pull a data set with simulated using cubic polynomial
#and ARMA(1,1) autocorrelation

sim_grid %>% 
  filter(Sim_Poly_Order == "Cubic", 
         Sim_Corr_Struct == "ARMA(1,1)", 
         Subject_n == 50,
         Time_n == 8)

#Let's using condition 1351
sim_condition_id <- 1265

setwd("/Users/Tanner/Desktop/Keith_and_Tanner_Correlation_Specification_Paper/Mixed_Model_Model_Grids")
files <- list.files()

#Pulling in the file with the condition we want
sim_condition <- qread(files[str_ends(files, pattern = paste0("_",sim_condition_id,".qs"))])


#Pulling a data set from this simulated condition
sim_data <- sim_condition$Sim_Data[[1]]


#Trying out different modeling functions and assessing runtime

# #Standard lme
# tic()
# profvis(
# lme(fixed = Sim_Outcome ~ 1 + Time + I(Time^2) + I(Time^3), 
#     random = ~ 1 + Time + I(Time^2) + I(Time^3) | Subject_ID,
#     correlation = corARMA(form = ~ Time | Subject_ID, p = 1, q = 1), 
#     control = cl_1,
#     data = sim_data)
# )
# toc()

#sim_model_grid <-qread(sim_data_chunk_) 

test_model_dir <- "/Users/Tanner/Desktop/Keith_and_Tanner_Correlation_Specification_Paper/Single_Condition_Model_Test"
setwd(test_model_dir)

#Plan for parallel processing
plan(multicore, workers = detectCores() - 1)
options(future.rng.onMisuse = "ignore")

tic()
#Computing Models
future_pwalk(list(sim_condition$Model_Poly_Order, sim_condition$Model_Corr_Struct, 
                                     sim_condition$Sim_Data, rownames(sim_condition)),
                                ~ model_function(poly_order = ..1, corr_struct = ..2, data_set = ..3) %>% 
                                  qsave(file = paste0("Test_Sim_Condition_", sim_condition_id,"_Mod_", ..4, ".qs"),
                                        preset = "uncompressed")
                                )
toc()
