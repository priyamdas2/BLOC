rm(list=ls())
setwd("U:/Sparse_Covariance_with_BLOC/Simulation study")
source("SpCov_supp_funs.R")
library(spcov)

### Reading data ###############################################################

# User-provided values
p <- 20
n <- 50
Ctype <- "Blockdiag"    # "Blockdiag", ("SparseUniform" not needed)
threshold <- 10^(-2)
num_dataset_rep <- 10

### Lambda grid 
log10lambdaLB <- -2
log10lambdaUB <- 1
NumLambdas <- 20
CVNumFolds <- 5
lambda_grid <- 10^seq(log10lambdaLB, log10lambdaUB, length.out = NumLambdas)
lambda_array <- lambda_grid
### Read data 
file_name <- sprintf("C_p_%d_n_%d_C_%s.csv", p, n, Ctype)
file_path <- file.path("Simulation data", file_name)
C_true <- as.matrix(read.csv(file_path, header = FALSE))

# penalty mask
file_name <- sprintf("P_cover_p_%d_n_%d_%s.csv", p, n, Ctype)
file_path <- file.path("Simulation data", file_name)
P <- as.matrix(read.csv(file_path, header = FALSE))

################################################################################
set.seed(1)


metric_names <- c("TPR", "FPR", "MCC", "RMSE", "MAD", "Comp_time")

# Create storage for final results and selected correlation matrices
results_final_all <- matrix(NA, nrow = num_dataset_rep, ncol = length(metric_names))
colnames(results_final_all) <- metric_names
C_final_list <- vector("list", num_dataset_rep)
All_best_lambda <- rep(NA, num_dataset_rep)


for (DataRep in 1:num_dataset_rep) {
  set.seed(DataRep)
  cat("Performing experiment number:", DataRep, "\n")
  # --- Load dataset ---
  file_name <- sprintf("X_p_%d_n_%d_C_%s_DataRep_%d.csv", p, n, Ctype, DataRep)
  file_path <- file.path("Simulation data", file_name)
  X <- read.csv(file_path, header = FALSE)
  Sigma_start <- var(X)
  p <- ncol(X)
  
  step.size <- 100
  
  # Initialize per-replicate storage
  results_mat <- matrix(NA, nrow = NumLambdas, ncol = length(metric_names[1:5]))
  colnames(results_mat) <- metric_names[1:5]
  EBICs <- rep(NA, NumLambdas)
  C_est_list <- vector("list", NumLambdas)
  
  # --- Model Fitting ---
  ptm <- proc.time()
  # --- Select best model for this replicate ---
  spcov_fit <- cv_loglikelihood(X, lambda_array, P, nfolds = CVNumFolds, step.size = 100, seed = DataRep)
  best_lambda <- spcov_fit$best_lambda
  All_best_lambda[DataRep] <- best_lambda
  # --- Fitting with best lambda ---
  S_whole <- cov(X)
  spcov_fit_best <- spcov(S_whole, S = S_whole, lambda = best_lambda * P, step.size = step.size)
  C_final_list[[DataRep]] <- spcov_fit_best$Sigma
  Comp_time <- (proc.time() - ptm)[3]
  C_est <- as.matrix(cov_to_cor(C_final_list[[DataRep]], threshold = threshold))
  Results <- evaluate_corr_metrics(C_est, C_true, threshold = threshold)
  results_final_all[DataRep, ] <- c(Results$TPR, Results$FPR, Results$MCC, Results$RMSE, Results$MAD, Comp_time)
}



output_folder <- "Simulation output"
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# Construct the file name
file_name <- sprintf("TPR_FPR_MCC_RMSE_MAD_time_custom_P_p_%d_n_%d_C_%s_NumReps_%d_SpCov.csv", 
                     p, n, Ctype, num_dataset_rep)
file_path <- file.path(output_folder, file_name)
write.csv(results_final_all, file = file_path, row.names = FALSE)
colMeans(results_final_all)
All_best_lambda