
cov_to_cor <- function(cov_mat, threshold = 1e-3) {
  # Validate input
  if (!is.matrix(cov_mat)) stop("Input must be a matrix.")
  if (nrow(cov_mat) != ncol(cov_mat)) stop("Covariance matrix must be square.")
  
  # Thresholding: zero out small entries
  cov_mat[abs(cov_mat) < threshold] <- 0
  
  # Compute standard deviations
  std_vec <- sqrt(diag(cov_mat))
  
  # Avoid division by zero
  if (any(std_vec == 0)) stop("Zero variance detected in diagonal.")
  
  # Compute correlation matrix
  cor_mat <- cov_mat / (std_vec %o% std_vec)  # Outer product
  
  return(cor_mat)
}

evaluate_corr_metrics <- function(C_est, C_true, threshold = 1e-3) {
  # Check that input matrices are square and of same size
  if (!is.matrix(C_est) || !is.matrix(C_true)) stop("Inputs must be matrices.")
  if (!all(dim(C_est) == dim(C_true))) stop("Matrices must have same dimensions.")
  
  p <- nrow(C_true)
  
  # Off-diagonal mask
  off_diag_mask <- matrix(TRUE, p, p)
  diag(off_diag_mask) <- FALSE
  
  # Extract off-diagonal elements
  C_true_off <- C_true[off_diag_mask]
  C_est_off  <- C_est[off_diag_mask]
  
  # Binary support structure based on threshold
  structure_true <- abs(C_true_off) > threshold
  structure_est  <- abs(C_est_off)  > threshold
  
  # Confusion matrix components
  TP <- sum(structure_true & structure_est)
  TN <- sum(!structure_true & !structure_est)
  FP <- sum(!structure_true & structure_est)
  FN <- sum(structure_true & !structure_est)
  
  eps <- .Machine$double.eps  # small constant
  
  # Performance metrics
  TPR  <- TP / (TP + FN + eps)
  FPR  <- FP / (FP + TN + eps)
  MCC  <- (as.numeric(TP) * TN - as.numeric(FP) * FN) / sqrt(as.numeric((TP + FP)) * (TP + FN) * (TN + FP) * (TN + FN) + eps)
  RMSE <- sqrt(mean((C_est_off - C_true_off)^2))
  MAD  <- mean(abs(C_est_off - C_true_off))
  
  return(list(TPR = TPR, FPR = FPR, MCC = MCC, RMSE = RMSE, MAD = MAD))
}



cv_loglikelihood <- function(X, lambda_array, P, nfolds = 5, step.size = 100, seed = 1) {
  set.seed(seed)
  n <- nrow(X)
  p <- ncol(X)
  folds <- sample(rep(1:nfolds, length.out = n))
  
  loglik_mat <- matrix(NA, nrow = length(lambda_array), ncol = nfolds)
  
  for (k in 1:nfolds) {
    test_idx <- which(folds == k)
    train_idx <- setdiff(1:n, test_idx)
    X_train <- X[train_idx, , drop = FALSE]
    X_test  <- X[test_idx, , drop = FALSE]
    
    S_test <- cov(X_test) * (length(test_idx) - 1) / length(test_idx)
    S_train <- cov(X_train) * (length(train_idx) - 1) / length(train_idx)
    
    for (i in seq_along(lambda_array)) {
      lambda <- lambda_array[i]
      Sigma_hat <- spcov(S_train, S = S_train, lambda = lambda * P, step.size = step.size)
      if (!is.null(Sigma_hat$Sigma)) {
        log_det <- tryCatch(log(det(Sigma_hat$Sigma)), error = function(e) NA)
        Sigma_inv <- tryCatch(solve(Sigma_hat$Sigma), error = function(e) NA)
        
        if (is.finite(log_det) && all(is.finite(Sigma_inv))) {
          ll <- -log_det - sum(Sigma_inv * S_test)
          loglik_mat[i, k] <- ll
        }
      }
    }
  }
  
  mean_loglik <- rowMeans(loglik_mat, na.rm = TRUE)
  best_lambda <- lambda_array[which.max(mean_loglik)]
  
  return(list(best_lambda = best_lambda, loglik = mean_loglik, matrix = loglik_mat))
}



compute_EBIC_cov <- function(Sigma_est, X, threshold = 1e-3, gamma = 0.5) {
  # Check input
  if (!is.matrix(Sigma_est) || !is.matrix(X)) stop("Inputs must be matrices.")
  if (nrow(Sigma_est) != ncol(Sigma_est)) stop("Sigma_est must be square.")
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Sample covariance with 1/n normalization
  S <- cov(X) * (n - 1) / n
  
  # Try Cholesky to compute log-determinant and inverse
  chol_result <- tryCatch({
    R <- chol(Sigma_est)
    log_det <- 2 * sum(log(diag(R)))
    Sigma_inv <- solve(Sigma_est)
    log_likelihood <- -n/2 * (log_det + sum(Sigma_inv * S))
    log_likelihood
  }, error = function(e) {
    return(Inf)
  })
  
  if (!is.finite(chol_result)) {
    return(Inf)
  }
  
  # Count non-zero off-diagonal elements
  off_diag_mask <- matrix(TRUE, p, p)
  diag(off_diag_mask) <- FALSE
  num_edges <- sum(abs(Sigma_est[off_diag_mask]) > threshold) / 2  # account for symmetry
  
  # EBIC formula
  ebic <- -2 * chol_result + num_edges * log(n) + 4 * gamma * log(choose(p, 2))
  
  return(ebic)
}