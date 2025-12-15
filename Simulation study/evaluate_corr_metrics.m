function [TPR, FPR, MCC, RMSE, MAD] = evaluate_corr_metrics(C_est, C_true, threshold)
% EVALUATE_CORR_METRICS computes TPR, FPR, MCC, RMSE, and MAD for the
% off-diagonal entries of estimated and true correlation matrices.
%
% Inputs:
%   C_est     : Estimated correlation matrix (p x p)
%   C_true    : Ground truth correlation matrix (p x p)
%   threshold : (optional) threshold to define nonzero structure (default: 1e-6)
%
% Outputs:
%   TPR  : True Positive Rate
%   FPR  : False Positive Rate
%   MCC  : Matthews Correlation Coefficient
%   RMSE : Root Mean Squared Error of off-diagonal entries
%   MAD  : Mean Absolute Deviation of off-diagonal entries

    if nargin < 3
        threshold = 1e-6; % small value to treat as structural zero
    end

    % Get off-diagonal indices
    p = size(C_true, 1);
    off_diag_mask = ~eye(p);

    % Extract off-diagonal elements
    C_true_off = C_true(off_diag_mask);
    C_est_off = C_est(off_diag_mask);

    % Binary structure: 1 = edge present, 0 = no edge
    structure_true = abs(C_true_off) > threshold;
    structure_est  = abs(C_est_off)  > threshold;

    % Confusion matrix components
    TP = sum(structure_true & structure_est);
    TN = sum(~structure_true & ~structure_est);
    FP = sum(~structure_true & structure_est);
    FN = sum(structure_true & ~structure_est);

    % Metrics
    TPR = TP / (TP + FN + eps);  % Avoid division by zero
    FPR = FP / (FP + TN + eps);
    MCC = (TP * TN - FP * FN) / sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN) + eps);

    % Errors (numerical)
    diff = C_est_off - C_true_off;
    RMSE = sqrt(mean(diff.^2));
    MAD  = mean(abs(diff));
end
