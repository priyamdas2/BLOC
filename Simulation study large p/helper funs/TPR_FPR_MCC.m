function [TPR, FPR, MCC] = TPR_FPR_MCC(X, SIGMA, tol)
% TPR_FPR_MCC - Compute True Positive Rate, False Positive Rate, and MCC
%               comparing X vs true matrix SIGMA (off-diagonal only).
%
% Inputs:
%   X       : estimated matrix (p x p)
%   SIGMA   : true matrix (p x p)
%   tol     : threshold for defining nonzero (default = 1e-3)
%
% Outputs:
%   TPR     : True Positive Rate
%   FPR     : False Positive Rate
%   MCC     : Matthews Correlation Coefficient

    if nargin < 3
        tol = 1e-3; % default tolerance
    end

    % Ensure same size
    if ~isequal(size(X), size(SIGMA))
        error('X and SIGMA must have same dimension');
    end

    p = size(X,1);

    % Binary support matrices based on tolerance
    est_support  = abs(X) >= tol;
    true_support = abs(SIGMA) >= tol;

    % Ignore diagonals
    est_support(1:p+1:end)  = 0;
    true_support(1:p+1:end) = 0;

    % Flatten
    est_vec  = est_support(:);
    true_vec = true_support(:);

    % Confusion matrix counts
    TP = sum(est_vec == 1 & true_vec == 1);
    TN = sum(est_vec == 0 & true_vec == 0);
    FP = sum(est_vec == 1 & true_vec == 0);
    FN = sum(est_vec == 0 & true_vec == 1);

    % Metrics
    if (TP+FN) == 0
        TPR = NaN;
    else
        TPR = TP / (TP + FN);
    end

    if (FP+TN) == 0
        FPR = NaN;
    else
        FPR = FP / (FP + TN);
    end

    denom = sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    if denom == 0
        MCC = NaN;
    else
        MCC = (TP*TN - FP*FN) / denom;
    end
end
