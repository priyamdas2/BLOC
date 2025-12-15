function [CorrEst_stack, SigmaEst_stack] = BLOC_SpCovMCplus(X_raw, P, lambda_array, gamma_array, doparallel, varargin)

% Standardize the matrix
X_means = mean(X_raw, 1);   % 1 x p
X_sds   = std(X_raw, 0, 1); % 1 x p, normalize by (n-1)
X_scaled = (X_raw - X_means) ./ X_sds;

% Starting point
Corr_start = corr(X_scaled);
NumLambdas = length(lambda_array);
NumGammas = length(gamma_array);
p = size(X_raw,2);
CorrEst_stack = zeros(p, p, NumLambdas, NumGammas);
SigmaEst_stack = zeros(p, p, NumLambdas, NumGammas);
DD = diag(X_sds);
start_vals = nan(NumLambdas,1);

for kk = 1:NumLambdas
    lambda = lambda_array(kk);
    for jj = 1:NumGammas
        fprintf('MC+: Performing CV: (%d,%d) / (%d,%d) \n', kk, jj, NumLambdas, NumGammas);
        gamma = gamma_array(jj);
        objFun = @(Sigma)-log_likelihood_gaussian_MCplus(X_scaled, Sigma, lambda, gamma, P);
        start_vals(kk) = objFun(Corr_start);
        if(doparallel == 1)
            [C_opt, fval, comp_time] = BLOCparallel(objFun, Corr_start, varargin{:});
        else
            [C_opt, fval, comp_time] = BLOC(objFun, Corr_start, varargin{:});
        end
        CorrEst_stack(:,:,kk,jj) = C_opt;
        SigmaEst_stack(:,:,kk,jj)= DD *C_opt* DD;
    end
end
end
