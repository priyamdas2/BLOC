clear all; clc;
clearvars;
addpath('./pdlc/');
addpath('./WeiZhaoMCP');
addpath('./BLOC/');
addpath('./helper funs/');
% addpath('./Simulation data/');

p = 200; % covariance dimension 100/200
method_num = 3; % 1 = Block matrix, 2 = Toeplitz, 3 = Banded
SIGMA = cov_model(p, method_num);

Ns = 100; % [50,100]
Num_exps = 10;
tol = 1e-3; % sparsity detection

% ============================================================
% Wei-Zhao MCP covariance-estimator settings
% ============================================================
%
% Direct covariance-matrix version of Wei-Zhao MCP.
%
% For comparison:
%   - use the SAME lambda grid,
%   - use the SAME repeated train/test CV splitting,
%   - use tau = 0 so the smooth loss is the Frobenius loss,
%   - optimize directly over the covariance matrix.
%
num_CV = 5;
num_lambda = 15;
lambdas = logspace(-3, 0, num_lambda);

tau = 0;
a_val = 2;
phi0 = 1;
gamma_bt = 2;

for iN = 1:length(Ns)

    output_WeiZhaoMCP_Cov = nan(Num_exps, 8);

    for k = 1:Num_exps

        n = Ns(iN);
        [n, k]

        % Exact same seed-controlled data generation structure
        rng(k)
        r = mvnrnd(zeros(p,1), SIGMA, n);

        % Full-sample covariance matrix
        S = cov(r);
        
        tStart = tic;

        % ------------------------------------------------------------
        % Lambda selection using the same CV structure as BLOC
        % ------------------------------------------------------------
        n2 = round(n/log(n)); % test
        n1 = n - n2;          % train

        FroErr = nan(num_CV, num_lambda);

        for kk = 1:num_CV

            J = randperm(n);

            % Covariance matrices for training and testing
            S1 = cov(r(J(1:n1),:));
            S2 = cov(r(J(n1+1:end),:));

            for l = 1:num_lambda

                lambda_val = lambdas(l);

                fprintf(['=> WeiZhaoMCP_Cov: exp %d, CV %d, ' ...
                         'lambda grid %d.\n'], k, kk, l);

                X_cv = WeiZhaoMCP_Cov( ...
                    S1, tau, phi0, gamma_bt, ...
                    lambda_val, a_val);

                FroErr(kk,l) = norm(X_cv - S2, 'fro');
            end
        end

        [~, mi] = min(sum(FroErr,1));
        lambda_opt = lambdas(mi);

        fprintf(['=> WeiZhaoMCP_Cov exp %d, ' ...
                 'optimal lambda: %.6f.\n'], ...
                 k, lambda_opt);

        % ------------------------------------------------------------
        % Final fit on the full-sample covariance matrix
        % ------------------------------------------------------------

        X = WeiZhaoMCP_Cov( ...
            S, tau, phi0, gamma_bt, ...
            lambda_opt, a_val);

        comp_time = toc(tStart);

        % ------------------------------------------------------------
        % Performance measures
        % X is the Wei-Zhao covariance estimate
        % ------------------------------------------------------------
        ErrFro = norm(X - SIGMA, 'fro');
        ErrSpe = norm(X - SIGMA);
        err_abs = mean(abs(X(:) - SIGMA(:)));

        % Convert to correlation only for support-recovery metrics,
        % consistent with the other covariance-estimator comparisons.
        X_corr = cov2corr(X);

        [TPR, FPR, MCC] = ...
            TPR_FPR_MCC(X_corr, SIGMA, tol);

        PD = min(real(eig(X))) > -1e-8;

        output = [ErrFro, ErrSpe, err_abs, ...
                  TPR, FPR, MCC, PD, comp_time];

        output_WeiZhaoMCP_Cov(k,:) = output;

    end

    filename = sprintf( ...
        ['outputs/output_FroSpecAbsTprFprMccPDTime_method_%d_' ...
         'p_%d_n_%d_NumExp_%d_WeiZhaoMCP_Cov.csv'], ...
         method_num, p, n, Num_exps);

    csvwrite(filename, output_WeiZhaoMCP_Cov);

end