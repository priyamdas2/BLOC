clear all; clc;
clearvars;
addpath('./pdlc/');
addpath('./WeiZhaoMCP');
addpath('./BLOC/');
addpath('./helper funs/');
% addpath('./Simulation data/');

p = 200; % covariance dimension 100/200
method_num = 1; % 1 = Block matrix, 2 = Toeplitz, 3 = Banded
SIGMA = cov_model(p, method_num);

Ns = [100]; % [50,100]
Num_exps = 10;
tol = 1e-3; % sparsity detection

% ============================================================
% Wei-Zhao MCP settings
% ============================================================
%
% The paper uses MCP with a = 2.
%
% For the direct head-to-head with BLOC below:
%   - use the SAME BLOC lambda grid,
%   - use the SAME BLOC repeated train/test CV splitting,
%   - use tau = 0 so the smooth loss is the Frobenius loss and the
%     comparison isolates optimization of the same MCP-penalized
%     correlation objective under positive definiteness.
%
% The paper states that lambda and tau were selected by five-fold CV,
% but it does not report the actual lambda/tau grids. The public Python
% class supplied to WeiZhaoMCP also takes phi0, gamma, lambda and tau as
% inputs rather than defining experiment-level defaults.
%
num_CV = 5;
num_lambda = 15;
lambdas = logspace(-3, 0, num_lambda);

tau = 0;       % matched-objective comparison with BLOC
a_val = 2;     % Wei-Zhao paper: MCP shape parameter
phi0 = 1;      % proximal-gradient initial quadratic coefficient
gamma_bt = 2;  % backtracking multiplier (>1)

for iN = 1:length(Ns)

    output_WeiZhaoMCP = nan(Num_exps, 8);

    for k = 1:Num_exps

        n = Ns(iN);
        [n, k]

        % Exact same seed-controlled data generation structure
        rng(k)
        r = mvnrnd(zeros(p,1), SIGMA, n);

        S = cov(r);
        STD_S = diag(S).^(0.5);
        S_corr = diag(1./STD_S) * S * diag(1./STD_S);
        
        tStart = tic;

        % ------------------------------------------------------------
        % Lambda selection using the same CV structure as BLOC
        % ------------------------------------------------------------
        n2 = round(n/log(n)); % test
        n1 = n - n2;          % train

        FroErr = nan(num_CV, num_lambda);

        for kk = 1:num_CV

            J = randperm(n);

            S1 = corr(r(J(1:n1),:));
            S2 = corr(r(J(n1+1:end),:));

            for l = 1:num_lambda

                lambda_val = lambdas(l);

                fprintf(['=> WeiZhaoMCP: exp %d, CV %d, ' ...
                         'lambda grid %d.\n'], k, kk, l);

                X_cv = WeiZhaoMCP(S1, tau, phi0, gamma_bt, ...
                                  lambda_val, a_val);

                FroErr(kk,l) = norm(X_cv - S2, 'fro');
            end
        end

        [~, mi] = min(sum(FroErr,1));
        lambda_opt = lambdas(mi);

        fprintf('=> WeiZhaoMCP exp %d, optimal lambda: %.6f.\n', ...
                k, lambda_opt);

        % ------------------------------------------------------------
        % Final fit on the full-sample correlation matrix
        % ------------------------------------------------------------

        X = WeiZhaoMCP(S_corr, tau, phi0, gamma_bt, ...
                       lambda_opt, a_val);
        comp_time = toc(tStart);

        % ------------------------------------------------------------
        % Performance measures
        % X is the Wei-Zhao correlation estimate Gamma-hat.
        % SIGMA has unit diagonal under the simulation models.
        % ------------------------------------------------------------
        ErrFro = norm(X - SIGMA, 'fro');
        ErrSpe = norm(X - SIGMA);
        err_abs = mean(abs(X(:) - SIGMA(:)));

        [TPR, FPR, MCC] = TPR_FPR_MCC(X, SIGMA, tol);

        PD = min(real(eig(X))) > -1e-8;

        output = [ErrFro, ErrSpe, err_abs, ...
                  TPR, FPR, MCC, PD, comp_time];

        output_WeiZhaoMCP(k,:) = output;
    end

    filename = sprintf( ...
        ['outputs/output_FroSpecAbsTprFprMccPDTime_method_%d_' ...
         'p_%d_n_%d_NumExp_%d_WeiZhaoMCP.csv'], ...
         method_num, p, n, Num_exps);

    csvwrite(filename, output_WeiZhaoMCP);
end
