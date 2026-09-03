clear all; clc;%close all;
clearvars;
addpath('./BLOC/');
%addpath('./pdlc/');
addpath(genpath('./manopt/'));
addpath('./supp funs/');

% ============================================================
% Simulation design
% ============================================================

p = 20; % low-dimensional discontinuous-loss experiment
method_num = 3; % 1 = Block matrix, 2 = Toeplitz, 3 = Banded
SIGMA = cov_model_matDepthsim(p, method_num);

penalty = 2; % 1 = SCAD, 2 = MCP
n = 200;

reps_to_do = 1:10; % any subset of {1,2,...,10}
tol = 1e-5; % sparsity detection

% Gross contamination
contam_prop = 0.10;
num_outliers = round(contam_prop*n);
outlier_scale = 3;

% ============================================================
% Fixed random directions U
% ============================================================
% IMPORTANT:
% Later competitor scripts should copy this block EXACTLY.

num_U = 500; % e.g. 100 or 500
U_seed = 246810;

rng(U_seed,'twister');

U = randn(p, num_U);
U = U ./ sqrt(sum(U.^2,1)); % each column has norm 1

% Gaussian conversion factor:
% beta = [Phi^{-1}(0.75)]^2
beta = norminv(0.75)^2;

% ============================================================
% Cross-validation settings
% ============================================================

n2 = round(n/log(n)); % Test
n1 = n - n2;          % Train

num_CV = 5;
num_lambda = 15;

lamdas = logspace(-3, 0, num_lambda);

C0 = eye(p);

% ============================================================
% BLOC parameters (CV)
% ============================================================

params.MaxTime = 200;
params.MaxRuns = 1;
params.MaxIter = 1000;
params.TolFun1 = 1e-6;
params.TolFun2 = 1e-4;
params.phi = 1e-6;
params.DisplayUpdate = 1;
params.DisplayEvery = 10;

% ============================================================
% BLOC parameters (final)
% ============================================================

paramsFinal.MaxTime = 3600;
paramsFinal.MaxRuns = 20;
paramsFinal.MaxIter = 1000;
paramsFinal.TolFun1 = 1e-6;
paramsFinal.TolFun2 = 1e-4;
paramsFinal.phi = 1e-6;
paramsFinal.DisplayUpdate = 1;
paramsFinal.DisplayEvery = 10;

% ============================================================
% Separate reproducibility seeds
% ============================================================
%
% These are kept separate so optimizer randomness cannot alter
% the generated data or CV splits.
%
% IMPORTANT:
% Later competitor scripts should use the SAME data seed and
% SAME CV_seed_base. Optimizer-specific seeds can be different.
%

CV_seed_base = 100000;

BLOC_CV_seed_base = 200000;
BLOC_Final_seed_base = 300000;

if ~exist('outputs','dir')
    mkdir('outputs');
end

% ============================================================
% Simulation replications
% ============================================================

parfor k = reps_to_do

    fprintf('=> Performing exp number: %d.\n', k);

    % --------------------------------------------------------
    % DATA GENERATION
    % KEEP THIS BLOCK IDENTICAL FOR ALL METHODS
    % --------------------------------------------------------

    rng(k,'twister');

    % Clean Gaussian sample from the sparse truth
    r = mvnrnd(zeros(p,1), SIGMA, n);

    % Replace exactly 10% of observations by gross outliers.
    %
    % IMPORTANT:
    % Keep the ordering
    %
    %       mvnrnd -> randperm -> contamination
    %
    % identical in all future competitor scripts.
    %
    outlier_idx = randperm(n, num_outliers);

    
    

    r(outlier_idx,:) = ...
        mvnrnd(zeros(p,1), ...
        outlier_scale^2 * eye(p), ...
        num_outliers);

    % Precompute squared projections on the common fixed
    % directions for the final full-data fit.
    proj2_full = (r*U).^2;

    % --------------------------------------------------------
    % CV SPLITS
    % KEEP THIS BLOCK IDENTICAL FOR ALL METHODS
    % --------------------------------------------------------

    rng(CV_seed_base + k,'twister');

    CV_indices = zeros(num_CV,n);

    for kk = 1:num_CV
        CV_indices(kk,:) = randperm(n);
    end

    % --------------------------------------------------------
    % Lambda selection
    % --------------------------------------------------------
    %
    % The estimator itself is still obtained by minimizing
    %
    %   - matrix depth + SCAD/MCP penalty.
    %
    % However, lambda is selected using the Frobenius distance
    % between the training estimator and a ROBUST correlation
    % estimate obtained from the held-out observations.
    %
    % The robust held-out correlation is based on Kendall's tau:
    %
    %   R_jk = sin((pi/2) * tau_jk).
    %
    % Smaller CV error is better.
    %

    FroErr = nan(num_CV,num_lambda);

    for kk = 1:num_CV

        J = CV_indices(kk,:);

        % Training and held-out samples
        r1 = r(J(1:n1),:);
        r2 = r(J(n1+1:end),:);

        % Squared projections needed by the matrix-depth
        % objective on the training sample.
        proj2_train = (r1*U).^2;

        % ----------------------------------------------------
        % Robust correlation target from held-out observations
        % ----------------------------------------------------

        Tau2 = corr(r2, 'Type', 'Kendall');

        S2_robust = sin((pi/2) * Tau2);

        % Enforce exact unit diagonal numerically.
        S2_robust(1:p+1:end) = 1;

        % ----------------------------------------------------
        % Fit across lambda grid
        % ----------------------------------------------------

        for l = 1:num_lambda

            fprintf(['=> Exp %d, CV %d, ' ...
                     'lambda grid %d.\n'], ...
                     k, kk, l);

            lambda = lamdas(l);

            if penalty == 1

                % SCAD
                a = 3.7;

                objFun = @(Theta) ...
                    obj_depth_scad_corr_avg( ...
                    Theta, ...
                    proj2_train, ...
                    U, ...
                    beta, ...
                    lambda, ...
                    a);

            else

                % MCP
                gamma = 3;

                objFun = @(Theta) ...
                    obj_depth_mcp_corr_avg( ...
                    Theta, ...
                    proj2_train, ...
                    U, ...
                    beta, ...
                    lambda, ...
                    gamma);

            end

            % Reproducible BLOC-side randomness;
            % independent of data and CV generation.
            rng(BLOC_CV_seed_base + ...
                10000*k + 100*kk + l, ...
                'twister');

            [C_opt, fval, comp_time] = ...
                BLOCparallel_v2( ...
                objFun, C0, params); %#ok<ASGLU>

            % ------------------------------------------------
            % Robust held-out CV criterion
            % ------------------------------------------------
            %
            % Compare the fitted correlation matrix with the
            % Kendall-based robust correlation estimated from
            % the held-out observations.
            %

            FroErr(kk,l) = ...
                norm(C_opt - S2_robust, 'fro');

        end
    end

    % --------------------------------------------------------
    % Select lambda minimizing mean robust held-out error
    % --------------------------------------------------------

    mean_FroErr = mean(FroErr,1);

    [~, mi] = min(mean_FroErr);

    lamda_opt = lamdas(mi);

    fprintf(['=> Exp %d, optimal lambda: %.6f, ' ...
             'mean robust CV error: %.6f.\n'], ...
             k, ...
             lamda_opt, ...
             mean_FroErr(mi));

    % Optional diagnostic display
    disp('      lambda      mean robust CV error');
    disp([lamdas(:), mean_FroErr(:)]);

    % --------------------------------------------------------
    % Final fit on the full contaminated sample
    % --------------------------------------------------------

    if penalty == 1

        % SCAD
        a = 3.7;

        objFun = @(Theta) ...
            obj_depth_scad_corr_avg( ...
            Theta, ...
            proj2_full, ...
            U, ...
            beta, ...
            lamda_opt, ...
            a);

    else

        % MCP
        gamma = 3;

        objFun = @(Theta) ...
            obj_depth_mcp_corr_avg( ...
            Theta, ...
            proj2_full, ...
            U, ...
            beta, ...
            lamda_opt, ...
            gamma);

    end

    rng(BLOC_Final_seed_base + k,'twister');

    [C_final, fval, comp_time_final] = ...
        BLOCparallel_v2( ...
        objFun, C0, paramsFinal); %#ok<ASGLU>

    % --------------------------------------------------------
    % Performance measures
    % --------------------------------------------------------

    ErrFro = norm(C_final-SIGMA,'fro');

    ErrSpe = norm(C_final-SIGMA);

    err_abs = ...
        mean(abs(C_final(:)-SIGMA(:)));

    [TPR,FPR,MCC] = ...
        TPR_FPR_MCC( ...
        C_final, SIGMA, tol);

    PD = ...
        min(real(eig(C_final))) > 1e-8;

    % Achieved UNPENALIZED empirical matrix depth
    %
    % This is useful later when comparing BLOC with
    % patternsearch/fmincon on the same discontinuous loss.
    Depth = ...
        matrix_depth_fixedU_avg( ...
        C_final, ...
        proj2_full, ...
        U, ...
        beta);

    % --------------------------------------------------------
    % Output
    % --------------------------------------------------------
    %
    % Columns:
    %
    % 1  Frobenius error
    % 2  Spectral error
    % 3  Mean absolute error
    % 4  TPR
    % 5  FPR
    % 6  MCC
    % 7  PD
    % 8  Achieved matrix depth
    % 9  Final optimization time
    %

    output = ...
        [ErrFro, ...
         ErrSpe, ...
         err_abs, ...
         TPR, ...
         FPR, ...
         MCC, ...
         PD, ...
         Depth, ...
         comp_time_final];

    filename = sprintf( ...
        ['outputs/' ...
         'output_FroSpecAbsTprFprMccPDDepthTime_' ...
         'method_%d_p_%d_n_%d_contam_%02d_K_%d_' ...
         'BLOC_penalty_%d_REP_%d.csv'], ...
         method_num, ...
         p, ...
         n, ...
         round(100*contam_prop), ...
         num_U, ...
         penalty, ...
         k);

    csvwrite(filename,output);
    output

end