clear all; clc;
clearvars;

addpath('./supp funs/');
addpath(genpath('./manopt/'));
% ============================================================
% Simulation design
% MUST MATCH BLOC
% ============================================================

p  = 20;
method_num = 2; % 1 = Block, 2 = Toeplitz, 3 = Banded

SIGMA = cov_model_matDepthsim(p, method_num);

penalty = 2; % 1 = SCAD, 2 = MCP

n = 200;

reps_to_do = 1:10; % any subset of {1,2,...,10}

tol = 1e-5;

% ============================================================
% Gross contamination
% MUST MATCH BLOC
% ============================================================

contam_prop = 0.10;

num_outliers = round(contam_prop*n);

outlier_scale = 3;

% ============================================================
% Fixed random directions
% MUST MATCH BLOC EXACTLY
% ============================================================

num_U = 500;

U_seed = 246810;

rng(U_seed,'twister');

U = randn(p,num_U);

U = U ./ sqrt(sum(U.^2,1));

% Gaussian conversion factor
beta = norminv(0.75)^2;

% ============================================================
% Cross-validation settings
% MUST MATCH BLOC
% ============================================================

n2 = round(n/log(n)); % Test
n1 = n - n2;          % Train

num_CV = 5;

num_lambda = 15;

lamdas = logspace(-3,0,num_lambda);

C0 = eye(p);

% ============================================================
% Manopt time limits
%
% fmincon is intentionally kept in the same form as the
% previous benchmark code.
% ============================================================

maxtime_CV = 200;

maxtime_final = 3600;

% ============================================================
% Reproducibility
% ============================================================

CV_seed_base = 100000;

% Only relevant for Manopt/random competitor behavior.
% Data/CV seeds are kept completely separate.
COMP_CV_seed_base = 400000;

COMP_Final_seed_base = 500000;

% ============================================================
% Competitors
%
% ROW ORDER IN OUTPUT:
%
% 1 = fmincon active-set
% 2 = fmincon interior-point
% 3 = fmincon SQP
% 4 = Manopt Barzilai-Borwein
% 5 = Manopt conjugate-gradient
% 6 = Manopt steepest-descent
% 7 = Manopt trust-region
% ============================================================

num_methods = 7;

method_names = { ...
    'fmincon_active_set'
    'fmincon_interior_point'
    'fmincon_SQP'
    'Manopt_BB'
    'Manopt_CG'
    'Manopt_SD'
    'Manopt_TR'
    };

% ============================================================
% Skip flags
%
% 0 = run
% 1 = skip
% ============================================================

skip_activeset      = 0;
skip_interiorpoint  = 0;
skip_sqp            = 0;

skip_BB             = 0;
skip_CG             = 0;
skip_SD             = 0;
skip_TR             = 1;

skip_method = [ ...
    skip_activeset
    skip_interiorpoint
    skip_sqp
    skip_BB
    skip_CG
    skip_SD
    skip_TR
    ];

% ============================================================

if ~exist('outputs','dir')
    mkdir('outputs');
end

% ============================================================
% Simulation replications
% ============================================================

parfor k = reps_to_do

    fprintf('\n');
    fprintf('====================================================\n');
    fprintf('       COMPETITOR SIMULATION -- REP %d\n',k);
    fprintf('====================================================\n');

    % ========================================================
    % DATA GENERATION
    %
    % KEEP IDENTICAL TO BLOC
    % ========================================================

    rng(k,'twister');

    % Clean Gaussian sample
    r = mvnrnd(zeros(p,1),SIGMA,n);

    % Select exactly 10% for contamination
    outlier_idx = randperm(n,num_outliers);

    % Replace with N(0, 9 I)
    r(outlier_idx,:) = ...
        mvnrnd( ...
        zeros(p,1), ...
        outlier_scale^2*eye(p), ...
        num_outliers);

    % Squared projections for full-data fit
    proj2_full = (r*U).^2;

    % ========================================================
    % CV SPLITS
    %
    % KEEP IDENTICAL TO BLOC
    % ========================================================

    rng(CV_seed_base+k,'twister');

    CV_indices = zeros(num_CV,n);

    for kk = 1:num_CV

        CV_indices(kk,:) = randperm(n);

    end

    % ========================================================
    % Store CV errors
    %
    % Dimensions:
    %
    % method x CV fold x lambda
    % ========================================================

    FroErr = nan(num_methods,num_CV,num_lambda);

    % ========================================================
    % CROSS VALIDATION
    % ========================================================

    for kk = 1:num_CV

        fprintf('\n');
        fprintf('REP %d -- CV fold %d/%d\n', ...
            k,kk,num_CV);

        J = CV_indices(kk,:);

        r1 = r(J(1:n1),:);

        r2 = r(J(n1+1:end),:);

        % ----------------------------------------------------
        % Training projections
        % ----------------------------------------------------

        proj2_train = (r1*U).^2;

        % ----------------------------------------------------
        % Robust held-out correlation target
        %
        % SAME criterion as BLOC
        % ----------------------------------------------------

        Tau2 = corr(r2,'Type','Kendall');

        S2_robust = sin((pi/2)*Tau2);

        S2_robust(1:p+1:end) = 1;

        % ====================================================
        % Lambda grid
        % ====================================================

        for l = 1:num_lambda

            lambda = lamdas(l);

            fprintf('\n');
            fprintf(['REP %d | CV %d | lambda %d/%d ' ...
                     '= %.6g\n'], ...
                     k,kk,l,num_lambda,lambda);

            % ------------------------------------------------
            % Define SAME objective for every optimizer
            % ------------------------------------------------

            if penalty == 1

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

            % =================================================
            % Run all competitor optimizers
            % =================================================

            for m = 1:num_methods

                if skip_method(m) == 1
                    continue;
                end

                fprintf('   -> %s\n',method_names{m});

                rand_seed = ...
                    COMP_CV_seed_base + ...
                    100000*m + ...
                    10000*k + ...
                    100*kk + ...
                    l;

                try

                    [C_opt,~,~] = ...
                        fit_competitor_matDepth( ...
                        m, ...
                        objFun, ...
                        C0, ...
                        p, ...
                        maxtime_CV, ...
                        rand_seed);

                    % -----------------------------------------
                    % Same robust held-out CV criterion
                    % as BLOC
                    % -----------------------------------------

                    FroErr(m,kk,l) = ...
                        norm( ...
                        C_opt-S2_robust, ...
                        'fro');

                catch ME

                    warning( ...
                        ['REP %d, CV %d, lambda %.6g, ' ...
                         '%s failed:\n%s'], ...
                        k, ...
                        kk, ...
                        lambda, ...
                        method_names{m}, ...
                        ME.message);

                    FroErr(m,kk,l) = NaN;

                end

            end

        end

    end

    % ========================================================
    % Select lambda SEPARATELY for each competitor
    % ========================================================

    lambda_opt = nan(num_methods,1);

    CV_min = nan(num_methods,1);

    fprintf('\n');
    fprintf('====================================================\n');
    fprintf('        SELECTED LAMBDAS -- REP %d\n',k);
    fprintf('====================================================\n');

    for m = 1:num_methods

        if skip_method(m) == 1
            continue;
        end

        tempErr = squeeze(FroErr(m,:,:));

        % tempErr:
        % num_CV x num_lambda

        meanErr = mean(tempErr,1,'omitnan');

        % Require all CV folds to have completed successfully
        valid_lambda = ...
            all(isfinite(tempErr),1);

        meanErr(~valid_lambda) = Inf;

        [bestErr,mi] = min(meanErr);

        if isfinite(bestErr)

            lambda_opt(m) = lamdas(mi);

            CV_min(m) = bestErr;

            fprintf('\n%s\n',method_names{m});

            fprintf( ...
                'Selected lambda = %.6f\n', ...
                lambda_opt(m));

            fprintf( ...
                'Mean robust CV error = %.6f\n', ...
                bestErr);

            disp('      lambda      mean robust CV error');

            disp([lamdas(:),meanErr(:)]);

        else

            warning( ...
                '%s has no valid lambda.', ...
                method_names{m});

        end

    end

    % ========================================================
    % FINAL FITS
    %
    % output_competitors:
    %
    % row = optimizer
    %
    % columns:
    %
    % 1 = Frobenius error
    % 2 = Spectral error
    % 3 = Mean absolute error
    % 4 = TPR
    % 5 = FPR
    % 6 = MCC
    % 7 = PD
    % 8 = Average directional depth
    % 9 = Final optimization time
    % ========================================================

    output_competitors = nan(num_methods,9);

    % --------------------------------------------------------
    % ONE FILE PER REP containing ALL competitors
    % --------------------------------------------------------

    filename = sprintf( ...
        ['outputs/' ...
         'output_FroSpecAbsTprFprMccPDDepthTime_' ...
         'method_%d_p_%d_n_%d_contam_%02d_K_%d_' ...
         'COMPETITORS_penalty_%d_REP_%d.csv'], ...
         method_num, ...
         p, ...
         n, ...
         round(100*contam_prop), ...
         num_U, ...
         penalty, ...
         k);

    % ========================================================
    % Final optimization
    % ========================================================

    for m = 1:num_methods

        if skip_method(m) == 1
            continue;
        end

        if ~isfinite(lambda_opt(m))
            continue;
        end

        fprintf('\n');
        fprintf('====================================================\n');
        fprintf('FINAL FIT\n');
        fprintf('REP    : %d\n',k);
        fprintf('METHOD : %s\n',method_names{m});
        fprintf('LAMBDA : %.6f\n',lambda_opt(m));
        fprintf('====================================================\n');

        % ----------------------------------------------------
        % Full-data objective
        % ----------------------------------------------------

        if penalty == 1

            a = 3.7;

            objFun = @(Theta) ...
                obj_depth_scad_corr_avg( ...
                Theta, ...
                proj2_full, ...
                U, ...
                beta, ...
                lambda_opt(m), ...
                a);

        else

            gamma = 3;

            objFun = @(Theta) ...
                obj_depth_mcp_corr_avg( ...
                Theta, ...
                proj2_full, ...
                U, ...
                beta, ...
                lambda_opt(m), ...
                gamma);

        end

        rand_seed = ...
            COMP_Final_seed_base + ...
            100000*m + ...
            k;

        try

            [C_final,fval_final,comp_time_final] = ...
                fit_competitor_matDepth( ...
                m, ...
                objFun, ...
                C0, ...
                p, ...
                maxtime_final, ...
                rand_seed);

            % =================================================
            % Performance measures
            % =================================================

            ErrFro = ...
                norm(C_final-SIGMA,'fro');

            ErrSpe = ...
                norm(C_final-SIGMA);

            err_abs = ...
                mean( ...
                abs(C_final(:)-SIGMA(:)));

            [TPR,FPR,MCC] = ...
                TPR_FPR_MCC( ...
                C_final, ...
                SIGMA, ...
                tol);

            PD = ...
                min(real(eig(C_final))) > 1e-8;

            Depth = ...
                matrix_depth_fixedU_avg( ...
                C_final, ...
                proj2_full, ...
                U, ...
                beta);

            output_competitors(m,:) = ...
                [ ...
                ErrFro, ...
                ErrSpe, ...
                err_abs, ...
                TPR, ...
                FPR, ...
                MCC, ...
                PD, ...
                Depth, ...
                comp_time_final ...
                ];

            fprintf('\n');
            fprintf( ...
                'Final objective = %.8f\n', ...
                fval_final);

            fprintf( ...
                ['Fro %.4f | Spec %.4f | Abs %.4f | ' ...
                 'TPR %.4f | FPR %.4f | MCC %.4f | ' ...
                 'PD %d | AvgDepth %.4f | Time %.4f\n'], ...
                 output_competitors(m,:));

        catch ME

            warning( ...
                'FINAL %s failed in REP %d:\n%s', ...
                method_names{m}, ...
                k, ...
                ME.message);

        end

        % ----------------------------------------------------
        % SAVE AFTER EACH COMPETITOR
        %
        % This means that if competitor 6 crashes, rows 1--5
        % have already been written to the REP file.
        % ----------------------------------------------------

        writematrix( ...
            output_competitors, ...
            filename);

    end

    % ========================================================
    % Final write
    % ========================================================

    writematrix( ...
        output_competitors, ...
        filename);

    fprintf('\n');
    fprintf('Saved REP %d output:\n%s\n', ...
        k,filename);

    fprintf('\nROW ORDER:\n');

    for m = 1:num_methods

        fprintf('%d = %s\n', ...
            m,method_names{m});

    end

    fprintf('\nFINAL OUTPUT MATRIX:\n');

    disp(output_competitors);

end