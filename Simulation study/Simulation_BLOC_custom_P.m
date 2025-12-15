clearvars;
warning('off', 'all');
addpath('./BLOC/');
addpath('./Simulation data/');
num_dataset_rep = 10;

rng(1)
% ------------------- User Inputs -------------------
p = 50;
n = 200;
threshold = 1e-2;
Ctype = 'SparseUniform';   %  'Blockdiag','SparseUniform'
methods_to_perform = [4];  % 1 = LASSO, 2 = Ridge, 3 = ENet, 4 = SCAD, 5 = MCP

% ------------------- BLOC Inputs -------------------
doparallel = 0; % set 0 if doing parfor over data reps
%maxNumCompThreads(1);
MaxTime = 3600;
MaxRuns = 5;
MaxIter = 1000;
TolFun1 = 1e-4;
TolFun2 = 1e-2;
phi = 1e-6;
DisplayUpdate = 1;
DisplayEvery = 2;
% ------------------- Lambda grid ---------------------
log10lambdaLB = -2;
log10lambdaUB = 1;
P_filename = sprintf('P_cover_p_%d_n_%d_%s.csv', p, n, Ctype);
P = readmatrix(P_filename);
NumLambdas = 20;
lambda_array = 10.^linspace(log10lambdaLB, log10lambdaUB, NumLambdas);

a_array_SCAD = linspace(3, 5, 5); % 3.7
gamma_array_MCP = linspace(2, 4, 5); %3

% ------------ Read correlation matrix -----------------
C_true_filename = sprintf('C_p_%d_n_%d_C_%s.csv', p, n, Ctype);
C_true = readmatrix(C_true_filename);
% ------------------------------------------------------

ALL_LASSO  = nan(num_dataset_rep, 6);
ALL_Ridge  = nan(num_dataset_rep, 6);
ALL_ENet   = nan(num_dataset_rep, 6);
ALL_SCAD   = nan(num_dataset_rep, 6);
ALL_MCP    = nan(num_dataset_rep, 6);
Best_lambda_LASSO = nan(num_dataset_rep, 1);
Best_lambda_Ridge = nan(num_dataset_rep, 1);
Best_lambda_ENet = nan(num_dataset_rep, 1);
Best_lambda_SCAD = nan(num_dataset_rep, 1);
Best_lambda_MCP = nan(num_dataset_rep, 1);

parfor data_rep = 1:num_dataset_rep
    fprintf('Performing data rep number: %d\n', data_rep);
    folder = 'Simulation data';
    filename = sprintf('X_p_%d_n_%d_C_%s_DataRep_%d.csv', p, n, Ctype,data_rep);
    filepath = fullfile(folder, filename);
    if isfile(filepath)
        X_raw = readmatrix(filepath);
        fprintf('Loaded data from: %s\n', filepath);
    else
        error('File not found: %s', filepath);
    end
    
    
    
    
    % -------- LASSO ------------------------------------------------------
    
    if ismember(1, methods_to_perform)
        tStart = tic;
        
        [CorrEst, SigmaEst] = BLOC_SpCovLasso(X_raw, P, lambda_array, doparallel, ...
            'MaxTime', MaxTime, 'MaxRuns', MaxRuns, 'MaxIter', MaxIter,...
            'TolFun1', TolFun1, 'TolFun2', TolFun2, 'phi', phi, 'DisplayUpdate', DisplayUpdate);
        
        EBICs = zeros(1, NumLambdas);
        for ii = 1:NumLambdas
            EBICs(ii) = compute_EBIC_cov(CorrEst(:,:,ii), X_raw, threshold);
        end
        [~, best_idx] = min(EBICs);
        C_est = CorrEst(:,:,best_idx);
        [TPR, FPR, MCC, RMSE, MAD] = evaluate_corr_metrics(C_est, C_true, threshold);
        time_taken = toc(tStart);
        ALL_LASSO(data_rep, :) = [TPR, FPR, MCC, RMSE, MAD, time_taken]';
        Best_lambda_LASSO(data_rep) = lambda_array(best_idx);
    end
    
    % ---------------------------------------------------------------------
    
    % -------- RIDGE ------------------------------------------------------
    
    if ismember(2, methods_to_perform)
        tStart = tic;
        [CorrEst, SigmaEst] = BLOC_SpCovRidge(X_raw, P, lambda_array, doparallel, ...
            'MaxTime', MaxTime, 'MaxRuns', MaxRuns, 'MaxIter', MaxIter,...
            'TolFun1', TolFun1, 'TolFun2', TolFun2, 'phi', phi, 'DisplayUpdate', DisplayUpdate);
        
        EBICs = nan(1, NumLambdas);
        for ii = 1:NumLambdas
            EBICs(ii) = compute_EBIC_cov(CorrEst(:,:,ii), X_raw, threshold);
        end
        [~, best_idx] = min(EBICs);
        C_est = CorrEst(:,:,best_idx);
        [TPR, FPR, MCC, RMSE, MAD] = evaluate_corr_metrics(C_est, C_true, threshold);
        time_taken = toc(tStart);
        ALL_Ridge(data_rep, :) = [TPR, FPR, MCC, RMSE, MAD, time_taken]';
        Best_lambda_Ridge(data_rep) = lambda_array(best_idx);
    end
    
    % -------- ENet -------------------------------------------------------
    
    if ismember(3, methods_to_perform)
        tStart = tic;
        alpha_array = linspace(0,1,5);
        [CorrEst, SigmaEst] = BLOC_SpCovElasticNet(X_raw, P, lambda_array, alpha_array, doparallel, ...
            'MaxTime', MaxTime, 'MaxRuns', MaxRuns, 'MaxIter', MaxIter,...
            'TolFun1', TolFun1, 'TolFun2', TolFun2, 'phi', phi, 'DisplayUpdate', DisplayUpdate);
        
        NumAlphas = length(alpha_array);
        EBICs = nan(NumLambdas, NumAlphas);
        for ii = 1:NumLambdas
            for jj = 1:NumAlphas
                EBICs(ii,jj) = compute_EBIC_cov(CorrEst(:,:,ii,jj), X_raw, threshold);
            end
        end
        [min_val, linear_idx] = min(EBICs(:));  % Get linear index of min
        [best_i, best_j] = ind2sub(size(EBICs), linear_idx);
        
        C_est = CorrEst(:,:,best_i, best_j);
        [TPR, FPR, MCC, RMSE, MAD] = evaluate_corr_metrics(C_est, C_true, threshold);
        time_taken = toc(tStart);
        ALL_ENet(data_rep, :) = [TPR, FPR, MCC, RMSE, MAD, time_taken]';
        Best_lambda_ENet(data_rep) = lambda_array(best_i);
    end
    
    % -------- SCAD -------------------------------------------------------
    
    if ismember(4, methods_to_perform)
        tStart = tic;
        [CorrEst, SigmaEst] = BLOC_SpCovSCAD(X_raw, P, lambda_array, a_array_SCAD, doparallel, ...
            'MaxTime', MaxTime, 'MaxRuns', MaxRuns, 'MaxIter', MaxIter,...
            'TolFun1', TolFun1, 'TolFun2', TolFun2, 'phi', phi, 'DisplayUpdate', DisplayUpdate);
        
        NumAs = length(a_array_SCAD);
        EBICs = nan(NumLambdas, NumAs);
        for ii = 1:NumLambdas
            for jj = 1:NumAs
                EBICs(ii,jj) = compute_EBIC_cov(CorrEst(:,:,ii,jj), X_raw, threshold);
            end
        end
        [min_val, linear_idx] = min(EBICs(:));  % Get linear index of min
        [best_i, best_j] = ind2sub(size(EBICs), linear_idx);
        
        C_est = CorrEst(:,:,best_i, best_j);
        [TPR, FPR, MCC, RMSE, MAD] = evaluate_corr_metrics(C_est, C_true, threshold);
        time_taken = toc(tStart);
        ALL_SCAD(data_rep, :) = [TPR, FPR, MCC, RMSE, MAD, time_taken]';
        Best_lambda_SCAD(data_rep) = lambda_array(best_i);
    end
    % -------- MC+ --------------------------------------------------------
    
    if ismember(5, methods_to_perform)
        tStart = tic;
        
        [CorrEst, SigmaEst] = BLOC_SpCovMCplus(X_raw, P, lambda_array, gamma_array_MCP, doparallel, ...
            'MaxTime', MaxTime, 'MaxRuns', MaxRuns, 'MaxIter', MaxIter,...
            'TolFun1', TolFun1, 'TolFun2', TolFun2, 'phi', phi, 'DisplayUpdate', DisplayUpdate);
        
        NumGammas = length(gamma_array_MCP);
        EBICs = nan(NumLambdas, NumGammas);
        for ii = 1:NumLambdas
            for jj = 1:NumGammas
                EBICs(ii,jj) = compute_EBIC_cov(CorrEst(:,:,ii,jj), X_raw, threshold);
            end
        end
        [min_val, linear_idx] = min(EBICs(:));  % Get linear index of min
        [best_i, best_j] = ind2sub(size(EBICs), linear_idx);
        
        C_est = CorrEst(:,:,best_i, best_j);
        [TPR, FPR, MCC, RMSE, MAD] = evaluate_corr_metrics(C_est, C_true, threshold);
        time_taken = toc(tStart);
        ALL_MCP(data_rep, :) = [TPR, FPR, MCC, RMSE, MAD, time_taken]';
        Best_lambda_MCP(data_rep) = lambda_array(best_i);
    end
    
    % ---------------------------------------------------------------------
end

metric_names = {'TPR', 'FPR', 'MCC', 'RMSE', 'MAD', 'Comp_time'};
output_folder = 'Simulation output';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Save each method matrix
method_names = {'LASSO', 'Ridge', 'ENet', 'SCAD', 'MCP'};
all_results = {ALL_LASSO, ALL_Ridge, ALL_ENet, ALL_SCAD, ALL_MCP};



for i = 1:length(method_names)
    if ismember(i, methods_to_perform)
        fprintf('Best lambdas: \n');
        switch i
            case 1
                disp(Best_lambda_LASSO');
            case 2
                disp(Best_lambda_Ridge');
            case 3
                disp(Best_lambda_ENet');
            case 4
                disp(Best_lambda_SCAD');
            case 5
                disp(Best_lambda_MCP');
            otherwise
                disp('Invalid value of i.');
        end
        fprintf('Method performed: %s\n', method_names{i});
        fprintf('TPR, FPR, MCC, RMSE, MAD, Comp_time: %.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n', mean(all_results{i}));

        T = array2table(all_results{i}, 'VariableNames', metric_names);
        filename = sprintf('TPR_FPR_MCC_RMSE_MAD_time_custom_P_p_%d_n_%d_C_%s_NumReps_%d_BLOC_%s.csv', ...
            p, n, Ctype, num_dataset_rep, method_names{i});
        writetable(T, fullfile(output_folder, filename), 'WriteRowNames', false);
    end
end
