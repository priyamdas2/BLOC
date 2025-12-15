clearvars;
warning('off', 'all');
addpath('./BLOC/');
addpath('./Simulation data/');
reps_to_do = [6,7,8]; % Perform any rep <subset of> {1,2,3,4,5,6,7,8,9,10}

rng(1)
% ------------------- User Inputs -------------------
p = 100;
n = 150;
threshold = 1e-2; % don't change
Ctype = 'Blockdiag';   %  'Blockdiag','SparseUniform'
methods_to_perform = [5];  % 1 = LASSO, 2 = Ridge, 3 = ENet, 4 = SCAD, 5 = MCP

% ------------------- BLOC Inputs -------------------
doparallel = 1; % set 0 if doing parfor over data reps
%maxNumCompThreads(1);
MaxTime = 2000;
MaxRuns = 10;
MaxIter = 500;
TolFun1 = 1e-4;
TolFun2 = 1e-2;
phi = 1e-6;
DisplayUpdate = 1;
DisplayEvery = 10;
% ------------------- Lambda grid ---------------------
log10lambdaLB = -2;
log10lambdaUB = 1;
P = ones(p,p) - eye(p);
NumLambdas = 20;
lambda_array = 10.^linspace(log10lambdaLB, log10lambdaUB, NumLambdas);
lambda_array = lambda_array(9:13);

a_array_SCAD = linspace(3, 5, 5); 
gamma_array_MCP = linspace(2, 4, 5); 

% ------------ Read correlation matrix -----------------
C_true_filename = sprintf('C_p_%d_n_%d_C_%s.csv', p, n, Ctype);
C_true = readmatrix(C_true_filename);
% ------------------------------------------------------



for data_rep = reps_to_do
 
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
    
    
    
    
   
    % -------- SCAD -------------------------------------------------------
    
    if ismember(4, methods_to_perform)
        tStart = tic;
        [CorrEst, SigmaEst] = BLOC_SpCovSCAD(X_raw, P, lambda_array, a_array_SCAD, doparallel, ...
            'MaxTime', MaxTime, 'MaxRuns', MaxRuns, 'MaxIter', MaxIter,...
            'TolFun1', TolFun1, 'TolFun2', TolFun2, 'phi', phi, 'DisplayUpdate', DisplayUpdate, 'DisplayEvery', DisplayEvery);
        
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
        ALL_SCAD = [TPR, FPR, MCC, RMSE, MAD, time_taken, lambda_array(best_i), a_array_SCAD(best_j)];

        filename = sprintf('Simulation output/TPR_FPR_MCC_RMSE_MAD_time_BestParams_BLOC_SCAD_p_%d_n_%d_C_%s_RepNum_%d.csv', ...
            p, n, Ctype, data_rep);
        writematrix(ALL_SCAD, filename);
        
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
        ALL_MCP = [TPR, FPR, MCC, RMSE, MAD, time_taken, lambda_array(best_i), gamma_array_MCP(best_j)];
        filename = sprintf('Simulation output/TPR_FPR_MCC_RMSE_MAD_time_BestParams_BLOC_MCP_p_%d_n_%d_C_%s_RepNum_%d.csv', ...
            p, n, Ctype, data_rep);
        writematrix(ALL_MCP, filename);
    end
   
end





