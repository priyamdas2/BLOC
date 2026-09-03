clearvars;
warning('off', 'all');
% Execute the line below if running on VCU HPC cluster
% run('/lustre/home/dasp4/MATLAB_cvx_linux/cvx/cvx_startup.m')

addpath('./Simulation data/');
addpath('./supp funs/');
addpath('./LLA_MM supp funs/');

num_dataset_rep = 10;

rng(1);

% ================================================================
% User Inputs
% ================================================================

p = 20;
n = 50;

threshold = 1e-2;        % DON'T CHANGE

Ctype = 'SparseUniform';     % 'Blockdiag','SparseUniform'

% Keep original BLOC numbering:
% 4 = SCAD
% 5 = MCP
methods_to_perform = [4];


% ================================================================
% LLA/MM Inputs
% ================================================================

opts.maxOuter = 20;
opts.maxInit  = 10;

opts.tol = 1e-5;

% Match the numerical PD threshold used in the BLOC
% Gaussian likelihood functions.
opts.eig_floor = 1e-4;

% L1/MM contraction initialization before SCAD/MCP tightening.
opts.doL1Init = true;

% Set false for large simulation runs if output becomes excessive.
opts.verbose = false;


% ================================================================
% Lambda / shape-parameter grids
%
% EXACTLY SAME AS BLOC
% ================================================================

log10lambdaLB = -2;
log10lambdaUB = 1;

P = ones(p,p) - eye(p);

NumLambdas = 20;
if p > 50
    NumLambdas = 10;
end

lambda_array = ...
    10.^linspace( ...
    log10lambdaLB, ...
    log10lambdaUB, ...
    NumLambdas);


% Same grids used in BLOC
a_array_SCAD = 3.7;

gamma_array_MCP = 3;


% ================================================================
% Read true correlation matrix
% ================================================================

C_true_filename = ...
    sprintf( ...
    'C_p_%d_n_%d_C_%s.csv', ...
    p,n,Ctype);

C_true = ...
    readmatrix(C_true_filename);


% ================================================================
% Storage
%
% Columns:
%
% TPR, FPR, MCC, RMSE, MAD, Comp_time
% ================================================================

ALL_SCAD = ...
    nan(num_dataset_rep,6);

ALL_MCP = ...
    nan(num_dataset_rep,6);


Best_lambda_SCAD = ...
    nan(num_dataset_rep,1);

Best_lambda_MCP = ...
    nan(num_dataset_rep,1);


% Also save the selected SCAD/MCP shape parameter.
% BLOC did not previously save these, but it is useful for
% documenting this matched comparison.
Best_a_SCAD = ...
    nan(num_dataset_rep,1);

Best_gamma_MCP = ...
    nan(num_dataset_rep,1);


% ================================================================
% Data replications
%
% NOTE:
% Use FOR rather than PARFOR because each LLA/MM fit calls CVX.
% Data, tuning and metric calculations remain identical to BLOC.
% ================================================================

for data_rep = 1:num_dataset_rep
    
    fprintf('\n');
    fprintf('====================================================\n');
    fprintf('Performing data rep number: %d\n',data_rep);
    fprintf('====================================================\n');
    
    
    % ============================================================
    % Read EXACT SAME simulated dataset used by BLOC
    % ============================================================
    
    folder = 'Simulation data';
    
    filename = ...
        sprintf( ...
        'X_p_%d_n_%d_C_%s_DataRep_%d.csv', ...
        p,n,Ctype,data_rep);
    
    filepath = ...
        fullfile(folder,filename);
    
    
    if isfile(filepath)
        
        X_raw = ...
            readmatrix(filepath);
        
        fprintf( ...
            'Loaded data from: %s\n', ...
            filepath);
        
    else
        
        error( ...
            'File not found: %s', ...
            filepath);
        
    end
    
    
    % ============================================================
    % SCAD
    % ============================================================
    
    if ismember(4,methods_to_perform)
        
        fprintf('\n');
        fprintf('---------------- SCAD LLA/MM ----------------\n');
        
        tStart = tic;
        
        
        NumAs = ...
            length(a_array_SCAD);
        
        
        % --------------------------------------------------------
        % Store estimated correlation matrix for every
        % lambda x a combination.
        %
        % Same dimensional organization as BLOC:
        %
        % p x p x NumLambdas x NumAs
        % --------------------------------------------------------
        
        CorrEst = ...
            nan( ...
            p,p, ...
            NumLambdas, ...
            NumAs);
        
        
        % --------------------------------------------------------
        % Fit entire tuning grid
        % --------------------------------------------------------
        
        for ii = 1:NumLambdas
            
            lambda = lambda_array(ii);
            
            for jj = 1:NumAs
                
                a = a_array_SCAD(jj);
                
                fprintf( ...
                    ['Data rep %d | SCAD | p = %d, n = %d | ' ...
                    'lambda %d/%d = %.6g | a = %.3f\n'], ...
                    data_rep, p, n, ...
                    ii, NumLambdas, lambda, a);
                
                try
                    
                    C0 = eye(p);
                    
                    % -----------------------------------------
                    % Time this individual lambda/a fit
                    % -----------------------------------------
                    tLambda = tic;
                    
                    [C_tmp,info_tmp] = ...
                        LLA_MM_GaussCorr( ...
                        X_raw, ...
                        lambda, ...
                        'SCAD', ...
                        a, ...
                        P, ...
                        C0, ...
                        opts); %#ok<NASGU>
                    
                    lambda_time = toc(tLambda);
                    
                    CorrEst(:,:,ii,jj) = C_tmp;
                    
                    fprintf( ...
                        ['>>> p = %d, n = %d | SCAD | ' ...
                        'lambda %d/%d | time taken = %.2f sec ' ...
                        '(%.2f min)\n'], ...
                        p, n, ii, NumLambdas, ...
                        lambda_time, lambda_time/60);
                    
                catch ME
                    
                    lambda_time = toc(tLambda);
                    
                    fprintf( ...
                        ['WARNING: SCAD LLA/MM failed | ' ...
                        'p = %d, n = %d | lambda %d/%d | ' ...
                        'time = %.2f sec\n'], ...
                        p, n, ii, NumLambdas, lambda_time);
                    
                    fprintf('Reason: %s\n', ME.message);
                    
                    CorrEst(:,:,ii,jj) = nan(p);
                    
                end
                
            end
        end
        % --------------------------------------------------------
        % EBIC tuning
        %
        % EXACT SAME EBIC function as BLOC.
        % --------------------------------------------------------
        
        EBICs = ...
            inf( ...
            NumLambdas, ...
            NumAs);
        
        
        for ii = 1:NumLambdas
            
            for jj = 1:NumAs
                
                C_candidate = ...
                    CorrEst(:,:,ii,jj);
                
                
                if all(isfinite(C_candidate(:)))
                    
                    EBICs(ii,jj) = ...
                        compute_EBIC_cov( ...
                        C_candidate, ...
                        X_raw, ...
                        threshold);
                    
                end
                
            end
            
        end
        
        
        % --------------------------------------------------------
        % Select best lambda and a
        % --------------------------------------------------------
        
        [min_val,linear_idx] = ...
            min(EBICs(:)); %#ok<ASGLU>
        
        
        if ~isfinite(min_val)
            
            error( ...
                ['All SCAD LLA/MM tuning combinations failed ' ...
                'for data replication %d.'], ...
                data_rep);
            
        end
        
        
        [best_i,best_j] = ...
            ind2sub( ...
            size(EBICs), ...
            linear_idx);
        
        
        C_est = ...
            CorrEst(:,:,best_i,best_j);
        
        
        % --------------------------------------------------------
        % SAME performance metrics as BLOC
        % --------------------------------------------------------
        
        [TPR,FPR,MCC,RMSE,MAD] = ...
            evaluate_corr_metrics( ...
            C_est, ...
            C_true, ...
            threshold);
        
        
        time_taken = ...
            toc(tStart);
        
        
        ALL_SCAD(data_rep,:) = ...
            [ ...
            TPR, ...
            FPR, ...
            MCC, ...
            RMSE, ...
            MAD, ...
            time_taken ...
            ];
        
        
        Best_lambda_SCAD(data_rep) = ...
            lambda_array(best_i);
        
        
        Best_a_SCAD(data_rep) = ...
            a_array_SCAD(best_j);
        
        
        fprintf('\n');
        
        fprintf( ...
            ['SCAD rep %d selected: ' ...
            'lambda = %.6g, a = %.3f\n'], ...
            data_rep, ...
            Best_lambda_SCAD(data_rep), ...
            Best_a_SCAD(data_rep));
        
        fprintf( ...
            ['TPR = %.4f | FPR = %.4f | MCC = %.4f | ' ...
            'RMSE = %.4f | MAD = %.4f | Time = %.2f\n'], ...
            TPR,FPR,MCC,RMSE,MAD,time_taken);
        
    end
    
    
    % ============================================================
    % MCP / MC+
    % ============================================================
    
    if ismember(5,methods_to_perform)
        
        fprintf('\n');
        fprintf('---------------- MCP LLA/MM -----------------\n');
        
        tStart = tic;
        
        
        NumGammas = ...
            length(gamma_array_MCP);
        
        
        % --------------------------------------------------------
        % p x p x NumLambdas x NumGammas
        % --------------------------------------------------------
        
        CorrEst = ...
            nan( ...
            p,p, ...
            NumLambdas, ...
            NumGammas);
        
        
        % --------------------------------------------------------
        % Fit entire tuning grid
        % --------------------------------------------------------
        
        for ii = 1:NumLambdas
            
            lambda = lambda_array(ii);
            
            for jj = 1:NumGammas
                
                gamma = gamma_array_MCP(jj);
                
                fprintf( ...
                    ['Data rep %d | MCP | p = %d, n = %d | ' ...
                    'lambda %d/%d = %.6g | gamma = %.3f\n'], ...
                    data_rep, p, n, ...
                    ii, NumLambdas, lambda, gamma);
                
                try
                    
                    C0 = eye(p);
                    
                    % -----------------------------------------
                    % Time this individual lambda/gamma fit
                    % -----------------------------------------
                    tLambda = tic;
                    
                    [C_tmp,info_tmp] = ...
                        LLA_MM_GaussCorr( ...
                        X_raw, ...
                        lambda, ...
                        'MCP', ...
                        gamma, ...
                        P, ...
                        C0, ...
                        opts); %#ok<NASGU>
                    
                    lambda_time = toc(tLambda);
                    
                    CorrEst(:,:,ii,jj) = C_tmp;
                    
                    fprintf( ...
                        ['>>> p = %d, n = %d | MCP | ' ...
                        'lambda %d/%d | time taken = %.2f sec ' ...
                        '(%.2f min)\n'], ...
                        p, n, ii, NumLambdas, ...
                        lambda_time, lambda_time/60);
                    
                catch ME
                    
                    lambda_time = toc(tLambda);
                    
                    fprintf( ...
                        ['WARNING: MCP LLA/MM failed | ' ...
                        'p = %d, n = %d | lambda %d/%d | ' ...
                        'time = %.2f sec\n'], ...
                        p, n, ii, NumLambdas, lambda_time);
                    
                    fprintf('Reason: %s\n', ME.message);
                    
                    CorrEst(:,:,ii,jj) = nan(p);
                    
                end
                
            end
        end
        
        
        % --------------------------------------------------------
        % EBIC tuning
        % --------------------------------------------------------
        
        EBICs = ...
            inf( ...
            NumLambdas, ...
            NumGammas);
        
        
        for ii = 1:NumLambdas
            
            for jj = 1:NumGammas
                
                C_candidate = ...
                    CorrEst(:,:,ii,jj);
                
                
                if all(isfinite(C_candidate(:)))
                    
                    EBICs(ii,jj) = ...
                        compute_EBIC_cov( ...
                        C_candidate, ...
                        X_raw, ...
                        threshold);
                    
                end
                
            end
            
        end
        
        
        % --------------------------------------------------------
        % Select best lambda and gamma
        % --------------------------------------------------------
        
        [min_val,linear_idx] = ...
            min(EBICs(:)); %#ok<ASGLU>
        
        
        if ~isfinite(min_val)
            
            error( ...
                ['All MCP LLA/MM tuning combinations failed ' ...
                'for data replication %d.'], ...
                data_rep);
            
        end
        
        
        [best_i,best_j] = ...
            ind2sub( ...
            size(EBICs), ...
            linear_idx);
        
        
        C_est = ...
            CorrEst(:,:,best_i,best_j);
        
        
        % --------------------------------------------------------
        % SAME performance metrics as BLOC
        % --------------------------------------------------------
        
        [TPR,FPR,MCC,RMSE,MAD] = ...
            evaluate_corr_metrics( ...
            C_est, ...
            C_true, ...
            threshold);
        
        
        time_taken = ...
            toc(tStart);
        
        
        ALL_MCP(data_rep,:) = ...
            [ ...
            TPR, ...
            FPR, ...
            MCC, ...
            RMSE, ...
            MAD, ...
            time_taken ...
            ];
        
        
        Best_lambda_MCP(data_rep) = ...
            lambda_array(best_i);
        
        
        Best_gamma_MCP(data_rep) = ...
            gamma_array_MCP(best_j);
        
        
        fprintf('\n');
        
        fprintf( ...
            ['MCP rep %d selected: ' ...
            'lambda = %.6g, gamma = %.3f\n'], ...
            data_rep, ...
            Best_lambda_MCP(data_rep), ...
            Best_gamma_MCP(data_rep));
        
        fprintf( ...
            ['TPR = %.4f | FPR = %.4f | MCC = %.4f | ' ...
            'RMSE = %.4f | MAD = %.4f | Time = %.2f\n'], ...
            TPR,FPR,MCC,RMSE,MAD,time_taken);
        
    end
    
end


% ================================================================
% Save outputs
% ================================================================

metric_names = ...
    {'TPR','FPR','MCC','RMSE','MAD','Comp_time'};


output_folder = ...
    'Simulation output';


if ~exist(output_folder,'dir')
    mkdir(output_folder);
end


% ================================================================
% SCAD output
% ================================================================

if ismember(4,methods_to_perform)
    
    fprintf('\n');
    fprintf('====================================================\n');
    fprintf('LLA/MM SCAD RESULTS\n');
    fprintf('====================================================\n');
    
    
    fprintf('Best lambdas:\n');
    disp(Best_lambda_SCAD');
    
    
    fprintf('Best SCAD a values:\n');
    disp(Best_a_SCAD');
    
    
    fprintf( ...
        ['TPR, FPR, MCC, RMSE, MAD, Comp_time: ' ...
        '%.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n'], ...
        mean(ALL_SCAD,1,'omitnan'));
    
    
    T = ...
        array2table( ...
        ALL_SCAD, ...
        'VariableNames', ...
        metric_names);
    
    
    filename = ...
        sprintf( ...
        ['TPR_FPR_MCC_RMSE_MAD_time_' ...
        'p_%d_n_%d_C_%s_NumReps_%d_' ...
        'LLA_MM_SCAD.csv'], ...
        p,n,Ctype,num_dataset_rep);
    
    
    writetable( ...
        T, ...
        fullfile(output_folder,filename), ...
        'WriteRowNames',false);
    
    
    % ------------------------------------------------------------
    % Save selected tuning parameters as separate CSV
    % ------------------------------------------------------------
    
    T_tune = ...
        table( ...
        Best_lambda_SCAD, ...
        Best_a_SCAD, ...
        'VariableNames', ...
        {'Best_lambda','Best_a'});
    
    
    tune_filename = ...
        sprintf( ...
        ['BestTuning_' ...
        'p_%d_n_%d_C_%s_NumReps_%d_' ...
        'LLA_MM_SCAD.csv'], ...
        p,n,Ctype,num_dataset_rep);
    
    
    writetable( ...
        T_tune, ...
        fullfile(output_folder,tune_filename), ...
        'WriteRowNames',false);
    
end


% ================================================================
% MCP output
% ================================================================

if ismember(5,methods_to_perform)
    
    fprintf('\n');
    fprintf('====================================================\n');
    fprintf('LLA/MM MCP RESULTS\n');
    fprintf('====================================================\n');
    
    
    fprintf('Best lambdas:\n');
    disp(Best_lambda_MCP');
    
    
    fprintf('Best MCP gamma values:\n');
    disp(Best_gamma_MCP');
    
    
    fprintf( ...
        ['TPR, FPR, MCC, RMSE, MAD, Comp_time: ' ...
        '%.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n'], ...
        mean(ALL_MCP,1,'omitnan'));
    
    
    T = ...
        array2table( ...
        ALL_MCP, ...
        'VariableNames', ...
        metric_names);
    
    
    filename = ...
        sprintf( ...
        ['TPR_FPR_MCC_RMSE_MAD_time_' ...
        'p_%d_n_%d_C_%s_NumReps_%d_' ...
        'LLA_MM_MCP.csv'], ...
        p,n,Ctype,num_dataset_rep);
    
    
    writetable( ...
        T, ...
        fullfile(output_folder,filename), ...
        'WriteRowNames',false);
    
    
    % ------------------------------------------------------------
    % Save selected tuning parameters
    % ------------------------------------------------------------
    
    T_tune = ...
        table( ...
        Best_lambda_MCP, ...
        Best_gamma_MCP, ...
        'VariableNames', ...
        {'Best_lambda','Best_gamma'});
    
    
    tune_filename = ...
        sprintf( ...
        ['BestTuning_' ...
        'p_%d_n_%d_C_%s_NumReps_%d_' ...
        'LLA_MM_MCP.csv'], ...
        p,n,Ctype,num_dataset_rep);
    
    
    writetable( ...
        T_tune, ...
        fullfile(output_folder,tune_filename), ...
        'WriteRowNames',false);
    
end