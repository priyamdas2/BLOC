clear all;  clc;%close all;
clearvars;
addpath('./pdlc/');
addpath('./BLOC/');
addpath('./helper funs/');
%addpath('./Simulation data/');

p = 50; % covariance dimension 20, 50, 100
method_num = 3;
SIGMA = cov_model(p, method_num);
penalty = 1; % 1 = SCAD, 2 = MCP
n = 50;


Num_exps = 1;
tol = 1e-3; % sparsity detection


n2 = round(n/log(n)); % Test
n1 = n - n2;          % Train
num_CV = 5; % 5
num_lambda = 15; % 15
lamdas = logspace(-3, 0, num_lambda);
C0 = eye(p); % starting point


%%% BLOC parameters (for CV)

params.MaxTime = 3;
params.MaxRuns = 1;
params.MaxIter = 1000;
params.TolFun1 = 1e-4;
params.TolFun2 = 1e-2;
params.phi = 1e-6;
params.DisplayUpdate = 1;
params.DisplayEvery = 10;

%%% BLOC parameters (for Final)

paramsFinal.MaxTime = 100;
paramsFinal.MaxRuns = 20;
paramsFinal.MaxIter = 1000;
paramsFinal.TolFun1 = 1e-4;
paramsFinal.TolFun2 = 1e-2;
paramsFinal.phi = 1e-6;
paramsFinal.DisplayUpdate = 1;
paramsFinal.DisplayEvery = 10;

%%%%%%%%%%%%%%%%%%%%

output_BLOC = nan(Num_exps, 7);

for k=1:Num_exps
    rng(k)
    fprintf('=> Performing exp number : %d. \n', k);
    r = mvnrnd(zeros(p,1), SIGMA, n);
    S = cov(r);%
    STD_S = diag(S).^(0.5); %std
    S_corr = diag(1./STD_S) * S * diag(1./STD_S);
    

    Theta = eye(p);
    FroErr = nan(num_CV,num_lambda);
    for kk=1:num_CV
        J  = randperm(n);     % m randomly chosen indices
        S1 = corr(r(J(1:n1),:));
        S2 = corr(r(J(n1+1:end),:));
        parfor l=1:num_lambda
            fprintf('=> Performing exp number : %d, CV number: %d, lambda grid number: %d. \n', k, kk, l);
            lambda = lamdas(l);
            if (penalty == 1)
                a = 3.7;
                objFun = @(Theta) obj_scad_corr(Theta, S1, lambda, a);
            else
                gamma = 3;
                objFun = @(Theta) obj_mcp_corr(Theta, S1, lambda, gamma);
            end

            [C_opt, fval, comp_time] = BLOC(objFun, C0, params);
            FroErr(kk,l) = norm(C_opt-S2,'fro');
        end
    end
    [~, mi] = min(sum(FroErr));
    lamda_opt = lamdas(mi);

    if (penalty == 1)
        a = 3.7;
        objFun = @(Theta) obj_scad_corr(Theta, S1, lamda_opt, a);
    else
        gamma = 3;
        objFun = @(Theta) obj_mcp_corr(Theta, S1, lamda_opt, gamma);
    end
    
    [C_final, fval, comp_time_final] = BLOCparallel(objFun, C0, paramsFinal);

    ErrFro = norm(C_final-SIGMA,'fro');
    ErrSpe = norm(C_final-SIGMA);
    err_abs = mean(abs(C_final(:) - SIGMA(:)));

    [TPR, FPR, MCC] = TPR_FPR_MCC(C_final, SIGMA, tol);
    output = [ErrFro, ErrSpe, err_abs, TPR, FPR, MCC, comp_time_final];
    output_BLOC(k,:) = output;
 
end


filename = sprintf(['outputs/output_FroSpecAbsTprFprMcc_method_%d_p_%d_n_%d_' ...
                    'NumExp_%d_BLOC_penalty_%d.csv'], method_num, p, n, Num_exps, penalty);
                
                
csvwrite(filename, output_BLOC);
%     filename = sprintf(['outputs/output_FroSpecAbsTprFprMcc_method_%d_p_%d_n_%d_' ...
%                     'NumExp_%d_2Lq_BCD_half.csv'], ...
%                     method_num, p, n, Num_exps);
%     csvwrite(filename, output_Lq_BCD_half);


