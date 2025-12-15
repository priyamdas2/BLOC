clearvars;
warning('off', 'all');
addpath('./BLOC/');
addpath('./Supp functions/');

case_id = 5;   % 1=BRCA, 2=CESC, 3=OV, 4=UCEC, 5=UCS
case_names = {'BRCA','CESC','OV','UCEC','UCS'};
filename = ['proteomics_' case_names{case_id} '.csv'];
X_raw = readmatrix(filename);

threshold = 1e-2;
n = size(X_raw,1);
p = size(X_raw,2);

T = readtable('protein_pathways_selected.csv');

% Take the 2nd column (pathway names) as strings
pathways = string(T{:,2});
nn = numel(pathways);

% Initialize block-diagonal 0/1 matrix
P = ones(nn,nn);

% Find start indices of each consecutive run (block)
run_starts = [1; find(pathways(1:end-1) ~= pathways(2:end)) + 1];
run_ends   = [run_starts(2:end)-1; nn];

% Fill each block with ones
for k = 1:numel(run_starts)
    ss = run_starts(k);
    ee = run_ends(k);
    P(ss:ee, ss:ee) = 0;
end

% ------------------- BLOC Inputs -------------------
doparallel = 0; % set 0 if doing parfor over data reps
%maxNumCompThreads(1);
MaxTime = 3600;
MaxRuns = 10;
MaxIter = 5000;
TolFun1 = 1e-4;
TolFun2 = 1e-2;
phi = 1e-6;
DisplayUpdate = 1;
DisplayEvery = 2;

% ------------------- Lambda grid ---------------------
log10lambdaLB = -2;
log10lambdaUB = 1;
NumLambdas = 10;
lambda_array = 10.^linspace(log10lambdaLB, log10lambdaUB, NumLambdas);
a_array_SCAD = 3.7;


% ------------------- SCAD with BLOC ---------------------
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


% ------------------- saving output ---------------------
filename = ['corr_' case_names{case_id} '.csv'];
% Save C_est as CSV (matrix only, no row/col names)
writematrix(C_est, filename);
