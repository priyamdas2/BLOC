function [col_mean, col_stderr] = summary_stats(mat)
% summary_stats - compute column-wise mean and standard error
%
% Inputs:
%   mat         : matrix (rows = replicates, cols = metrics)
%
% Outputs:
%   col_mean    : 1 x numCols vector of column means
%   col_stderr  : 1 x numCols vector of column standard errors

    % Number of replicates
    n = size(mat, 1);

    % Column means
    col_mean = mean(mat, 1);

    % Standard error = std / sqrt(n)
    col_stderr = std(mat, 0, 1) / sqrt(n);
end
