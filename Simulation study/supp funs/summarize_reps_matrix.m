function summarize_reps_matrix(p, n, method, cond, numReps, outputFile)
% Summarize replicate-wise raw CSVs into a 10x6 matrix with headers:
% Columns: TPR, FPR, MCC, RMSE, MAD, Comp_time
%
% Parameters
%   p, n        : dimensions
%   method      : e.g., 'BLOC_SCAD', 'BLOC_MCP'
%   cond        : condition string, e.g., 'SparseUniform', 'Blockdiag'
%   numReps     : number of replicates (e.g., 10)
%   outputFile  : path for the summary CSV to write
%
% Expected raw filename pattern (adjust if needed):
%   TPR_FPR_MCC_RMSE_MAD_time_BestParams_<method>_p_<p>_n_<n>_C_<cond>_RepNum_<r>.csv
%
% Notes:
%   - Ignores the last two columns in each raw CSV (uses only first 6).
%   - Writes a CSV with headers and exactly numReps rows (one per replicate).


    all_results = nan(numReps, 6);

    for r = 1:numReps
        fname = sprintf('Simulation output/TPR_FPR_MCC_RMSE_MAD_time_BestParams_%s_p_%d_n_%d_C_%s_RepNum_%d.csv', ...
            method, p, n, cond, r);

        if ~isfile(fname)
            error('File not found: %s', fname);
        end

        % Read as numeric; tolerate headers/no-headers
        data = readmatrix(fname);

        if isempty(data)
            error('Empty or unreadable data in file: %s', fname);
        end

        % Use only the first 6 columns (ignore last two)
        if size(data, 2) < 6
            error('File %s has fewer than 6 columns.', fname);
        end

        row = data(1, 1:6);  % assume 1 row per file; take first row if multiple
        all_results(r, :) = row;
    end

    % Turn into a table to preserve headers in the CSV
    summary_table = array2table(all_results, ...
        'VariableNames', {'TPR','FPR','MCC','RMSE','MAD','Comp_time'});

    % Write CSV with column headers, no row names
    writetable(summary_table, outputFile);
end
