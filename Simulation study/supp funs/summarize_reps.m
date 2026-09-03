function summarize_reps(p, n, method, cond, numReps, outputFile)
    % p, n       : dimensions
    % method     : e.g. 'BLOC_SCAD', 'BLOC_MCP'
    % cond       : condition string (e.g. 'SparseUniform', 'Blockdiag')
    % numReps    : number of replicates (e.g. 10)
    % inputDir   : folder where raw CSV files are located
    % outputFile : name of summary CSV to save

    all_results = [];

    for r = 1:numReps
        % Construct file name (adjust pattern if needed)
        fname = sprintf('Simulation output/TPR_FPR_MCC_RMSE_MAD_time_BestParams_%s_p_%d_n_%d_C_%s_RepNum_%d.csv', ...
                         method, p, n, cond, r);

        % Read the file (ignore last 2 columns)
        data = readmatrix(fname);
        data = data(:,1:6);

        all_results = [all_results; data];
    end

    % Compute mean across replicates
    summary_mean = mean(all_results,1);

    % Optional: also compute standard error
    summary_se = std(all_results,0,1) / sqrt(numReps);

    % Combine into one table
    summary_table = array2table([summary_mean; summary_se], ...
        'VariableNames', {'TPR','FPR','MCC','RMSE','MAD','Comp_time'}, ...
        'RowNames', {'Mean','SE'});

    % Save
    writetable(summary_table, outputFile, 'WriteRowNames', true);
end
