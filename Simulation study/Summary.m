% -------- CONFIG --------
baseDir   = 'Simulation output';  % folder containing the CSVs
p_list    = [20, 50, 100];
n_list    = [50, 100, 150];
conds     = {'Blockdiag','SparseUniform'};
methods   = {'BLOC_MCP','BLOC_SCAD','SpCov'};
headers   = {'TPR','FPR','MCC','RMSE','MAD'};

% -------- SCRIPT --------
for idx = 1:numel(p_list)
    p = p_list(idx);
    n = n_list(idx);

    for c = 1:numel(conds)
        cond = conds{c};

        % Container for the 3x5 summary (strings "mean (se)")
        out = strings(numel(methods), numel(headers));

        for m = 1:numel(methods)
            method = methods{m};
            fname = fullfile(baseDir, sprintf( ...
                'TPR_FPR_MCC_RMSE_MAD_time_p_%d_n_%d_C_%s_NumReps_10_%s.csv', ...
                p, n, cond, method));

            if ~isfile(fname)
                % If missing file, fill the entire row with NaN
                out(m, :) = "NaN";
                continue;
            end

            data = readmatrix(fname);
            if isempty(data) || size(data,2) < 5
                % Bad or empty data -> fill with NaN
                out(m, :) = "NaN";
                continue;
            end

            % Keep only the first 5 columns: TPR, FPR, MCC, RMSE, MAD
            X = data(:, 1:5);

            % Compute mean and standard error over reps
            nreps = size(X,1);
            mu = mean(X, 1, 'omitnan');
            se = std(X, 0, 1, 'omitnan') ./ sqrt(nreps);

            % Format: TPR/FPR/MCC -> 2 decimals; RMSE/MAD -> 3 decimals; SE -> 3 decimals
            row_str = strings(1,5);
            row_str(1) = sprintf('%.2f (%.3f)', mu(1), se(1)); % TPR
            row_str(2) = sprintf('%.2f (%.3f)', mu(2), se(2)); % FPR
            row_str(3) = sprintf('%.2f (%.3f)', mu(3), se(3)); % MCC
            row_str(4) = sprintf('%.3f (%.3f)', mu(4), se(4)); % RMSE
            row_str(5) = sprintf('%.3f (%.3f)', mu(5), se(5)); % MAD

            out(m, :) = row_str;
        end

        % Convert to table for clean headers & row names
        T = array2table(out, 'VariableNames', headers, 'RowNames', methods);

        % Show in console
        fprintf('\nSummary for p=%d, n=%d, C=%s\n', p, n, cond);
        disp(T)

        % Write CSV
        outName = sprintf('Summary_ALL_p_%d_n_%d_C_%s.csv', p, n, cond);
        writetable(T, outName, 'WriteRowNames', true);
    end
end
