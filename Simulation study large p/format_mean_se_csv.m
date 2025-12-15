function outPath = format_mean_se_csv(method, p, n, baseDir)
% FORMAT_MEAN_SE_CSV  Combine MEAN and SE CSVs into a formatted CSV.
%
%   outPath = format_mean_se_csv(method, p, n)
%   outPath = format_mean_se_csv(method, p, n, baseDir)
%
% Reads:
%   Summary_ALL_methods_MEAN_method_<method>_p_<p}_n_<n>.csv
%   Summary_ALL_methods_SE_method_<method>_p_<p}_n_<n>.csv
%
% Writes:
%   Summary_ALL_methods_FORMATTED_method_<method>_p_<p}_n_<n>.csv
%
% Rules:
%   - Drop the LAST column from both inputs
%   - Means: round to 2 decimals for all columns EXCEPT the 3rd column -> 3 decimals
%   - SEs: round to 3 decimals, shown in parentheses
%   - Inputs have NO column names
%   - Output includes row names:
%       BLOC_SCAD, BLOC_MCP, L1_ADMM, Lq_BCD_half, Lq_BCD_zero,
%       SCAD_BCD, IRW_ADMM, Trad_IRW, soft_thres, hard_thres, SCAD_thres

    if nargin < 4 || isempty(baseDir)
        baseDir = pwd;
    end

    % ----- Filenames
    meanFile = fullfile(baseDir, sprintf('outputs/Summary_ALL_methods_MEAN_method_%d_p_%d_n_%d.csv', method, p, n));
    seFile   = fullfile(baseDir, sprintf('outputs/Summary_ALL_methods_SE_method_%d_p_%d_n_%d.csv',   method, p, n));
    outPath  = fullfile(baseDir, sprintf('outputs/Summary_ALL_methods_FORMATTED_method_%d_p_%d_n_%d.csv', method, p, n));

    if ~isfile(meanFile)
        error('Mean file not found: %s', meanFile);
    end
    if ~isfile(seFile)
        error('SE file not found: %s', seFile);
    end

    % ----- Read (no headers)
    M = readmatrix(meanFile);
    S = readmatrix(seFile);

    if isempty(M) || isempty(S)
        error('One of the input files is empty or could not be parsed.');
    end

    % ----- Drop last column from both
    if size(M,2) >= 2, M = M(:, 1:end-1); end
    if size(S,2) >= 2, S = S(:, 1:end-1); end

    % ----- Align by position (min common width & height)
    nRows = min(size(M,1), size(S,1));
    nCols = min(size(M,2), size(S,2));
    M = M(1:nRows, 1:nCols);
    S = S(1:nRows, 1:nCols);

    % ----- Row names (as requested)
    rowNames = {
        'BLOC_SCAD'
        'BLOC_MCP'
        'L1_ADMM'
        'Lq_BCD_half'
        'Lq_BCD_zero'
        'SCAD_BCD'
        'IRW_ADMM'
        'Trad_IRW'
        'soft_thres'
        'hard_thres'
        'SCAD_thres'
    };

    % Trim or pad row names to match nRows
    if numel(rowNames) < nRows
        % pad with generic names if there are more rows than expected
        rowNames(end+1:nRows) = arrayfun(@(k) sprintf('Method_%d', k), numel(rowNames)+1:nRows, 'UniformOutput', false);
    elseif numel(rowNames) > nRows
        rowNames = rowNames(1:nRows);
    end

    % ----- Build formatted cell array: first col = row names, others = "mean (SE)"
    outCell = cell(nRows, nCols + 1);
    outCell(:,1) = rowNames;

    for j = 1:nCols
        for i = 1:nRows
            m = M(i,j);
            s = S(i,j);

            % Rounding/format rules
            if j == 3
                mStr = sprintf('%.3f', m);
            else
                mStr = sprintf('%.2f', m);
            end
            sStr = sprintf('%.3f', s);

            % Handle NaNs gracefully
            if ~isfinite(m) && ~isfinite(s)
                outCell{i, j+1} = '';
            elseif isfinite(m) && ~isfinite(s)
                outCell{i, j+1} = mStr;
            elseif ~isfinite(m) && isfinite(s)
                outCell{i, j+1} = ['(' sStr ')'];
            else
                outCell{i, j+1} = [mStr ' (' sStr ')'];
            end
        end
    end

    % ----- Write CSV (no column headers)
    writecell(outCell, outPath);
    fprintf('Wrote formatted CSV: %s\n', outPath);
end
