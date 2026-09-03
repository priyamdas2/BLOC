function outPath = format_mean_se_csv(method, p, n, baseDir)
% FORMAT_MEAN_SE_CSV
% Combine MEAN and SE CSVs into a formatted CSV.
%
% Reads:
%   Summary_ALL_methods_MEAN_method_<method>_p_<p>_n_<n>.csv
%   Summary_ALL_methods_SE_method_<method>_p_<p>_n_<n>.csv
%
% Writes:
%   Summary_ALL_methods_FORMATTED_method_<method>_p_<p>_n_<n>.csv
%
% -------------------------------------------------------------------------
% ROW ORDER MUST MATCH summary_mean / summary_stderr construction:
%
%   1.  BLOC_MCP
%   2.  BLOC_SCAD
%   3.  WeiZhaoMCP
%   4.  WeiZhaoMCP_Cov
%   5.  L1_ADMM
%   6.  Lq_BCD_half
%   7.  Lq_BCD_zero
%   8.  SCAD_BCD
%   9.  IRW_ADMM
%   10. Trad_IRW
%   11. soft_thres
%   12. hard_thres
%   13. SCAD_thres
%
% -------------------------------------------------------------------------
% COLUMN ORDER in input:
%
%   1. Frobenius norm
%   2. Spectral norm
%   3. Absolute error
%   4. TPR
%   5. FPR
%   6. MCC
%   7. PD
%   8. Computation time
%
% Formatting currently retained:
%   - Means: 2 decimals, except column 3 -> 3 decimals
%   - SEs:   3 decimals
% -------------------------------------------------------------------------

    if nargin < 4 || isempty(baseDir)
        baseDir = pwd;
    end


    % =====================================================================
    % Filenames
    % =====================================================================

    meanFile = fullfile(baseDir, sprintf( ...
        'Summary_ALL_methods_MEAN_method_%d_p_%d_n_%d.csv', ...
        method, p, n));

    seFile = fullfile(baseDir, sprintf( ...
        'Summary_ALL_methods_SE_method_%d_p_%d_n_%d.csv', ...
        method, p, n));

    outPath = fullfile(baseDir, sprintf( ...
        'Summary_ALL_methods_FORMATTED_method_%d_p_%d_n_%d.csv', ...
        method, p, n));


    if ~isfile(meanFile)
        error('Mean file not found: %s', meanFile);
    end

    if ~isfile(seFile)
        error('SE file not found: %s', seFile);
    end


    % =====================================================================
    % Read input matrices
    % =====================================================================

    M = readmatrix(meanFile);
    S = readmatrix(seFile);


    if isempty(M) || isempty(S)
        error('One of the input files is empty or could not be parsed.');
    end


    % =====================================================================
    % Align dimensions by position
    % =====================================================================

    nRows = min(size(M,1), size(S,1));
    nCols = min(size(M,2), size(S,2));

    M = M(1:nRows, 1:nCols);
    S = S(1:nRows, 1:nCols);


    % =====================================================================
    % EXACT ROW ORDER FROM summary_mean / summary_stderr
    % =====================================================================

    rowNames = {
        'BLOC_MCP'
        'BLOC_SCAD'
        'WeiZhaoMCP'
        'WeiZhaoMCP_Cov'
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


    % Sanity check
    if nRows ~= numel(rowNames)

        warning( ...
            ['Number of rows in summary files (%d) does not match ' ...
             'expected number of methods (%d).'], ...
             nRows, numel(rowNames));

    end


    % Trim or pad only if necessary
    if numel(rowNames) < nRows

        for k = numel(rowNames)+1:nRows
            rowNames{k,1} = sprintf('Method_%d', k);
        end

    elseif numel(rowNames) > nRows

        rowNames = rowNames(1:nRows);

    end


    % =====================================================================
    % Build formatted output
    %
    % Column 1 = method
    % Remaining columns = "mean (SE)"
    % =====================================================================

    outCell = cell(nRows, nCols + 1);

    outCell(:,1) = rowNames;


    for j = 1:nCols

        for i = 1:nRows

            m = M(i,j);
            s = S(i,j);


            % -------------------------------------------------------------
            % Mean formatting
            %
            % Column order:
            %
            % 1 Fro
            % 2 Spec
            % 3 Abs
            % 4 TPR
            % 5 FPR
            % 6 MCC
            % 7 PD
            % 8 Time
            % -------------------------------------------------------------

            if j == 3
                mStr = sprintf('%.3f', m);
            else
                mStr = sprintf('%.2f', m);
            end


            % SE always 3 decimals
            sStr = sprintf('%.3f', s);


            % -------------------------------------------------------------
            % NaN handling
            % -------------------------------------------------------------

            if ~isfinite(m) && ~isfinite(s)

                outCell{i,j+1} = '';

            elseif isfinite(m) && ~isfinite(s)

                outCell{i,j+1} = mStr;

            elseif ~isfinite(m) && isfinite(s)

                outCell{i,j+1} = ['(' sStr ')'];

            else

                outCell{i,j+1} = ...
                    [mStr ' (' sStr ')'];

            end

        end

    end


    % =====================================================================
    % Write CSV
    %
    % No column headers, preserving your current output convention.
    % =====================================================================

    writecell(outCell, outPath);

    fprintf('Wrote formatted CSV: %s\n', outPath);

    disp(outCell);

end