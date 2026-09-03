clear; clc;

% ============================================================
% SUMMARIZE BLOC + COMPETITOR SIMULATION OUTPUTS
% ============================================================
%
% Each table entry is reported as:
%
%       mean (standard error)
%
%
% BLOC output files contain one 1 x 9 row per replication.
% COMPETITORS output files contain seven 1 x 9 rows per replication.
%
% ============================================================


% ============================================================
% USER SETTINGS
% ============================================================

input_dir  = 'outputs';
output_dir = 'outputs';

% Simulation setting whose files should be summarized
method_num  = 2;
p           = 20;
n           = 200;
contam_prop = 0.10;
K           = 500;       % same quantity called num_U in simulation code

% Replications to summarize
reps_to_read = 1:10;

% ------------------------------------------------------------
% Choose ANY or ALL of the following four cases
% true  = include
% false = exclude
% ------------------------------------------------------------

include_BLOC_SCAD = true;
include_BLOC_MCP  = true;

include_COMPETITORS_SCAD = true;
include_COMPETITORS_MCP  = true;


% ============================================================
% FIXED DEFINITIONS -- SHOULD MATCH THE SIMULATION FILES
% ============================================================

measure_names = { ...
    'Fro', ...
    'Spec', ...
    'Abs', ...
    'TPR', ...
    'FPR', ...
    'MCC', ...
    'PD', ...
    'Depth', ...
    'Time'};

% Number of decimal places used for each measure
% Order: Fro, Spec, Abs, TPR, FPR, MCC, PD, Depth, Time

mean_digits = [3,3,3,3,3,3,3,3,1];
se_digits   = [4,4,4,3,3,3,3,4,2];

% Must match the row order in the COMPETITORS simulation code
competitor_names = { ...
    'fmincon_active_set', ...
    'fmincon_interior_point', ...
    'fmincon_SQP', ...
    'Manopt_BB', ...
    'Manopt_CG', ...
    'Manopt_SD', ...
    'Manopt_TR'};

num_measures    = numel(measure_names);
num_competitors = numel(competitor_names);
num_reps        = numel(reps_to_read);

contam_pct = round(100*contam_prop);


% ============================================================
% CHECK THAT AT LEAST ONE CASE WAS REQUESTED
% ============================================================

if ~(include_BLOC_SCAD || include_BLOC_MCP || ...
     include_COMPETITORS_SCAD || include_COMPETITORS_MCP)

    error('At least one of the four include flags must be true.');

end


% ============================================================
% READ REQUESTED OUTPUT FILES
% ============================================================
%
% BLOC_SCAD_data / BLOC_MCP_data:
%       replication x measure
%
% COMP_SCAD_data / COMP_MCP_data:
%       replication x competitor x measure
% ============================================================

BLOC_SCAD_data = nan(num_reps,num_measures);
BLOC_MCP_data  = nan(num_reps,num_measures);

COMP_SCAD_data = nan(num_reps,num_competitors,num_measures);
COMP_MCP_data  = nan(num_reps,num_competitors,num_measures);


for rr = 1:num_reps

    rep = reps_to_read(rr);

    % --------------------------------------------------------
    % BLOC + SCAD
    % --------------------------------------------------------

    if include_BLOC_SCAD

        filename = fullfile( ...
            input_dir, ...
            sprintf( ...
            ['output_FroSpecAbsTprFprMccPDDepthTime_' ...
             'method_%d_p_%d_n_%d_contam_%02d_K_%d_' ...
             'BLOC_penalty_1_REP_%d.csv'], ...
             method_num,p,n,contam_pct,K,rep));

        BLOC_SCAD_data(rr,:) = ...
            read_bloc_file(filename,num_measures);

    end


    % --------------------------------------------------------
    % BLOC + MCP
    % --------------------------------------------------------

    if include_BLOC_MCP

        filename = fullfile( ...
            input_dir, ...
            sprintf( ...
            ['output_FroSpecAbsTprFprMccPDDepthTime_' ...
             'method_%d_p_%d_n_%d_contam_%02d_K_%d_' ...
             'BLOC_penalty_2_REP_%d.csv'], ...
             method_num,p,n,contam_pct,K,rep));

        BLOC_MCP_data(rr,:) = ...
            read_bloc_file(filename,num_measures);

    end


    % --------------------------------------------------------
    % COMPETITORS + SCAD
    % --------------------------------------------------------

    if include_COMPETITORS_SCAD

        filename = fullfile( ...
            input_dir, ...
            sprintf( ...
            ['output_FroSpecAbsTprFprMccPDDepthTime_' ...
             'method_%d_p_%d_n_%d_contam_%02d_K_%d_' ...
             'COMPETITORS_penalty_1_REP_%d.csv'], ...
             method_num,p,n,contam_pct,K,rep));

        COMP_SCAD_data(rr,:,:) = ...
            read_competitor_file( ...
            filename,num_competitors,num_measures);

    end


    % --------------------------------------------------------
    % COMPETITORS + MCP
    % --------------------------------------------------------

    if include_COMPETITORS_MCP

        filename = fullfile( ...
            input_dir, ...
            sprintf( ...
            ['output_FroSpecAbsTprFprMccPDDepthTime_' ...
             'method_%d_p_%d_n_%d_contam_%02d_K_%d_' ...
             'COMPETITORS_penalty_2_REP_%d.csv'], ...
             method_num,p,n,contam_pct,K,rep));

        COMP_MCP_data(rr,:,:) = ...
            read_competitor_file( ...
            filename,num_competitors,num_measures);

    end

end


% ============================================================
% BUILD SUMMARY TABLE
% ============================================================
%
% Requested row order:
%
%   BLOC rows first
%   then competitors in the same order as the simulation code
%
% If both penalties are included for a method, SCAD is placed first,
% followed immediately by MCP.
% ============================================================

row_names   = {};
summary_cell = cell(0,num_measures);


% ------------------------------------------------------------
% BLOC rows first
% ------------------------------------------------------------

if include_BLOC_SCAD

    row_names{end+1,1} = 'BLOC-SCAD';

    summary_cell(end+1,:) = ...
        summarize_matrix(BLOC_SCAD_data,mean_digits,se_digits);

end


if include_BLOC_MCP

    row_names{end+1,1} = 'BLOC-MCP';

    summary_cell(end+1,:) = ...
        summarize_matrix(BLOC_MCP_data,mean_digits,se_digits);

end


% ------------------------------------------------------------
% Competitor rows
% ------------------------------------------------------------

for m = 1:num_competitors

    if include_COMPETITORS_SCAD

        x = reshape( ...
            COMP_SCAD_data(:,m,:), ...
            num_reps,num_measures);

        row_names{end+1,1} = ...
            sprintf('%s-SCAD',competitor_names{m});

        summary_cell(end+1,:) = ...
            summarize_matrix(x,mean_digits,se_digits);

    end


    if include_COMPETITORS_MCP

        x = reshape( ...
            COMP_MCP_data(:,m,:), ...
            num_reps,num_measures);

        row_names{end+1,1} = ...
            sprintf('%s-MCP',competitor_names{m});

        summary_cell(end+1,:) = ...
            summarize_matrix(x,mean_digits,se_digits);

    end

end


% ============================================================
% CONVERT TO TABLE
% ============================================================
%
% p, n, method_num, contamination, K, etc. are NOT written
% as columns. They appear only in the output filename.
% ============================================================

Summary = cell2table( ...
    summary_cell, ...
    'VariableNames',measure_names);

Summary.Properties.RowNames = row_names;


% ============================================================
% CREATE OUTPUT FILENAME
% ============================================================
%
% Examples:
%
% Only BLOC-SCAD:
% ALL_summary_method_3_p_20_n_200_contam_10_K_500_BLOC_penalty_1.csv
%
% If multiple cases are selected, every selected case is appended
% to the filename.
% ============================================================

selected_tags = {};

if include_BLOC_SCAD
    selected_tags{end+1} = 'BLOC_penalty_1';
end

if include_BLOC_MCP
    selected_tags{end+1} = 'BLOC_penalty_2';
end

if include_COMPETITORS_SCAD
    selected_tags{end+1} = 'COMPETITORS_penalty_1';
end

if include_COMPETITORS_MCP
    selected_tags{end+1} = 'COMPETITORS_penalty_2';
end

case_tag = strjoin(selected_tags,'_');

summary_filename = fullfile( ...
    output_dir, ...
    sprintf( ...
    ['ALL_summary_method_%d_p_%d_n_%d_contam_%02d_K_%d_%s.csv'], ...
    method_num,p,n,contam_pct,K,case_tag));


% ============================================================
% WRITE SUMMARY
% ============================================================

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

writetable( ...
    Summary, ...
    summary_filename, ...
    'WriteRowNames',true);


fprintf('\n====================================================\n');
fprintf('SUMMARY WRITTEN TO:\n%s\n',summary_filename);
fprintf('====================================================\n\n');

disp(Summary);


% ============================================================
% LOCAL FUNCTIONS
% ============================================================

function out = read_bloc_file(filename,num_measures)

    out = nan(1,num_measures);

    if ~isfile(filename)

        warning('Missing BLOC file:\n%s',filename);
        return;

    end

    A = readmatrix(filename);

    if isempty(A)

        warning('Empty BLOC file:\n%s',filename);
        return;

    end

    if size(A,2) < num_measures

        warning( ...
            ['BLOC file has only %d columns; ' ...
             '%d were expected:\n%s'], ...
             size(A,2),num_measures,filename);
        return;

    end

    % BLOC should contain one row.
    % If more than one row is present, use the first row and warn.
    if size(A,1) > 1

        warning( ...
            ['BLOC file has %d rows; expected 1. ' ...
             'Using the first row:\n%s'], ...
             size(A,1),filename);

    end

    out = A(1,1:num_measures);

end


function out = read_competitor_file( ...
    filename,num_competitors,num_measures)

    out = nan(num_competitors,num_measures);

    if ~isfile(filename)

        warning('Missing COMPETITORS file:\n%s',filename);
        return;

    end

    A = readmatrix(filename);

    if isempty(A)

        warning('Empty COMPETITORS file:\n%s',filename);
        return;

    end

    if size(A,2) < num_measures

        warning( ...
            ['COMPETITORS file has only %d columns; ' ...
             '%d were expected:\n%s'], ...
             size(A,2),num_measures,filename);
        return;

    end

    rows_available = min(size(A,1),num_competitors);

    out(1:rows_available,:) = ...
        A(1:rows_available,1:num_measures);

    if size(A,1) < num_competitors

        warning( ...
            ['COMPETITORS file has only %d rows; ' ...
             '%d were expected. Missing rows remain NaN:\n%s'], ...
             size(A,1),num_competitors,filename);

    elseif size(A,1) > num_competitors

        warning( ...
            ['COMPETITORS file has %d rows; ' ...
             'only the first %d will be used:\n%s'], ...
             size(A,1),num_competitors,filename);

    end

end


function summary_row = summarize_matrix(X,mean_digits,se_digits)

    num_measures = size(X,2);

    summary_row = cell(1,num_measures);

    for j = 1:num_measures

        x = X(:,j);

        % If ANY replication is NaN or Inf,
        % keep the summary for this measure as NaN.
        if any(~isfinite(x))

            summary_row{j} = 'NaN (NaN)';
            continue;

        end

        mu = mean(x);

        if numel(x) >= 2
            se = std(x,0) / sqrt(numel(x));
        else
            se = NaN;
        end

        % Metric-specific formatting
        format_string = sprintf( ...
            '%%.%df (%%.%df)', ...
            mean_digits(j), ...
            se_digits(j));

        summary_row{j} = ...
            sprintf(format_string,mu,se);

    end

end