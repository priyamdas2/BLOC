clc; clear;
addpath('./Outputs/');

funcNames = {'ackley', 'griewank', 'rosenbrock', 'rastrigin'};
M_values = [5, 10, 20, 50];
num_reps = 10;
methods = {'BLOC', 'BLOCparallel', 'fmincon:active-set', 'fmincon:interior-point',...
    'fmincon:sqp', 'Manopt:barzilai-borwein', 'Manopt:conjugate-gradient',...
    'Manopt:steepest-descent', 'Manopt:trust-region'};
colNames = {'min value', 'value s.e.', 'mean time', 'time s.e.'};

% Generate the filenames
for i = 1:length(funcNames)
    for j = 1:length(M_values)
        fname = sprintf('All_funvals_%s_M_%d_reps_10.csv', funcNames{i}, M_values(j));
        funvals = readmatrix(fname);
        fname2 = sprintf('All_comp_times_%s_M_%d_reps_10.csv', funcNames{i}, M_values(j));
        comp_times = readmatrix(fname2);
        
        
        min_funvals = min(funvals);
        stderr_funvals = std(funvals) / sqrt(num_reps);
        mean_times = mean(comp_times);
        stderr_times = std(comp_times) / sqrt(num_reps);
        
        % Combine into a 9 x 4 matrix
        summary_matrix = [min_funvals' stderr_funvals' mean_times' stderr_times'];
        
        formatted_min_funvals = arrayfun(@(x) sprintf('%.2e', x), summary_matrix(:, 1), 'UniformOutput', false);
        formatted_stderr_funvals = arrayfun(@(x) sprintf('%.2e', x), summary_matrix(:, 2), 'UniformOutput', false);
        mean_times = summary_matrix(:, 3);
        truncated_mean_times = fix(mean_times * 100) / 100;
        stderr_times = summary_matrix(:, 4);
        truncated_stderr_times = fix(stderr_times * 1000) / 1000;
        
        formatted_time = arrayfun(@(m, s) sprintf('%.2f (%.3f)', m, s), mean_times, stderr_times, 'UniformOutput', false);

        % Create table with three columns
        summary_table = table(formatted_min_funvals, formatted_stderr_funvals, formatted_time, ...
            'VariableNames', {'min value', 'value s.e.', 'time [mean (se)]'}, ...
            'RowNames', methods);
        

        output_fname = sprintf('Outputs/Summary_%s_M_%d.csv', funcNames{i}, M_values(j));
        writetable(summary_table, output_fname, 'WriteRowNames', true);
        
    end
end


%%% High-dimensional


files = {
    'All_highdim_funvalNcomptimes_ackley_M_100_reps_10.csv', ...
    'All_highdim_funvalNcomptimes_griewank_M_100_reps_10.csv', ...
    'All_highdim_funvalNcomptimes_rosenbrock_M_100_reps_10.csv',...
    'All_highdim_funvalNcomptimes_rastrigin_M_100_reps_10.csv'
};

% Function names for row labels
funcs = {'ackley', 'griewank', 'rosenbrock', 'rastrigin'};

% Initialize output cell array
summary_cell = cell(4, 3);

for i = 1:4
    % Read data
    data = readmatrix(files{i});
    funvals = data(:,1);
    comptimes = data(:,2);

    % Calculate statistics
    min_val = min(funvals);
    se_val = std(funvals) / sqrt(length(funvals));
    mean_time = mean(comptimes);
    se_time = std(comptimes) / sqrt(length(comptimes));

    % Format values
    formatted_min_val = sprintf('%.2e', min_val);
    formatted_se_val = sprintf('%.2e', se_val);
    formatted_time = sprintf('%.2f (%.3f)', mean_time, se_time);

    % Store in cell array
    summary_cell{i,1} = formatted_min_val;
    summary_cell{i,2} = formatted_se_val;
    summary_cell{i,3} = formatted_time;
end

% Convert to table
summary_table = cell2table(summary_cell, ...
    'VariableNames', {'min_value', 'value_se', 'time_mean_se'}, ...
    'RowNames', funcs);

% Display result
disp(summary_table);

% Optional: save as CSV
writetable(summary_table, 'Outputs/Summary_high_dim.csv', 'WriteRowNames', true);
