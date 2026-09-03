clear all;  clc;%close all;
clearvars;
addpath('./pdlc/');
addpath('./BLOC/');
addpath('./helper funs/');

method_num = 3;
p = 200;
n = 100;
Num_exps = 10;


filename = sprintf(['outputs/output_FroSpecAbsTprFprMccPDTime_method_%d_p_%d_n_%d_' ...
                    'NumExp_%d_WeiZhaoMCP.csv'], ...
                    method_num, p, n, Num_exps);

output_WeiZhaoMCP = csvread(filename);

[mean_vals_WeiZhaoMCP, se_vals_WeiZhaoMCP] = ...
    summary_stats(output_WeiZhaoMCP);


summary_mean = mean_vals_WeiZhaoMCP;
summary_stderr = se_vals_WeiZhaoMCP;


filename = sprintf(['outputs/Summary_WeiZhaoMCP_MEAN_method_%d_p_%d_n_%d.csv'], ...
                    method_num, p, n);

dlmwrite(filename, round(summary_mean, 3), 'precision', '%.3f');


filename = sprintf(['outputs/Summary_WeiZhaoMCP_SE_method_%d_p_%d_n_%d.csv'], ...
                    method_num, p, n);

dlmwrite(filename, round(summary_stderr, 3), 'precision', '%.3f');