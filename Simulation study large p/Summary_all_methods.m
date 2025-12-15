clear all;  clc;%close all;
clearvars;
addpath('./pdlc/');
addpath('./BLOC/');
addpath('./helper funs/');
addpath('./outputs/');

method_num = 3;
p = 100;
n = 50;
Num_exps = 10;


filename = sprintf(['output_FroSpecAbsTprFprMcc_method_%d_p_%d_n_%d_' ...
                    'NumExp_%d_BLOC_penaly_1.csv'], ...
                    method_num, p, n, Num_exps);
output_BLOC_SCAD = csvread(filename);
[mean_vals_01, se_vals_01] = summary_stats(output_BLOC_SCAD);

filename = sprintf(['output_FroSpecAbsTprFprMcc_method_%d_p_%d_n_%d_' ...
                    'NumExp_%d_BLOC_penaly_2.csv'], ...
                    method_num, p, n, Num_exps);
output_BLOC_MCP = csvread(filename);
[mean_vals_02, se_vals_02] = summary_stats(output_BLOC_MCP);

filename = sprintf(['output_FroSpecAbsTprFprMcc_method_%d_p_%d_n_%d_' ...
                    'NumExp_%d_1L1_ADMM.csv'], ...
                    method_num, p, n, Num_exps);
output_L1_ADMM = csvread(filename);
[mean_vals_1, se_vals_1] = summary_stats(output_L1_ADMM);

filename = sprintf(['output_FroSpecAbsTprFprMcc_method_%d_p_%d_n_%d_' ...
                    'NumExp_%d_2Lq_BCD_half.csv'], ...
                    method_num, p, n, Num_exps);
output_Lq_BCD_half = csvread(filename);
[mean_vals_2, se_vals_2] = summary_stats(output_Lq_BCD_half);

filename = sprintf(['output_FroSpecAbsTprFprMcc_method_%d_p_%d_n_%d_' ...
                    'NumExp_%d_3Lq_BCD_zero.csv'], ...
                    method_num, p, n, Num_exps);
output_Lq_BCD_zero = csvread(filename);
[mean_vals_3, se_vals_3] = summary_stats(output_Lq_BCD_zero);

filename = sprintf(['output_FroSpecAbsTprFprMcc_method_%d_p_%d_n_%d_' ...
                    'NumExp_%d_4SCAD_BCD.csv'], ...
                    method_num, p, n, Num_exps);
output_SCAD_BCD = csvread(filename);
[mean_vals_4, se_vals_4] = summary_stats(output_SCAD_BCD);

filename = sprintf(['output_FroSpecAbsTprFprMcc_method_%d_p_%d_n_%d_' ...
                    'NumExp_%d_5IRW_ADMM.csv'], ...
                    method_num, p, n, Num_exps);
output_IRW_ADMM = csvread(filename);
[mean_vals_5, se_vals_5] = summary_stats(output_IRW_ADMM);

filename = sprintf(['output_FroSpecAbsTprFprMcc_method_%d_p_%d_n_%d_' ...
                    'NumExp_%d_6Trad_IRW.csv'], ...
                    method_num, p, n, Num_exps);
output_Trad_IRW = csvread(filename);
[mean_vals_6, se_vals_6] = summary_stats(output_Trad_IRW);

filename = sprintf(['output_FroSpecAbsTprFprMcc_method_%d_p_%d_n_%d_' ...
                    'NumExp_%d_7soft_thres.csv'], ...
                    method_num, p, n, Num_exps);
output_soft_thres = csvread(filename);
[mean_vals_7, se_vals_7] = summary_stats(output_soft_thres);

filename = sprintf(['output_FroSpecAbsTprFprMcc_method_%d_p_%d_n_%d_' ...
                    'NumExp_%d_8hard_thres.csv'], ...
                    method_num, p, n, Num_exps);
output_hard_thres = csvread(filename);
[mean_vals_8, se_vals_8] = summary_stats(output_hard_thres);

filename = sprintf(['output_FroSpecAbsTprFprMcc_method_%d_p_%d_n_%d_' ...
                    'NumExp_%d_9SCAD_thres.csv'], ...
                    method_num, p, n, Num_exps);
output_SCAD_thres = csvread(filename);
[mean_vals_9, se_vals_9] = summary_stats(output_SCAD_thres);


summary_mean = [mean_vals_01; mean_vals_02; mean_vals_1; mean_vals_2; mean_vals_3; mean_vals_4; mean_vals_5;
    mean_vals_6; mean_vals_7; mean_vals_8; mean_vals_9];

summary_stderr = [se_vals_01; se_vals_02; se_vals_1; se_vals_2; se_vals_3; se_vals_4; se_vals_5;
    se_vals_6; se_vals_7; se_vals_8; se_vals_9];

filename = sprintf(['outputs/Summary_ALL_methods_MEAN_method_%d_p_%d_n_%d.csv'], ...
                    method_num, p, n);
dlmwrite(filename, round(summary_mean, 3), 'precision', '%.3f');


filename = sprintf(['outputs/Summary_ALL_methods_SE_method_%d_p_%d_n_%d.csv'], ...
                    method_num, p, n);
dlmwrite(filename, round(summary_stderr, 3), 'precision', '%.3f');

out = format_mean_se_csv(method_num, p, n);