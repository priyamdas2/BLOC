clear all;  clc;%close all;
clearvars;
addpath('./Simulation output/');


p = 50;
n = 100;
NumReps = 10;
method_name = 'BLOC_SCAD'; % 'BLOC_SCAD'/ 'BLOC_MCP'
Ctype = 'SparseUniform';   %'SparseUniform' / 'Blockdiag'


filename = sprintf('Simulation output/TPR_FPR_MCC_RMSE_MAD_time_p_%d_n_%d_C_%s_NumReps_%d_%s.csv', ...
            p, n, Ctype, NumReps, method_name);

summarize_reps_matrix(p, n, method_name, Ctype, NumReps, filename);