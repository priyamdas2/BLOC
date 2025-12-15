rng(1)

% -------------------- USER SETTINGS --------------------
p = 100;                      % number of variables
n = 100;                      % number of samples
num_data_rep = 10;
structure_options = [1,2];  % choose from: 1 = Blockdiag, 2= SparseUniform
threshold = 1e-2; % don't change

output_folder = 'Simulation data';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end


for i = 1:length(structure_options)
    opt = structure_options(i);
    rand_seed = 1;
    switch opt
        case 1
            C_true = generate_corr_Blockdiag(p, rand_seed);  % blocks of 5, p must be divisible by 5
            name = 'Blockdiag';
            P = covering_Blockdiag(C_true, threshold);
        case 2
            sparsity = 0.95; min_val = 0.3; max_val = 0.6; num_tries = 100;
            if(p > 20)
                sparsity = 0.98;
            end
            if(p > 50)
                sparsity = 0.99;
            end
            C_true = generate_corr_SparseUniform(p, sparsity, min_val, max_val, num_tries);
            name = 'SparseUniform';
            rand_seed_here = 123;
            P = covering_SparseUniform(C_true, threshold, rand_seed_here);
        otherwise
            error('Unknown structure option: %d', opt);
    end
    filename_C = fullfile(output_folder, sprintf('C_p_%d_n_%d_C_%s.csv', p, n, name));
    writematrix(C_true, filename_C);
    filename_P = fullfile(output_folder, sprintf('P_cover_p_%d_n_%d_%s.csv', p, n, name));
    writematrix(P, filename_P);
    
    for data_rep = 1:num_data_rep
        rng(data_rep)
        
        % Simulate data
        mu_true = zeros(p,1);
        X_raw = mvnrnd(mu_true, C_true, n);

        filename_X = fullfile(output_folder, sprintf('X_p_%d_n_%d_C_%s_DataRep_%d.csv', p, n, name, data_rep));
        writematrix(X_raw, filename_X);
    end
end



