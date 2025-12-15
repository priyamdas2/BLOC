function Corr = generate_corr_Blockdiag(p, rand_seed)
%GENERATE_CORR_BLOCKDIAG Generates a block-diagonal correlation matrix
% Each block is a 5x5 random correlation matrix.
%
% Inputs:
%   p         : total dimension (must be divisible by 5)
%   rand_seed: base random seed for reproducibility (optional)
%
% Output:
%   Corr      : p x p block-diagonal correlation matrix

    if nargin < 2
        rand_seed = 1;
    end

    block_size = 5;

    % Check divisibility
    if mod(p, block_size) ~= 0
        error('Total dimension p must be divisible by block size (%d).', block_size);
    end

    num_blocks = p / block_size;
    Corr = zeros(p);

    for b = 1:num_blocks
        row_start = (b-1)*block_size + 1;
        row_end = b*block_size;
        curr_seed = rand_seed + b;

        % Generate a 5x5 random correlation matrix
        C_block = randCorrMatrix(block_size, curr_seed);

        % Place block into the big matrix
        Corr(row_start:row_end, row_start:row_end) = C_block;
    end

    % Ensure symmetry just in case of numerical issues
    Corr = (Corr + Corr') / 2;
end
