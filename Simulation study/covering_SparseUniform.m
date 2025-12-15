function P = covering_SparseUniform(C_true, threshold, rand_seed_here)
% Constructs a covering matrix P based on C_true.
% Diagonal entries are zero.
% Among off-diagonals:
%   - If abs(C_true(i,j)) <= threshold => set to 1
%   - If abs(C_true(i,j)) > threshold => randomly assign 50% to 1, 50% to 0
%
% Inputs:
%   C_true    : square matrix (e.g., correlation matrix)
%   threshold : scalar threshold
%   rand_seed : random seed (optional, for reproducibility)

    if nargin > 2
        rng(rand_seed_here);
    end

    % Initialize
    p = size(C_true, 1);
    P = zeros(p);

    % Logical mask for off-diagonal
    off_diag = ~eye(p);

    % Extract upper triangular off-diagonals to avoid redundancy
    [row, col] = find(triu(off_diag, 1));
    n = length(row);

    for k = 1:n
        i = row(k);
        j = col(k);
        val = abs(C_true(i,j));

        if val <= threshold
            P(i,j) = 1;
            P(j,i) = 1;
        else
            % Randomly assign 1 or 0 with 50% probability
            if rand < 0.5
                P(i,j) = 1;
                P(j,i) = 1;
            end
        end
    end
end
