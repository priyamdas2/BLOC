function Corr = generate_corr_Sparsetoeplitz(p, rho, band, out_band_zero_prob, corr_threshold, rand_seed)
%GENERATE_CORR_SPARSETOEPLITZ Generates a sparse correlation matrix
% based on a Toeplitz structure with probabilistic and threshold-based sparsity.
%
% Inputs:
%   p                  : dimension of the correlation matrix
%   rho                : base correlation parameter for Toeplitz structure
%   band               : bandwidth to preserve (|i - j| <= band)
%   out_band_zero_prob : probability of zeroing out elements beyond the band
%   corr_threshold     : threshold to zero out small correlations
%   rand_seed          : random seed for reproducibility
%
% Output:
%   Corr               : p x p sparse correlation matrix

    if nargin < 6
        rand_seed = 1;
    end

    rng(rand_seed);

    % Step 1: Create base Toeplitz matrix
    toeprow = rho .^ (0:(p-1));      % 1, rho, rho^2, ...
    Corr = toeplitz(toeprow);

    % Step 2: Zero out entries with |i-j| > band with some probability
    for i = 1:p
        for j = (i+1):p
            if abs(i - j) > band
                if rand < out_band_zero_prob
                    Corr(i,j) = 0;
                    Corr(j,i) = 0;
                end
            end
        end
    end

    % Step 3: Apply correlation threshold
    Corr(abs(Corr) < corr_threshold) = 0;

    % Step 4: Ensure diagonal is exactly 1
    Corr(1:p+1:end) = 1;

    % Step 5: Symmetrize (for numerical stability)
    Corr = (Corr + Corr') / 2;
end
