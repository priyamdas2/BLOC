function X = generate_gaussian_samples(sample_size, C)
% Generates 'sample_size' Gaussian observations with mean 0 and correlation matrix C

    % Validate input
    if ~ismatrix(C) || size(C,1) ~= size(C,2)
        error('C must be a square matrix');
    end
    if ~isequal(C, C') || any(eig(C) < 0)
        error('C must be a symmetric positive semi-definite matrix');
    end

    d = size(C, 1);              % Dimension of the multivariate distribution
    mu = zeros(1, d);            % Mean vector (zero)
    
    % Generate samples
    X = mvnrnd(mu, C, sample_size);  % sample_size x d matrix
end