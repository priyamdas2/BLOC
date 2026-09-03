function C = generate_corr_SparseUniform(p, sparsity, min_val, max_val, num_tries)
% GENERATE_CORR_SPARSEUNIFORM generates a sparse random correlation matrix.
% Repeats up to num_tries times until the matrix is positive definite (PD).
%
%   p         : number of variables (dimension of the matrix)
%   sparsity  : proportion of off-diagonal entries set to zero (e.g., 0.9)
%   min_val   : minimum value for Uniform sampling of nonzero entries
%   max_val   : maximum value for Uniform sampling
%   num_tries : maximum number of attempts to generate a PD matrix
%
% Returns:
%   C : a p x p PD correlation matrix (or empty if unsuccessful)

    for attempt = 1:num_tries
        % Start with identity matrix
        C = eye(p);

        % Total number of off-diagonal elements in upper triangle
        num_offdiag = p * (p - 1) / 2;
        num_nonzero = round((1 - sparsity) * num_offdiag);

        % Get linear indices of upper triangle (excluding diagonal)
        upper_inds = find(triu(ones(p), 1));
        selected_inds = randsample(upper_inds, num_nonzero);

        % Generate values from Uniform[min_val, max_val]
        values = min_val + (max_val - min_val) * rand(num_nonzero, 1);

        % Fill symmetric off-diagonal values
        for k = 1:num_nonzero
            [i, j] = ind2sub([p, p], selected_inds(k));
            C(i, j) = values(k);
            C(j, i) = values(k);
        end

        % Check if the matrix is positive definite
        [~, flag] = chol(C);
        if flag == 0
            return;  % Successfully generated a PD matrix
        end
    end

    % If none of the attempts succeeded
    warning('Failed to generate a positive definite matrix after %d tries.', num_tries);
    C = [];  % Return empty if no valid matrix was found
end
