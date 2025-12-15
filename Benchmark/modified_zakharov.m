function val = modified_zakharov(C)
% MODIFIED_ZAKHAROV_CORR applies Zakharov function to non-diagonal elements of a pd matrix.
% Input:
%   C : N x N pd matrix 
% Output:
%   val : scalar Zakharov function value of non-diagonal elements

    if size(C,1) ~= size(C,2)
        error('Input matrix must be square.');
    end

    % Extract and scale non-diagonal elements
    N = size(C,1);
    non_diag_mask = ~eye(N);
    x = 5 * C(non_diag_mask);  % scale like typical Zakharov domain

    d = length(x);
    indices = 1:d;
    
    sum1 = sum(x.^2);
    sum2 = sum(0.5 * indices * x);

    val = sum1 + sum2^2 + sum2^4;
end