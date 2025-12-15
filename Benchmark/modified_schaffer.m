function val = modified_schaffer(C)
% MODIFIED_SCHAFFER_CORR applies Schaffer N.2 function to non-diagonal elements of a pd matrix.
% Input:
%   C : N x N pd matrix 
% Output:
%   val : scalar Schaffer N.2 value of non-diagonal elements

    if size(C,1) ~= size(C,2)
        error('Input matrix must be square.');
    end

    % Extract and scale non-diagonal elements
    N = size(C,1);
    non_diag_mask = ~eye(N);
    x = 10 * C(non_diag_mask);  % scale to [-10, 10] like usual domain

    % Apply Schaffer N.2 to adjacent pairs
    val = 0;
    for i = 1:length(x)-1
        xi = x(i);
        xj = x(i+1);
        num = sin(xi^2 - xj^2)^2 - 0.5;
        denom = (1 + 0.001*(xi^2 + xj^2))^2;
        val = val + 0.5 + num / denom;
    end
end
