function val = modified_rosenbrock(C)
% MODIFIED_ROSENBROCK_CORR applies Rosenbrock function to non-diagonal elements of a pd matrix.
% Input:
%   C : N x N pd matrix 
% Output:
%   val : scalar Rosenbrock function value of non-diagonal elements

    if size(C,1) ~= size(C,2)
        error('Input matrix must be square.');
    end


    % Extract non-diagonal elements
    N = size(C,1);
    mask = ~eye(N);
    x = 100*C(mask);  % Vector of length d = N^2 - N
    d = numel(x);

    % Modified Rosenbrock formula
    val = 0;
    for i = 1:(d-1)
        % Compute Rosenbrock for non-diagonal elements
        val = val + 100 * (x(i+1) - x(i)^2)^2 + (x(i)-1)^2;
    end
end