function P = covering_Blockdiag(C_true, threshold)
% Constructs a covering matrix P based on C_true.
% Diagonal entries are zero.
% Off-diagonal entries are:
%   1 if abs(C_true(i,j)) <= threshold
%   0 otherwise

    % Initialize P with zeros
    P = zeros(size(C_true));

    % Get logical mask for off-diagonal entries
    off_diag = ~eye(size(C_true));

    % Apply threshold condition to off-diagonal entries
    P(off_diag) = abs(C_true(off_diag)) <= threshold;
end
