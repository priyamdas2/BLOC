function [c, ceq] = correlationConstraint(x)
    p = (1 + sqrt(1 + 8*length(x)))/2; % Solve for p from length(x)
    R = corrMatrixFromVec(x, round(p));
    eigvals = eig(R);
    c = -eigvals + 1e-8; % ensure all eigenvalues >= 0
    ceq = []; % no equality constraints
end