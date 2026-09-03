function C = spd_to_corr(X)
% ================================================================
% SPD_TO_CORR
%
% Convert an SPD matrix X to its corresponding correlation matrix:
%
%       C = D^(-1/2) X D^(-1/2),
%
% where
%
%       D = diag(diag(X)).
%
% The resulting C has unit diagonal.
% ================================================================

% Numerical symmetrization
X = (X+X')/2;

d = diag(X);

if any(~isfinite(d))

    error( ...
        'Nonfinite diagonal encountered in SPD matrix.');

end

if any(d <= 0)

    error( ...
        'Nonpositive diagonal encountered in SPD matrix.');

end

dinv = 1 ./ sqrt(d);

% Equivalent to D^(-1/2) X D^(-1/2)
C = X .* (dinv*dinv');

% Numerical cleanup
C = (C+C')/2;

p = size(C,1);

C(1:p+1:end) = 1;

end