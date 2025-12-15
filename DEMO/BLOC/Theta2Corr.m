function C = Theta2Corr(Theta)
% Compute M from the length of Theta
len_Theta = length(Theta);
M = (1 + sqrt(1 + 8*len_Theta)) / 2;
if mod(M, 1) ~= 0
    error('Invalid length of Theta. It must be M*(M-1)/2 for some integer M.');
end
M = round(M);

% Initialize L
L = zeros(M, M);
idx = 1;

for m = 1:M
    if m == 1
        L(1,1) = 1;
    else
        theta_slice = Theta(idx : idx + m - 2);
        idx = idx + m - 1;
        
        row = zeros(1, m);
        for k = 1:m
            if k == 1
                row(k) = cos(theta_slice(1));
            elseif k < m
                prod_sin = prod(sin(theta_slice(1:k-1)));
                row(k) = prod_sin * cos(theta_slice(k));
            else
                row(k) = prod(sin(theta_slice));
            end
        end
        L(m,1:m) = flip(row);
    end
end

if min(L(1:M+1:end)) < 10^-3
    warning('C may be close to rank deficient (small diagonal in L).');
end
% Compute correlation matrix
C = L * L';
end