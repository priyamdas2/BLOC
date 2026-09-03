function Theta = Corr2Theta(C)
% C must be a positive definite correlation matrix of size M x M
M = size(C, 1);
L = chol(C, 'lower');

% Preallocate angle vector
Theta = zeros(M*(M-1)/2, 1);
idx = 1;

for m = 2:M
    v = flip(L(m,1:m));              % Flip the m entries of the m-th row
    v = v / norm(v);                  % Normalize to unit norm
    
    if(m == 2)
        angles = atan2(v(2), v(1));
    else
        angles = zeros(1, m-1);
        for k = 1:(m-1)
            if k < m-1
                denom = norm(v(k:end));
                if(denom == 0)
                    angles(k) = 0;
                else
                    angles(k) = acos(v(k) / denom);
                end
            else
                angles(k) = atan2(v(m), v(m-1)); % last angle
                if(angles(k) < 0)
                    angles(k) = angles(k) + 2*pi;
                end
            end
        end
    end
    
    % Store in phi vector
    Theta(idx:idx + m - 2) = angles(:);
    idx = idx + m - 1;
end

end