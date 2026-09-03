function X = cov_model_matDepthsim(n, model)
% cov_model_matDepthsim
%
% Covariance models for the matrix-depth simulation.
%
% model = 1 : Block matrix with fixed block size = 5
% model = 2 : Toeplitz matrix
% model = 3 : Banded matrix

X = zeros(n);

if model == 1 % Block matrix

    bs = 5;       % block size
    bn = n/bs;    % number of blocks

    if mod(n,bs) ~= 0
        error('For the block model, n must be divisible by block size 5.');
    end

    for ib = 1:bn

        idx = (1:bs) + (ib-1)*bs;

        X(idx,idx) = 0.8*ones(bs,bs);

    end

    X = X + 0.2*eye(n);


elseif model == 2 % Toeplitz matrix

    for k = 1:n
        for l = 1:n
            X(k,l) = 0.75^abs(k-l);
        end
    end


elseif model == 3 % Banded matrix

    for k = 1:n
        for l = 1:n
            X(k,l) = (1-abs(k-l)/10) * ...
                     (abs(k-l)<=10);
        end
    end

end

end