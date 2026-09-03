function indices = ThetaGroupStartIndicesRow2onward(M)
    % Returns the starting indices of each angle group in the Theta vector
    % (Row 2 onwards)
    % M: dimension of the correlation matrix
    
    indices = zeros(1,M-1);
    for m = 1:(M-1)
        indices(m) = (m-1)*m/2 + 1;
    end
    indices = indices(2:end);
end