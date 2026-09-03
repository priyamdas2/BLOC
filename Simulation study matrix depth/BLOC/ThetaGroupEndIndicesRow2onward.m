function indices = ThetaGroupEndIndicesRow2onward(M)
    % Returns the ending indices of each angle group in Theta,
    % starting from row 2 (i.e., skipping the first element)
    
    indices = zeros(1, M-2);
    for m = 2:(M-1)
        indices(m-1) = m*(m+1)/2;
    end
end