function x = vecFromCorrMatrix(R)
    p = size(R, 1);
    idx = tril(true(p), -1);
    x = R(idx);
end