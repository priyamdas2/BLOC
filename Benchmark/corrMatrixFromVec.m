function R = corrMatrixFromVec(x, p)
    R = eye(p);
    idx = tril(true(p), -1);
    R(idx) = x;
    R = R + R' - eye(p);
end