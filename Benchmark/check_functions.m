clear all;
M = 10;
no_checks = 100;
norm_diffs = nan(no_checks,1);
is_valid = nan(no_checks,1);
for jj = 1:no_checks
    C0 = randCorrMatrix(M, round(rand(1)*10000));
    Theta = Corr2Theta(C0);
    N = length(Theta);
    [domain_lb, domain_ub] = Domain_ub_lb(N);
    is_valid(jj) = all(Theta >= domain_lb & Theta <= domain_ub);
    C_est = Theta2Corr(Theta);
    norm_diffs(jj) = norm(C0 - C_est);
end

max(norm_diffs)
[min(is_valid), max(is_valid)]

M = 10;
C0 = randCorrMatrix(M, round(rand(1)*10000));
Theta0 = Corr2Theta(C0);
N = length(Theta0);

start_element_locs = ThetaGroupStartIndicesRow2onward(M)
end_element_locs = ThetaGroupEndIndicesRow2onward(M)


