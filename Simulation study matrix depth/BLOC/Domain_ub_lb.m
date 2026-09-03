function [domain_lb, domain_ub] = Domain_ub_lb(N)
M = (1 + sqrt(1 + 8*N)) / 2;
start_element_locs = ThetaGroupStartIndicesRow2onward(M);
end_element_locs = ThetaGroupEndIndicesRow2onward(M);

domain_lb = zeros(N,1);
domain_lb(1) = -pi/2;
domain_ub_temp = pi*ones(N,1);
domain_ub_temp(end_element_locs) = 2*pi;
domain_ub_temp(start_element_locs) = pi/2;
domain_ub = domain_ub_temp;
domain_ub(1) = pi/2;
end