function D_avg = matrix_depth_fixedU_avg(C, proj2, U, beta)

CU = C * U;

directional_scale = sum(U .* CU, 1);

threshold = beta * directional_scale;

prop_below = mean(proj2 <= threshold, 1);
prop_above = mean(proj2 >= threshold, 1);

direction_depth = min(prop_below, prop_above);

D_avg = mean(direction_depth);

end