function x_mapped = smooth_reflect_pi2(x)
    % Smooth reflection into (-pi/2, pi/2) using arcsin(sin(x))
    x_mapped = asin(sin(x));
end