function val = map_start_elements(angle)
% Computes f(phi) = pi/2 - |mod(phi, pi) - pi/2|

    val = pi/2 - abs(mod(angle, pi) - pi/2);
end