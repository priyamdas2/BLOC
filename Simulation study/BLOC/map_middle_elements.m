function val = map_middle_elements(angle)
% Computes f(phi) = pi - |mod(phi, 2*pi) - pi|

    val = pi - abs(mod(angle, 2*pi) - pi);
end