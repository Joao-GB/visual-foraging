function delta = angle_diff(a2, a1)
    delta = a2 - a1;
    delta = mod(delta + pi, 2*pi) - pi;  % wrap to [-pi, pi]
end
