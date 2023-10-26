function [value, isterminal, direction] = cross_xz(t,S)
    value = (S(2) >= 0);
    isterminal = 1;
    direction = 0;
end