function X0 = init_from_first_order_approximation(Lx, Lib_point_type, Orbit_type)
%% help
% Computes the initial values of a "linear" lyapunov orbit.
% NOT TESTED

global mu dconv

if ~strcmp(Orbit_type,'L')
    msg = 'Initialization from approximation implemented only for lyapunov';
    error(msg) 
end

Ax = 10000 / dconv;

d1 = abs(Lx+mu); 
d2 = abs(Lx-1+mu); 
Uxx = 1 - (1-mu)/(d1)^3 - mu/(d2)^3 + 3*(1-mu)*(d1)^2/(d1)^5 + 3*mu*(d2)^2/(d2)^5;
Uyy = 1 - (1-mu)/(d1)^3 - mu/(d2)^3; 
    
b1 = 2 - (Uxx+Uyy)/2;
b2_2 = -Uxx*Uyy;
s = sqrt(b1+sqrt(b1^2+b2_2));
b3 = (s^2 - Uxx) / (2*s); % could be s^2+Uxx / 2s

x0 = -Ax+Lx;
y0 = 0;
z0 = 0;
dx0 = 0;
dy0 = b3*Ax*s;
dz0 = 0;
X0 = [x0, y0, z0, dx0, dy0, dz0];
end