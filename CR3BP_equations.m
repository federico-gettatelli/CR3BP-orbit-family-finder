function dS = CR3BP_equations(t,S0)
%% help
% Equations of motion of an object in the CR3BP + Equations of the state transition matrix associated with it 

global mu
dS = zeros(42,1);

x = S0(1);
y = S0(2);
z = S0(3);
vx = S0(4);
vy = S0(5);
vz = S0(6);

r1 = sqrt((x+mu)^2+y^2+z^2);
r2 = sqrt((x+mu-1)^2+y^2+z^2);

PHI0 = reshape(S0(7:42),[6,6]);

dS(1) = vx;
dS(2) = vy;
dS(3) = vz;
dS(4) = 2*vy + x -(1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3;
dS(5) = -2*vx + y -(1-mu)*y/r1^3 - mu*y/r2^3;
dS(6) = -(1-mu)*z/r1^3 - mu*z/r2^3;

U11 = 1 -(1-mu)/r1^3 - mu/r2^3 + 3*(1-mu)*(x+mu)^2/r1^5 + 3*mu*(x+mu-1)^2/r2^5;
U12 = + 3*(1-mu)*(x+mu)*y/r1^5 + 3*mu*(x+mu-1)*y/r2^5;
U13 = + 3*(1-mu)*(x+mu)*z/r1^5 + 3*mu*(x+mu-1)*z/r2^5;

U22 = 1 -(1-mu)/r1^3 - mu/r2^3 + 3*(1-mu)*y^2/r1^5 + 3*mu*y^2/r2^5;
U23 = + 3*(1-mu)*y*z/r1^5 + 3*mu*y*z/r2^5;

U33 = -(1-mu)/r1^3 - mu/r2^3 + 3*(1-mu)*z^2/r1^5 + 3*mu*z^2/r2^5;

U = [U11 U12 U13; U12 U22 U23; U13 U23 U33];
OMEGA = [0 1 0; -1 0 0; 0 0 0];
F = [zeros(3,3) eye(3); U 2*OMEGA];

dPHI = F*PHI0;

dS(7:42) = reshape(dPHI, [36,1]);
end
