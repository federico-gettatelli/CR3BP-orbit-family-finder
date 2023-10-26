%% help
% Finds the orbit family related to a lagrangian point.
% Lib_point_type specifies the lagrangian point ('L1','L2','L3','L4','L5')
% Orbit_type specifies the type of orbit family to be found: L = Lyapunov/Planar, 
%   Hs = Southern Halo, Hn = Northern halo, V = Vertical; A = Axial; B = Butterfly 
% Init_condition_type specifies from which orbit to initialize the search
% If "Init_cond" is empty initial conditions are given by Init_condition_type
% Most parameters of the code are tuned for the Earth-Moon system

clc; clear all; close all; 
%% Choose problem
global dconv mu CHANGE MAX_ITER R_M2
Lib_point_type = 'L1';       
Orbit_type = 'Hs';            
Init_condition_type = 'M';    
Init_cond = [];             
MAX_ITER = 30;

%% Variables 
% Mass [kg]: Earth = 5.9724E24; Moon = 7.346E22;
M1 = 5.9724E24; 
M2 = 7.346E22; 
mu = M2/(M1+M2);

G = 6.6743*10^-11;
dconv = 338400;
vconv = sqrt(G*(M1+M2)/(dconv*10^9));
tconv = dconv/vconv;
CHANGE = 1000 / dconv; % Use low values for lyapunov orbits
R_M2 = 1740 / dconv;

[Lx,Ly] = libration_point_value(Lib_point_type, mu);

%% Initial Conditions
Init_cond = select_initial_conditions(Lib_point_type, Orbit_type, Init_condition_type, Init_cond);

%% Find Orbit Family
Orbit_Family = find_orbit_family_initial_conditions(Lib_point_type, Orbit_type, Init_cond, Lx);

%% Post Process
figure(3)
hold on

for i=1:length(Orbit_Family(:,1))-1
    options = odeset('RelTol',1e-11,'AbsTol',1e-11);
    [~,X] = ode45(@trecorpi_nonlineare_halo, [0 7], Orbit_Family(i,:), options);
    plot3(X(:,1), X(:,2), X(:,3),'b-')
end

grid on
plot3(1-mu, 0, 0, 'or', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
plot3(Lx, Ly, 0,'or', 'MarkerSize', 6, 'MarkerFaceColor', 'b')





