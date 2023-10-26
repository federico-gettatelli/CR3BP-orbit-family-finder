function Gorbit = find_halo_family(Init_cond, Lib_point_type, Orbit_type)
%% help 
% Finds the Halo orbits, related to a specified lagrangian point,
% through a correction scheme that laverages the x-z symmetry of the orbits
% (for L1,L2,L2).
% for further info: https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Masters/2006_Grebow.pdf
% and: "Three dimensional, periodic, halo, orbits" Kathleen Howell 

% GLOBAL PARAMETERS:
% CHANGE   = change in the x coordinate between two subsequent orbits (use
%            a small value if larget orbits are to be found)
% MAX_ITER = Maximum number of iteration used in the correction scheme 
%            (use a large value if larget orbits are to be found)  
% mu       = M1/(M2+M1) where M1 and M2 are the principal masses of the CR3BP
% R_M2     = Radius of the second body. Used to stop the process.

% INPUT: 
% Lib_point_type = Lagrangian point ('L1','L2','L3','L4','L5') 
% Init_cond      = Initial conditions to initialize the search 
% Orbit_type     = Specifies if Southern ('Hs') or Northern ('Hn') Halo family      

% OUTPUT:
% Gorbit = Nx6 matrix that contains the initial state of N orbits of the
%          Orbit_type family

%%
global mu CHANGE MAX_ITER R_M2
options = odeset('RelTol',1e-11,'AbsTol',1e-11, 'Events', 'cross_xz');
PHI0 = reshape(eye(6),[1,36]);

X0 = Init_cond;
X_SOL = X0;
t_end = 10;

is_finding = true;
is_finding_iter = 0;
norm_z_zn = 2;
norm_x_xn = 1;
Gorbit = [X0];

CHANGE_DZ_AFTER_DX = false;
STOP_CHANGE = false;
DZ = CHANGE;
if strcmp(Orbit_type,'Hs')
    DZ = -DZ;
end

DX = CHANGE;
if strcmp(Lib_point_type,'L2')
    DX = -DX;
end
MIN_Z = 0;
iter_2 = 0;
%%
while is_finding
    vxf = 1;
    vzf = 1;
    x0_old = X_SOL(1);
    z0_old = X_SOL(3);
    iter = 0;
    if norm_z_zn > norm_x_xn
        %%
        X_SOL(3) = X_SOL(3) + DZ;
        while (abs(vxf) > 1e-10 || abs(vzf) > 1e-10) && iter <= MAX_ITER
            [~,S]=ode113(@CR3BP_equations, [0 t_end], [X_SOL PHI0], options);
            % final values
            vxf = S(end,4);
            vzf = S(end,6);
            PHIf = reshape(S(end,7:end),[6,6]);
            
            x = S(end,1);
            y = S(end,2);
            z = S(end,3);
            r1 = sqrt((x+mu)^2+y^2+z^2);
            r2 = sqrt((x+mu-1)^2+y^2+z^2);
            
            dy = S(end,5);
            ddx = 2*dy + x -(1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3;
            ddz = -(1-mu)*z/r1^3 - mu*z/r2^3;
            M = [PHIf(4,1) PHIf(4,5); PHIf(6,1) PHIf(6,5)] - [ddx/dy ddz/dy]'.*[PHIf(2,1) PHIf(2,5)];
            
            % update
            delx0 = M\[vxf, vzf]';
            X_SOL(1) = X_SOL(1) - delx0(1);
            X_SOL(5) = X_SOL(5) - delx0(2);
            MIN_Z = min([abs(X_SOL(3)), abs(z)]);
            iter = iter + 1;
        end        
    elseif norm_x_xn >= norm_z_zn
        %%
        X_SOL(1) = X_SOL(1) + DX;
        if ~STOP_CHANGE
            CHANGE_DZ_AFTER_DX = true;
            STOP_CHANGE = true;
        end
        while (abs(vxf) > 1e-10 || abs(vzf) > 1e-10) && iter <= MAX_ITER
            [~,S]=ode113(@CR3BP_equations, [0 t_end], [X_SOL PHI0], options);
            % final values
            vxf = S(end,4);
            vzf = S(end,6);
            PHIf = reshape(S(end,7:end),[6,6]);
            
            x = S(end,1);
            y = S(end,2);
            z = S(end,3);
            r1 = sqrt((x+mu)^2+y^2+z^2);
            r2 = sqrt((x+mu-1)^2+y^2+z^2);
            
            dy = S(end,5);
            ddx = 2*dy + x -(1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3;
            ddz = -(1-mu)*z/r1^3 - mu*z/r2^3;
            M = [PHIf(4,3) PHIf(4,5); PHIf(6,3) PHIf(6,5)] - [ddx/dy ddz/dy]'.*[PHIf(2,3) PHIf(2,5)];
            
            % update
            delx0 = M\[vxf, vzf]';
            X_SOL(3) = X_SOL(3) - delx0(1);
            X_SOL(5) = X_SOL(5) - delx0(2);
            MIN_Z = min([abs(X_SOL(3)), abs(z)]);
            iter = iter + 1;
        end
    end
    
    if (abs(X_SOL(1)) < 4 && abs(X_SOL(3)) < 4)
        x0_new = X_SOL(1);
        z0_new = X_SOL(3);
        norm_x_xn = abs(x0_new - x0_old);
        norm_z_zn = abs(z0_new - z0_old);
        
        Gorbit = [Gorbit; X_SOL];
        is_finding_iter = 0;
        if CHANGE_DZ_AFTER_DX
            DZ = -DZ;
            CHANGE_DZ_AFTER_DX = false;
        end
        disp(iter)
        iter_2 = iter_2 + 1;
    else
        is_finding_iter = is_finding_iter + 1;
        X_SOL = Gorbit(end,:) + [0 0 2*DZ 0 0 0];
    end
    
    if is_finding_iter == 10
        is_finding = false;
    elseif (MIN_Z < R_M2+CHANGE && iter_2 > 10) || iter_2 > 120
        is_finding = false;
    end
end
end

% T_GUESS = 3.4155;
% options = optimoptions('fsolve','MaxIterations',1000,'MaxFunctionEvaluations',inf);
% z0_FIXED = X_SOL(3) + DZ;
% S_GUESS = [X_SOL(1), X_SOL(5), T_GUESS];
% S_SOL = fsolve(@Shooting_z_fixed, S_GUESS, options);
% X_SOL = [S_SOL(1), 0, z0_FIXED, 0, S_SOL(2), 0];
% T_GUESS = S_SOL(3);
% 
% x0_FIXED = X_SOL(1) + DX;
% S_GUESS = [X_SOL(3), X_SOL(5), T_GUESS];
% S_SOL = fsolve(@Shooting_x_fixed, S_GUESS, options);
% X_SOL = [x0_FIXED, 0, S_SOL(1), 0, S_SOL(2), 0];
% T_GUESS = S_SOL(3);
% 
% function[ERR] = Shooting_x_fixed(S0_GUESS)
% global x0_FIXED
% X0_GUESS = [x0_FIXED, 0, S0_GUESS(1), 0, S0_GUESS(2), 0];
% t_END = S0_GUESS(3);
% 
% options = odeset('RelTol',1e-11,'AbsTol',1e-11);
% [~,S] = ode113(@trecorpi_nonlineare_halo, [0 t_END], X0_GUESS, options);
% 
% ERR = [S(end,2), S(end,4), S(end,6)];
% 
% end
% 
% function[ERR] = Shooting_z_fixed(S0_GUESS)
% global z0_FIXED
% X0_GUESS = [S0_GUESS(1), 0, z0_FIXED, 0, S0_GUESS(2), 0];
% t_END = S0_GUESS(3);
% 
% options = odeset('RelTol',1e-11,'AbsTol',1e-11);
% [~,S] = ode113(@trecorpi_nonlineare_halo, [0 t_END], X0_GUESS, options);
% 
% ERR = [S(end,2), S(end,4), S(end,6)];
% 
% end