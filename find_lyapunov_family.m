function [Gorbit] = find_lyapunov_family(Init_cond, Lib_point_type, Lx)
%% help
% Finds the Lyapunov orbits, related to a specified lagrangian point,
% through a correction scheme that laverages the x-z symmetry of the orbits
% (for L1,L2,L2). L4 and L5 search scheme not implememted yet.
% for further info: https://engineering.purdue.edu/people/kathleen.howell.1/Publications/Masters/2006_Grebow.pdf

% GLOBAL PARAMETERS:
% CHANGE   = change in the x coordinate between two subsequent orbits (use
%            a small value if larget orbits are to be found)
% MAX_ITER = Maximum number of iteration used in the correction scheme 
%            (use a large value if larget orbits are to be found)  
% mu       = M1/(M2+M1) where M1 and M2 are the principal masses of the CR3BP

% INPUT: 
% Lib_point_type = Lagrangian point ('L1','L2','L3','L4','L5') 
% Init_cond      = Initial conditions to initialize the search 
% Lx             = x coordinate of the Lib_point_type

% OUTPUT:
% Gorbit = Nx6 matrix that contains the initial state of N orbits of the
%          Orbit_type family

%%
global CHANGE MAX_ITER mu 
options = odeset('RelTol',1e-11,'AbsTol',1e-11, 'Events', 'cross_xz');
PHI0 = reshape(eye(6),[1,36]);
X0 = Init_cond;
T_END = 5;

DX = CHANGE;
if strcmp(Lib_point_type, 'L2')
    DX = -DX;
end

TOL_ERR = 10^-6;
if strcmp(Lib_point_type, 'L3')
    TOL_ERR = 10^-3;
end

is_finding = true;
CHANGE_DX = false;

X_SOL = X0 - [DX 0 0 0 0 0];
Gorbit = [];
is_finding_iter = 0;
i_good = 0;
LEFT_X = 0;
RIGHT_X = 0;

while is_finding
    vxf = 1;
    iter = 0;    
    X_SOL(1) = X_SOL(1) + DX; 
    
    while (abs(vxf) > TOL_ERR && iter <= MAX_ITER)
        [~,S]=ode45(@CR3BP_equations, [0 T_END], [X_SOL PHI0], options);
        
        % final values
        vxf = S(end,4);
        PHIf = reshape(S(end,7:end),[6,6]);
        
        x = S(end,1);
        y = S(end,2);
        z = S(end,3);
        r1 = sqrt((x+mu)^2+y^2+z^2);
        r2 = sqrt((x+mu-1)^2+y^2+z^2);
        dy = S(end,5);
        ddy = -2*vxf + y -(1-mu)*y/r1^3 - mu*y/r2^3;
        
        deldy0 = vxf / (PHIf(4,5)-PHIf(2,5)*ddy/dy);
        
        X_SOL(5) = X_SOL(5) - deldy0;
        iter = iter + 1;
        RIGHT_X = max(X_SOL(1), x);
        LEFT_X = min(X_SOL(1), x);
    end
    
    if iter <= MAX_ITER
        i_good = i_good + 1;
        if mod(i_good,1) == 0
            Gorbit = [Gorbit; X_SOL];
        end
    else
        is_finding_iter = is_finding_iter + 1;
    end
    disp(iter);
    
    if ~CHANGE_DX
        switch Lib_point_type
            case {'L1','L2'}
                if min([abs(LEFT_X-Lx), abs(RIGHT_X-Lx), abs(LEFT_X-1+mu), abs(RIGHT_X-1+mu)]) < CHANGE
                    CHANGE_DX = true;
                    DX = -DX;
                    X_SOL = X0;
                end
            case 'L3'
                if min([abs(LEFT_X-Lx), abs(RIGHT_X-Lx), abs(RIGHT_X+mu)]) < CHANGE
                    CHANGE_DX = true;
                    DX = -DX;
                    X_SOL = X0;
                end
            otherwise
                msg = 'L4 and L5 lyapunov search not implemented yet';
                error(msg)
        end
    end
    
    if is_finding_iter == 5
        is_finding = false;                       
    end
end
end