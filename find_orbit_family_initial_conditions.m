function Gorbit = find_orbit_family_initial_conditions(Lib_point_type, Orbit_type, Init_cond, Lx)
%% help
% Wrapper of the functions used to find the different family of orbits 
% around the lagrangian points of the CR3BP 

% INPUT:
% Lib_point_type = Lagrangian point ('L1','L2','L3','L4','L5') 
% Orbit_type     = Which orbit family to find 
% Init_cond      = Initial conditions to initialize the search 
% Lx             = x coordinate of the Lib_point_type

% OUTPUT:
% Gorbit = Nx6 matrix that contains the initial state of N orbits of the
%          Orbit_type family

switch Orbit_type
    case {'Hn','Hs'}
        Gorbit = find_halo_family(Init_cond, Lib_point_type, Orbit_type);
    case 'L'
        Gorbit = find_lyapunov_family(Init_cond, Lib_point_type, Lx);
    otherwise
        msg = 'Orbit Family search algorithm not implemented';
        error(msg)
end