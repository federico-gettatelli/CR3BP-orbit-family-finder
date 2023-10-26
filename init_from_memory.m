function X0 = init_from_memory(Lib_point_type, Orbit_type)
%% help
% Loads the state of an orbit of the desired orbit family.
% Values valid for the Earth-Moon CR3BP.

% INPUT: 
% Lib_point_type = Lagrangian point ('L1','L2','L3','L4','L5') 
% Orbit_type     = Specifies orbit family type 

% OUTPUT:
% X0      = Initial state to initialize the search 

%%
switch Orbit_type
    case 'L'
        switch Lib_point_type
            case 'L1'
                X0 = [0.8234, 0, 0, 0, 0.1262379, 0];
            case 'L2'
                X0 = [1.1762, 0, 0, 0, -0.1228570, 0];
            case 'L3'
                X0 = [-1.6967, 0, 0, 0, 1.2795, 0];%[-1.0560, 0, 0, 0, 0.1017, 0];
            case {'L4','L5'}
                X0 = [0.4750, sqrt(3)/2, 0, 0.0697, -1.0915, 0];
        end
    case 'V'
        switch Lib_point_type
            case 'L1'
                X0 = [1.0118, 0, 0.1739, 0, -0.0799, 0];
            case 'L2'
                X0 = [1.1119, 0, 0, 0, -0.1812 0.4358];
            case 'L3'
                X0 = [-0.7916, 0, 0.6160, 0, -0.2129, 0];
            case {'L4','L5'}
                X0 = [0.6032, 0.42545, 0.6678, -0.3544, -0.3857, 0.5728];
        end
    case 'A'
        switch Lib_point_type
            case 'L1'
                X0 = [0.7816, 0, 0, 0, 0.4432, 0.0000];
            case 'L2'
                X0 = [1.2200, 0, 0, 0, -0.4275, 0.0000];
            case 'L3'
                X0 = [-1.8963, 0, 0, 0, 1.6715 0.0000];
            case {'L4','L5'}
                X0 = [0.9767, -0.0334, 0.1, -0.0737, -0.0626, 0.4671];
        end
    case 'B'
        switch Lib_point_type
            case 'L1'
                X0 = [];
            case 'L2'
                X0 = [1.0118, 0, 0.1739, 0, -0.0799, 0];
            case 'L3'
                X0 = [];
            case {'L4','L5'}
                X0 = [];
        end
    case {'Hs', 'Hn'}
        switch Lib_point_type
            case 'L1'
                X0 = [0.8234, 0, 0.0000, 0, 0.1262379, 0];
            case 'L2'
                X0 = [1.180928, 0, 0.0000, 0, -0.1560761, 0];
            case 'L3'
                X0 = [-1.6967, 0, 0.0000, 0, 1.2796, 0];
        end
end

if isempty(X0)
    msg = 'Orbit not in memory';
    error(msg)
end

end