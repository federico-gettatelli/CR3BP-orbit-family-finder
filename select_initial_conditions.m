function X0 = select_initial_conditions(Lib_point_type, Orbit_type, Init_condition_type, Init_cond)
%% help
%SCRIVERE

if ~isempty(Init_cond)
    X0 = Init_cond;
    return
end

switch Init_condition_type
    case 'M'
        X0 = init_from_memory(Lib_point_type, Orbit_type);  
    case 'A'
        X0 = init_from_first_order_approximation(Lx, Lib_point_type, Orbit_type);
    case 'R'
        msg = 'Initialization from random not implemented';
        error(msg) 
        %X0 = init_from_random_value(Lx, Ly, Lib_point_type, Orbit_type);
    otherwise
        msg = 'Selected initialization procedure not implemented';
        error(msg)        
end

end
