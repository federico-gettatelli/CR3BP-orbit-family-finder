function print_figure_orbits(Orbit_Family, STABILITY_INDEX, Lx, Ly, Lib_point_type, NRHO, STABILITY_THRESHOLD)
global mu dconv

ORBITS_TO_PRINT = Orbit_Family;
VALID_INDEX = [];

ITER_NRHO = 1;
if strcmp(NRHO, 'true')    
    ORBITS_TO_PRINT = [];
    for i=1:length(STABILITY_INDEX(:,1)) 
        if abs(STABILITY_INDEX(i,1:6)) < STABILITY_THRESHOLD
            ORBITS_TO_PRINT(ITER_NRHO, :) = Orbit_Family(i,:);
            VALID_INDEX(ITER_NRHO, :) = STABILITY_INDEX(i,:);
            ITER_NRHO = ITER_NRHO + 1;
        end
    end
    figure(2)
    title('Stability Indices')
    hold on
    plot(STABILITY_INDEX(:,7)*dconv, STABILITY_INDEX(:,1:6), 'k*')
    plot(VALID_INDEX(:,7)*dconv, VALID_INDEX(:,1:6), 'b*')
    grid on
    axis([0 2*10^4 -3 3])
end

figure(1)
title('Halo Orbits')
hold on
for i=1:length(ORBITS_TO_PRINT(:,1))    
    options = odeset('RelTol',1e-11,'AbsTol',1e-11);
    [~,X] = ode45(@non_linear_3body, [0 7], ORBITS_TO_PRINT(i,:), options);
    plot3(X(:,1), X(:,2), X(:,3),'b-')
end
grid on
plot3(1-mu, 0, 0, 'or', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
plot3(Lx, Ly, 0,'or', 'MarkerSize', 6, 'MarkerFaceColor', 'b')
text(1-mu, 0, 0,'Moon')
text(Lx, Ly, 0,Lib_point_type)

save('Orbits_initial_conditions.mat', 'ORBITS_TO_PRINT')
writematrix(ORBITS_TO_PRINT, 'Orbits_initial_conditions.txt')
end