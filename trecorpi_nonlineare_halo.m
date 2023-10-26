function dS = trecorpi_nonlineare_halo(t,S0)
    global mu
    dS = zeros(6,1);
    x = S0(1);
    y = S0(2);
    z = S0(3);
    vx = S0(4);
    vy = S0(5);
    vz = S0(6);

    r1 = sqrt((x+mu)^2+y^2+z^2);
    r2 = sqrt((x+mu-1)^2+y^2+z^2);

    dS(1) = vx;
    dS(2) = vy;
    dS(3) = vz;
    dS(4) = 2*vy + x -(1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3;
    dS(5) = -2*vx + y -(1-mu)*y/r1^3 - mu*y/r2^3;
    dS(6) = -(1-mu)*z/r1^3 - mu*z/r2^3;
end