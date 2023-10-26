function [Lx,Ly] = libration_point_value(Lib_type, mu)

switch Lib_type
    case 'L1'
        p = [1, (mu-3), (3-2*mu), -mu, 2*mu, -mu];
        r = roots(p);
        Lx = 1-mu-r(imag(r)==0); 
        Ly = 0;
    case 'L2'
        p = [1, -(mu-3), (3-2*mu), -mu, -2*mu, -mu];
        r = roots(p);
        Lx = 1-mu+r(imag(r)==0); 
        Ly = 0;
    case 'L3'
        %(L3 formula could be wrong)
        p = [1, (2-4*mu), (6*mu^2-6*mu+1), -4*mu^3+6*mu^2-2*mu-1, mu^4-2*mu^3+mu^2+4*mu-2, -3*mu^2+3*mu-1];
        r = roots(p);
        Lx = -r(imag(r)==0); 
        Ly = 0;
    case 'L4'
        Lx = 1/2 - mu;
        Ly = sqrt(3)/2;
    case 'L5'
        Lx = 1/2 - mu;
        Ly = -sqrt(3)/2;
    otherwise
        msg = 'Selected incorrect libration point';
        error(msg)
end

end