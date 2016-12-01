function [phi0, phi1] = basis_functions_1D(deg)

% get basis functions from polynomial interpolants

phi0 = {};
phi1 = {};

% get shape basis functions for element, Lagrange interpolation (-1,1)
% P =  degree interpolating polynomial
switch deg
    case 1  % linear
        phi0{1} = @(s) -0.5*s + 0.5;
        phi0{2} = @(s) 0.5*s + 0.5;
        phi1{1} = @(s) -0.5;  % derivative
        phi1{2} = @(s) 0.5;
    case 2  % quadratic
        phi0{1} = @(s) 0.5*s^2 - 0.5*s;
        phi0{2} = @(s) -s^2 + 1;
        phi0{3} = @(s) 0.5*s^2 + 0.5*s;
        phi1{1} = @(s) s - 0.5;
        phi1{2} = @(s) -2*s;
        phi1{3} = @(s) s + 0.5;
    case 3  % cubic
        phi0{1} = @(s) -27/48*(s^3 - s^2 - 1/9*s + 1/9);
        phi0{2} = @(s) 27/16*(s^3 - 1/3*s^2 - s + 1/3);
        phi0{3} = @(s) -27/16*(s^3 + 1/3*s^2 - s - 1/3);
        phi0{4} = @(s) 27/48*(s^3 + s^2 - 1/9*s - 1/9);
        phi1{1} = @(s) -27/48*(3*s^2 - 2*s - 1/9);
        phi1{2} = @(s) 27/16*(3*s^2 - 2/3*s - 1);
        phi1{3} = @(s) -27/16*(3*s^2 + 2/3*s - 1);
        phi1{4} = @(s) 27/48*(3*s^2 + 2*s - 1/9);
end