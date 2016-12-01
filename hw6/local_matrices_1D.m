function [A0, A1] = local_matrices_1D(phi0, phi1,wi,xi)

% [wi, xi] from gauss pts

% init local matrices
A0 = zeros(length(phi0));
A1 = zeros(length(phi1));

% create local mass (u) and stiffness (u') matrices
for i = 1:length(phi0)
    for j = i:length(phi0)
        
        f0 = @(s) phi0{i}(s)*phi0{j}(s);  % function to integrate
        f1 = @(s) phi1{i}(s)*phi1{j}(s);  % derivative function to integrate
        
        %                 % integration of f0 over -1,1 using 4 point quadrature
        %                 A0(i,j) = wi(1)*f0(xi(1)) + wi(2)*f0(xi(2)) + wi(3)*f0(xi(3)) + wi(4)*f0(xi(4));
        %                 % integration of f1 over -1,1 using 4 point quadrature
        %                 A1(i,j) =  wi(1)*f1(xi(1)) + wi(2)*f1(xi(2)) + wi(3)*f1(xi(3)) + wi(4)*f1(xi(4));
        
        for k = 1:length(wi)
            % integration of f0 over -1,1 using 4 point quadrature
            A0(i,j) = A0(i,j) + wi(k)*f0(xi(k));
            % integration of f1 over -1,1 using 4 point quadrature
            A1(i,j) = A1(i,j) + wi(k)*f1(xi(k));
        end
        
    end
end

% copy upper triangular part to lower triangular part for symmetry
A0 = A0 + triu(A0,1)';
A1 = A1 + triu(A1,1)';

end

%         Mass=(h/30)*[4 2 -1; 2 16 2; -1 2 4];
%         Stiff=(1/(3*h))*[7 -8 1; -8 16 -8; 1 -8 7];
%         A0 = Mass; A1 = Stiff;