% Numeric PDEs
% Programming Assignment #2
% Sam Friedman
% 9/27/2016


% ref: http://www.math.chalmers.se/~mohammad/teaching/PDEbok/Draft_I+II.pdf

poly = [1:3]';  %#ok<NBRAK> % polynomial degrees

% approx u(x) ~ U(x)
% Prob 1.1
data_dir = 'data/';
hw = 2;
prob = 1;
part = 1;
prbsfx = [data_dir,'hw',num2str(hw),'pb',num2str(prob),'pt',num2str(part)];

gauss_pts = 4;  % use 4 quadrature points for all polynomial basis functions

allA = cell(1,length(poly));
allAprm = cell(1,length(poly));

% Gauss Quadrature weights and abscissa
% ref: https://pomax.github.io/bezierinfo/legendre-gauss.html
switch gauss_pts
    case 2
        wi = [ 1.0000000000000000, 1.0000000000000000];
        xi = [-0.5773502691896257, 0.5773502691896257];
    case 3
        wi = [ 0.8888888888888888, 0.5555555555555556, 0.5555555555555556];
        xi = [ 0.0000000000000000,-0.7745966692414834, 0.7745966692414834];
    case 4
        wi = [ 0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538];
        xi = [-0.3399810435848563, 0.3399810435848563,-0.8611363115940526, 0.8611363115940526];
    case 5
        wi = [ 0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891];
        xi = [ 0.0000000000000000,-0.5384693101056831, 0.5384693101056831,-0.9061798459386640, 0.9061798459386640];
    case 6
        wi = [ 0.3607615730481386, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.1713244923791704, 0.1713244923791704];
        xi = [ 0.6612093864662645,-0.6612093864662645,-0.2386191860831969, 0.2386191860831969,-0.9324695142031521, 0.9324695142031521];
end

% nrms = table(poly);
% nrms.Properties.VariableNames = {'ply'};

for p = 1:length(poly);
    
    ply = poly(p);  % degree polynomial
    
    phi = {};
    phiprm = {};
    
    % get shape basis functions for element, Lagrange interpolation (-1,1)
    switch ply
        case 1  % linear
            phi{1} = @(s) -0.5*s + 0.5;
            phi{2} = @(s) 0.5*s + 0.5;
            phiprm{1} = @(s) -0.5;  % derivative
            phiprm{2} = @(s) 0.5;        
        case 2  % quadratic
%             phi{1} = @(s) (2*s - 1)*(s - 1);
%             phi{2} = @(s) 4*s*(s - 1);
%             phi{3} = @(s) s*(2*s - 1);
%             phiprm{1} = @(s) 4*s - 3;
%             phiprm{2} = @(s) -8*s + 4;
%             phiprm{3} = @(s) 4*s - 1;
            phi{1} = @(s) 0.5*s^2 - 0.5*s;
            phi{2} = @(s) -s^2 + 1;
            phi{3} = @(s) 0.5*s^2 + 0.5*s;
            phiprm{1} = @(s) s - 0.5;
            phiprm{2} = @(s) -2*s;
            phiprm{3} = @(s) s + 0.5;
        case 3  % cubic
            phi{1} = @(s) -27/48*(s^3 - s^2 - 1/9*s + 1/9);
            phi{2} = @(s) 1/3*(s^3 - 1/3*s^2 - s + 1/3);
            phi{3} = @(s) -1/3*(s^3 + 1/3*s^2 - s - 1/3);
            phi{4} = @(s) 27/48*(s^3 + s^2 - 1/9*s - 1/9);
            phiprm{1} = @(s) -27/48*(3*s^2 - 2*s - 1/9);
            phiprm{2} = @(s) 1/3*(3*s^2 - 2/3*s - 1);
            phiprm{3} = @(s) -1/3*(3*s^2 + 2/3*s - 1);
            phiprm{4} = @(s) 27/48*(3*s^2 + 2*s - 1/9);
    end
    
    jac = 1;  % jacobian
    
    A = zeros(length(phi));
    Aprm = zeros(length(phiprm));
    
    for i = 1:length(phi)
        for j = i:length(phi)
            
            f = @(s) phi{i}(s)*phi{j}(s);  % function to integrate
            % integration of f over -1,1 using 4 point quadrature
            A(i,j) = wi(1)*f(xi(1)) + wi(2)*f(xi(2)) + wi(3)*f(xi(3)) + wi(4)*f(xi(4));  
            
            fprm = @(s) phiprm{i}(s)*phiprm{j}(s);  % derivative function to integrate
            % integration of fprm over -1,1 using 4 point quadrature
            Aprm(i,j) =  wi(1)*fprm(xi(1)) + wi(2)*fprm(xi(2)) + wi(3)*fprm(xi(3)) + wi(4)*fprm(xi(4)); 
        
        end
    end
    
   % copy upper triangular part to lower triangular part for symmetry
   A = A + triu(A,1)';  
   Aprm = Aprm + triu(Aprm,1)';
   
   allA{p} = A;
   allAprm{p} = Aprm;
%    soln = table(x);  % init soln data table with FD points
    
%     % save solution data to table
%     solnsfix = ['ply',num2str(ply)];
%     solnvars = {['uext', solnsfix],['uapx', solnsfix],['e', solnsfix]};
%     tmp = table(uext,uapx,e);
%     tmp.Properties.VariableNames = solnvars;
%     soln = [soln, tmp]; %#ok<AGROW> % concatenate tmp table with solutions table
%     
end
%
%     % write solutions table to data file for each n value
%     writetable(soln,[prbsfx,'n',num2str(n),'soln','.dat'],'Delimiter','\t');
%
% end
%
% % export data to .dat file
%
% nrms = zeros(size(L2,1),size(L2,2)*2+1);
% nrmsvars = {'n'};
% nrms(:,1) = nn;
%
% for j = 1:size(L2,2)
%
%     nrms(:,j*2) = L2(:,j);
%     nrms(:,j*2+1) = Linf(:,j);
%
%     beta = betalist(j);
%     nrmssfix = ['be',num2str(beta)];
%     nrmsvars = [nrmsvars, {['L2', nrmssfix], ['Linf', nrmssfix]}]; %#ok<AGROW>
%
% end
%
% nrmtbl = array2table(nrms);
% nrmtbl.Properties.VariableNames = nrmsvars;
% writetable(nrmtbl,[prbsfx,'nrms','.dat'],'Delimiter','\t');

% s = sym(L);
% v = vpa(s,5); % assign numerical precision
% latex(v)




