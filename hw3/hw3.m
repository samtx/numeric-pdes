% Numeric PDEs
% Programming Assignment #2
% Sam Friedman
% 9/27/2016


% ref: http://www.math.chalmers.se/~mohammad/teaching/PDEbok/Draft_I+II.pdf

poly = [1:3]';  %#ok<NBRAK> % polynomial degrees
N = 10; % number of elements

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
    
    P = poly(p);  % degree polynomial
    
    phi0 = {};
    phi1 = {};
    phi2 = {};
    
    % get shape basis functions for element, Lagrange interpolation (-1,1)
    switch P
        case 1  % linear
            phi0{1} = @(s) -0.5*s + 0.5;
            phi0{2} = @(s) 0.5*s + 0.5;
            phi1{1} = @(s) -0.5;  % derivative
            phi1{2} = @(s) 0.5;
            phi2{1} = @(s) 0;
            phi2{2} = @(s) 0;
        case 2  % quadratic
            phi0{1} = @(s) 0.5*s^2 - 0.5*s;
            phi0{2} = @(s) -s^2 + 1;
            phi0{3} = @(s) 0.5*s^2 + 0.5*s;
            phi1{1} = @(s) s - 0.5;
            phi1{2} = @(s) -2*s;
            phi1{3} = @(s) s + 0.5;
            phi2{1} = @(s) 1;
            phi2{2} = @(s) -2;
            phi2{3} = @(s) 1;
        case 3  % cubic
            phi0{1} = @(s) -27/48*(s^3 - s^2 - 1/9*s + 1/9);
            phi0{2} = @(s) 1/3*(s^3 - 1/3*s^2 - s + 1/3);
            phi0{3} = @(s) -1/3*(s^3 + 1/3*s^2 - s - 1/3);
            phi0{4} = @(s) 27/48*(s^3 + s^2 - 1/9*s - 1/9);
            phi1{1} = @(s) -27/48*(3*s^2 - 2*s - 1/9);
            phi1{2} = @(s) 1/3*(3*s^2 - 2/3*s - 1);
            phi1{3} = @(s) -1/3*(3*s^2 + 2/3*s - 1);
            phi1{4} = @(s) 27/48*(3*s^2 + 2*s - 1/9);
            phi2{1} = @(s) -27/48*(6*s - 2);
            phi2{2} = @(s) 1/3*(6*s - 2/3);
            phi2{3} = @(s) -1/3*(6*s + 2/3);
            phi2{4} = @(s) 27/48*(6*s + 2);
    end
    
    jac = 1;  % jacobian
    
    A0 = zeros(length(phi0));
    A1 = zeros(length(phi1));
    A2 = zeros(length(phi2));
    
    % create local mass (u) and bending (u'') matrices
    for i = 1:length(phi0)
        for j = i:length(phi0)
            
            f0 = @(s) phi0{i}(s)*phi0{j}(s);  % function to integrate
            % integration of f0 over -1,1 using 4 point quadrature
            A0(i,j) = wi(1)*f0(xi(1)) + wi(2)*f0(xi(2)) + wi(3)*f0(xi(3)) + wi(4)*f0(xi(4));
            
            f1 = @(s) phi1{i}(s)*phi1{j}(s);  % derivative function to integrate
            % integration of f1 over -1,1 using 4 point quadrature
            A1(i,j) =  wi(1)*f1(xi(1)) + wi(2)*f1(xi(2)) + wi(3)*f1(xi(3)) + wi(4)*f1(xi(4));
            
            f2 = @(s) phi2{i}(s)*phi2{j}(s);  % second derivative function to integrate
            % integration of f2 over -1,1 using 4 point quadrature
            A2(i,j) =  wi(1)*f2(xi(1)) + wi(2)*f2(xi(2)) + wi(3)*f2(xi(3)) + wi(4)*f2(xi(4));
            
        end
    end
    
    % copy upper triangular part to lower triangular part for symmetry
    A0 = A0 + triu(A0,1)';
    A1 = A1 + triu(A1,1)';
    
    allA{p} = A0;
    allAprm{p} = A1;
    
    
    % init global matrices
    glbA0 = zeros(N*P+1);
    glbA1 = zeros(N*P+1);
    glbA2 = zeros(N*P+1);
    
    t0 = 0; tf = 1;
    k = 1; beta = 0;
    a = 0; b = 0;  % b.c. values
    
    % constants
    q = 200;    % lb/in2
    D = 8.8e7;  % lb*in
    L = 50;     % in
    
    
    
    sList = [100, 1000, 10000];  % S values,  [lb/in]
    
    
    L2 = zeros(length(nn),length(sList));
    Linf = zeros(length(nn),length(sList));
    
    nrms = table(nn);
    nrms.Properties.VariableNames = {'n'};
    
    for i = 1:length(nn);
        
        n = nn(i);
        
        h = (tf-t0)/(n+1);  % set step size
        t = [t0+h:h:tf-h]'; %#ok<NBRAK>
        soln = table(t);  % init soln data table with FD points
        
        for j = 1:length(sList)
            
            S = sList(j);
            
            gamma = S*(L^2)/D;
            
            Z = q*(L^4)/(2*D);
            
            % build f(t) vector
            ft = Z*t.*(1-t);
            
            f = ft*(h^2);  % init f vector
            
            %         % create coefficient vector
            %         c(1) = -k;
            %         c(2) = 2*k + gamma*h^2 ;
            %         c(3) = -k;
            %
            %         % create A matrix
            %         A = diag(c(1)*ones(n-1,1),-1) + ...
            %             diag(c(2)*ones(n,1)) + ...
            %             diag(c(3)*ones(n-1,1),+1);
            
            % boundary conditions
            f(1) = f(1) - c(1)*a;      % left
            f(end) = f(end) - c(3)*b;  % right
            
            % solve AU=f matrix equation
            uapx = A0\f;  % u approximate
            
            % exact solution
            u_exact = @(t) Z/gamma*(-t.^2 + t-2/gamma + ...
                2/(gamma * sinh(sqrt(gamma))) * ...
                (sinh(sqrt(gamma)*t) + sinh(sqrt(gamma)*(1-t))));
            uext = u_exact(t);
            
            % compute error
            e = uext - uapx;
            
            % save norms data to matrices
            L2(i,j) = norm(e,2);
            Linf(i,j) = norm(e,Inf);
            
            % take abs value and add eps to error for graphing purposes
            e = abs(e) + eps;
            
            % save solution data to table
            solnsfix = ['S',num2str(S)];
            solnvars = {['uext', solnsfix],['uapx', solnsfix],['e', solnsfix]};
            tmp = table(uext,uapx,e);
            tmp.Properties.VariableNames = solnvars;
            soln = [soln, tmp]; %#ok<AGROW> % concatenate tmp table with solutions table
            
        end
        
        % write solutions table to data file for each n value
        writetable(soln,[prbsfx,'n',num2str(n),'soln','.dat'],'Delimiter','\t');
        
    end
    
    % export data to .dat file
    
    nrms = zeros(size(L2,1),size(L2,2)*2+1);
    nrmsvars = {'n'};
    nrms(:,1) = nn;
    
    for j = 1:size(L2,2)
        
        nrms(:,j*2) = L2(:,j);
        nrms(:,j*2+1) = Linf(:,j);
        
        S = sList(j);
        nrmssfix = ['S',num2str(S)];
        nrmsvars = [nrmsvars, {['L2', nrmssfix], ['Linf', nrmssfix]}]; %#ok<AGROW>
        
    end
    
    nrmtbl = array2table(nrms);
    nrmtbl.Properties.VariableNames = nrmsvars;
    writetable(nrmtbl,[prbsfx,'nrms','.dat'],'Delimiter','\t');
    
    
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



