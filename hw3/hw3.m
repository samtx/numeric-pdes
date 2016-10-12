% Numeric PDEs
% Programming Assignment #3
% Sam Friedman
% 10/11/2016


% ref: http://www.math.chalmers.se/~mohammad/teaching/PDEbok/Draft_I+II.pdf

poly = [2:3]';  %#ok<NBRAK> % polynomial degrees
nn = [10,20,40]; % number of elements
sList = [100, 1000, 10000];  % S values,  [lb/in]

% constants
q = 200;    % lb/in2
D = 8.8e7;  % lb*in
L = 50;     % in
Z = q*(L^3)/(2*D);

% domain
t0 = 0; tf = 1;

% boundary conditions
u0 = 0; uf = 0;

% approx u(x) ~ U(x)
% Prob 1
data_dir = 'data/';
hw = 3;
prob = 1;
part = 1;
prbsfx = [data_dir,'hw',num2str(hw),'pb',num2str(prob),'pt',num2str(part)];


% Gauss Quadrature weights and abscissa
% ref: https://pomax.github.io/bezierinfo/legendre-gauss.html

gauss_pts = 4;  % use 4 quadrature points for all polynomial basis functions

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

for pp = 1:length(poly);
    
    P = poly(pp);  % degree polynomial
    
    phi0 = {};
    phi1 = {};
    
    
    
    % get shape basis functions for element, Lagrange interpolation (-1,1)
    switch P
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
    
    L2 = zeros(length(nn),length(sList));
    Linf = zeros(length(nn),length(sList));
    H1 = zeros(length(nn),length(sList));
    
    for nidx = 1:length(nn);
        
        N = nn(nidx);  % number of finite elements
        %         jac = 1;  % jacobian
        
        h = (tf-t0)/N;  % set step size
        t = [t0:h/P:tf]'; %#ok<NBRAK> % list of t points at all nodes %#ok<NBRAK>
        soln = table(t(2:end-1));  % init soln data table with FD points except boundary points
        
        % ------------------- Create global matrices ------------------
        
        % init global matrices
        glbA0 = zeros(N*P+1);
        glbA1 = zeros(N*P+1);
        
        % loop over the finite elements and create local mass and stiffness matrices.
        for n = 1:N
            
            % init local matrices
            A0 = zeros(length(phi0));
            A1 = zeros(length(phi1));
            
            % create local mass (u) and stiffness (u') matrices
            for i = 1:length(phi0) %#ok<FXSET>
                for j = i:length(phi0)
                    
                    f0 = @(s) phi0{i}(s)*phi0{j}(s);  % function to integrate
                    % integration of f0 over -1,1 using 4 point quadrature
                    A0(i,j) = wi(1)*f0(xi(1)) + wi(2)*f0(xi(2)) + wi(3)*f0(xi(3)) + wi(4)*f0(xi(4));
                    
                    f1 = @(s) phi1{i}(s)*phi1{j}(s);  % derivative function to integrate
                    % integration of f1 over -1,1 using 4 point quadrature
                    A1(i,j) =  wi(1)*f1(xi(1)) + wi(2)*f1(xi(2)) + wi(3)*f1(xi(3)) + wi(4)*f1(xi(4));
                    
                end
            end
            
            % copy upper triangular part to lower triangular part for symmetry
            A0 = A0 + triu(A0,1)';
            A1 = A1 + triu(A1,1)';
            
            % scale local matrices for step size
            A0 = A0 * h/2;
            A1 = A1 * 2/h;
            
            % add local matrix to global matrix
            idx1 = n*P-(P-1);% begin row/column = n*P-1
            idx2 = (n+1)*P-(P-1); % end row/column = (n+1)*P-1
            glbA0(idx1:idx2,idx1:idx2) = glbA0(idx1:idx2,idx1:idx2) + A0;
            glbA1(idx1:idx2,idx1:idx2) = glbA1(idx1:idx2,idx1:idx2) + A1;
            
        end
        
        % ---------------- Choose S value --------------------
        for sidx = 1:length(sList)
            
            % ----------    Create force vector -----------------
            % build f(t) vector
            f = Z*t.*(1-t);
            
            S = sList(sidx);
            gamma = S*(L^2)/D;
            
            % apply constant to mass matrix
            glbA0  = glbA0 * gamma;
            
            % add mass and stiffness matrices
            glbA = glbA0 + glbA1;
            
            % ----------------- apply boundary conditions -----------------
%             f = f - glbA(:,1)*u0;
%             f = f - glbA(:,end)*uf;
            glbA(1,:) = [];     % remove first row
            glbA(end,:) = [];   % remove last row
%             glbA(:,1) = [];     % remove first column
%             glbA(:,end) = [];   % remove last column
            f(1) = [];   % remove first entry in f vector
            f(end) = []; % remove last entry in f vector
            
            % -------------- solve AU = f matrix equation --------------
            uapx = glbA\f;  % u approximate
            
            % exact solution
            u_exact = @(t) Z/gamma*(-t.^2 + t-2/gamma + ...
                2/(gamma * sinh(sqrt(gamma))) * ...
                (sinh(sqrt(gamma)*t) + sinh(sqrt(gamma)*(1-t))));
            uext = u_exact(t(2:end-1));
            
            % ----------------- compute error -------------------
            e = uext - uapx;
            
            % save norms data to matrices
            L2(nidx,sidx) = norm(e,2);
            Linf(nidx,sidx) = norm(e,Inf);
            H1(nidx,sidx) = .0001;
            
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
        writetable(soln,[prbsfx,'ply',num2str(P),'n',num2str(n),'soln','.dat'],'Delimiter','\t');
        
    end
    
    % export data to .dat file
    
    nrms = zeros(size(L2,1),size(L2,2)*3);
    nrmsvars = {'N'};
    for j = 1:size(L2,2)
        
        nrms(:,(j-1)*3+1) = L2(:,j);
        nrms(:,(j-1)*3+2) = Linf(:,j);
        nrms(:,(j-1)*3+3) = H1(:,j);
        
        S = sList(j);
        nrmssfix = ['S',num2str(S)];
        nrmsvars = [nrmsvars, {['L2', nrmssfix], ['Linf', nrmssfix], ['H1', nrmssfix]}]; %#ok<AGROW>
        
    end
    
    nrms = [nn', nrms]; %#ok<AGROW>
    nrmtbl = array2table(nrms);
    nrmtbl.Properties.VariableNames = nrmsvars;
    writetable(nrmtbl,[prbsfx,'ply',num2str(P),'nrms','.dat'],'Delimiter','\t');
    
end

%    soln = table(x);  % init soln data table with FD points

%     % save solution data to table
%     solnsfix = ['ply',num2str(ply)];
%     solnvars = {['uext', solnsfix],['uapx', solnsfix],['e', solnsfix]};
%     tmp = table(uext,uapx,e);
%     tmp.Properties.VariableNames = solnvars;
%     soln = [soln, tmp]; %#ok<AGROW> % concatenate tmp table with solutions table
%
% end
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



