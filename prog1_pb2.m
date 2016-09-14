% Numeric PDEs
% Programming Assignment #1
% Sam Friedman
% 9/4/2016

% PDE
% -k*u''(x) + beta*u'(x) = f(x)
% 0 < x < 1, u(0) = a, u(1) = b

nn = [25, 50, 100, 200]';  % internal mesh points

% approx u(x) ~ U(x)
% Prob 1.1
data_dir = 'data/';
hw = 1;
prob = 2;
part = 1;
prbsfx = [data_dir,'hw',num2str(hw),'pb',num2str(prob),'pt',num2str(part)];

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
        
        % create coefficient vector
        c(1) = -k;
        c(2) = 2*k + gamma*h^2 ;
        c(3) = -k;       
        
        % create A matrix
        A = diag(c(1)*ones(n-1,1),-1) + ...
            diag(c(2)*ones(n,1)) + ...
            diag(c(3)*ones(n-1,1),+1);
        
        % boundary conditions
        f(1) = f(1) - c(1)*a;      % left
        f(end) = f(end) - c(3)*b;  % right
        
        % solve AU=f matrix equation
        uapx = A\f;  % u approximate
        
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


%         plot(x,uext,'-k',x,uapx,'--r','Linewidth',2)
%         legend('Exact Solution','Approx Solution');
%         xlabel('x');ylabel('y')








