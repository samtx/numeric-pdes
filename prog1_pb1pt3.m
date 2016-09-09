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
prob = 1;
part = 3;
prbsfx = [data_dir,'hw',num2str(hw),'pb',num2str(prob),'pt',num2str(part)];

x0 = 0; xf = 1;
k = 1; fx = 1; a = 0; b = 0;
betalist = 1;  % beta values
% for beta = betalist
%
% end

L2 = zeros(length(nn),length(betalist));
Linf = zeros(length(nn),length(betalist));

nrms = table(nn);
nrms.Properties.VariableNames = {'n'};

for i = 1:length(nn);
    
    n = nn(i);
    
    h = (xf-x0)/(n+1);  % set step size
    x = [x0+h:h:xf-h]'; %#ok<NBRAK>
    soln = table(x);  % init soln data table with FD points
    
    for j = 1:length(betalist)
        
        beta = betalist(j);
        
        f = ones(n,1)*fx*(h^2);  % init f vector
        
        % create coefficient vector
        c(1) = -k - beta*h/2;
        c(2) = 2*k ;
        c(3) = -k + beta*h/2;       
        
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
        z = exp(beta) - 1;
        u_exact = @(X) 1/(z*beta)*(1 + z*X - exp(beta*X));
        uext = u_exact(x);
        
        % compute error
        e = uext - uapx;
        
        % save norms data to matrices
        L2(i,j) = norm(e,2);
        Linf(i,j) = norm(e,Inf);
        
        % take abs value and add eps to error for graphing purposes
        e = abs(e) + eps;      
        
        % save solution data to table
        solnsfix = ['be',num2str(beta)];
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
    
    beta = betalist(j);
    nrmssfix = ['be',num2str(beta)];
    nrmsvars = [nrmsvars, {['L2', nrmssfix], ['Linf', nrmssfix]}]; %#ok<AGROW>
    
end

nrmtbl = array2table(nrms);
nrmtbl.Properties.VariableNames = nrmsvars;
writetable(nrmtbl,[prbsfx,'nrms','.dat'],'Delimiter','\t');


%         plot(x,uext,'-k',x,uapx,'--r','Linewidth',2)
%         legend('Exact Solution','Approx Solution');
%         xlabel('x');ylabel('y')








