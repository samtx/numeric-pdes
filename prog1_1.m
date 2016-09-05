% Numeric PDEs
% Programming Assignment #1
% Sam Friedman
% 9/4/2016

% PDE
% -k*u''(x) + beta*u'(x) = f(x)
% 0 < x < 1, u(0) = a, u(1) = b

% approx u(x) ~ U(x)
% Prob 1.1

prob = 1;
part = 1;
k = 1; f = 0; a = 0; b = 1;

for n = [25, 50, 100, 200]  % internal mesh points
    
    h = (b-a)/(n+1);  % set step size
    x = [a+h:h:b-h]'; %#ok<NBRAK>
    
    data = table(x);
    for beta = [1, 8, 64]  % beta values
        
        
        f = ones(n,1)*f;  % init f vector

        
        % create A matrix
        A = diag((-k - beta*h/2)*ones(n-1,1),-1) + ...
            diag(2*k*ones(n,1)) + ...
            diag((-k + beta*h/2)*ones(n-1,1),+1);
        
        % boundary conditions
        f(1) = f(1) - (-k - beta*h/2) * a;      % left
        f(end) = f(end) - (-k + beta*h/2) * b;  % right
        
        % solve AU=f matrix equation
        uapx = A\f;  % u approximate
        
        % exact solution
        u_exact = @(X) (1 - exp(beta*X))/(1 - exp(beta));
        uext = u_exact(x); 
        
        % compute error
        e = uext - uapx;
        
        plot(x,uext,'-k',x,uapx,'--r','Linewidth',2)
        legend('Exact Solution','Approx Solution');
        xlabel('x');ylabel('y')
        
        dlmwrite
        e_1p1b1n25 = e;
        uext_1p1b1n25 = u;
        uapx_1p1b1n25 = U;
        break
    end
    break
end

% export data to .dat file
