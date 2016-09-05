% Numeric PDEs
% Programming Assignment #1
% Sam Friedman
% 9/4/2016

% PDE
% -k*u''(x) + beta*u'(x) = f(x)
% 0 < x < 1, u(0) = a, u(1) = b

n = [25, 50, 100, 200];  % internal mesh points

% approx u(x) ~ U(x)
% Prob 1.1
k = 1; f = 0; beta = [1, 8, 64]; a = 0; b = 1;

i = 1; j = 1;
h = (1-0)/(n(i)+1);
f = ones(n,1)*f;
u = zeros(n,1);

% create A matrix
A = diag((-k - beta(j)*h/2)*ones(n(i)-1,1),-1) + ...
    diag(2*k*ones(n(i),1)) + ...
    diag((-k + beta(j)*h/2)*ones(n(i)-1,1),+1);

% boundary conditions
f(1) = f(1) - (-k - beta(j)*h/2) * a;      % left
f(end) = f(end) - (-k + beta(j)*h/2) * b;  % right

U = A\f;  % solve AU=f matrix equation

% exact solution
u_exact = @(x) (1 - exp(beta(j)*x))/(1 - exp(beta(j))); 
u = u_exact(a+h:h:b-h);
u = u';

% compute error
e = u - U;
x = a+h:h:b-h;
plot(x,u,'-k',x,U,'--r','Linewidth',2)
legend('Exact Solution','Approx Solution');
xlabel('x');ylabel('y')