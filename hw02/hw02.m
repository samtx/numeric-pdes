% Numeric PDEs
% Programming Assignment #2
% Sam Friedman
% 9/27/2016


nn = 10;  % internal mesh points

poly = [1:6]';  %#ok<NBRAK> % polynomial degrees

% approx u(x) ~ U(x)
% Prob 1.1
data_dir = 'data/';
hw = 2;
prob = 1;
part = 1;
prbsfx = [data_dir,'hw',num2str(hw),'pb',num2str(prob),'pt',num2str(part)];

k = 1; fx = 0; a = 0; b = 1;
% for beta = betalist
%
% end

nrms = table(nn);
nrms.Properties.VariableNames = {'p'};

for i = 1:length(nn);
    
    n = nn(i);
    
    h = (b-a)/(n+1);  % set step size
    x = [a+h:h:b-h]'; %#ok<NBRAK>
    soln = table(x);  % init soln data table with FD points
    
    
    % Gauss Quadrature weights and abscissa
    % ref: https://pomax.github.io/bezierinfo/legendre-gauss.html
    switch n
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
            wi = [ 0.3607615730481386, 0.3607615730481386];
            xi = [];
        case 7
            wi = [];
            xi = [];
        case 8
            wi = [];
            xi = [];
    end
    
    
    for j = 1:length(poly)
        
        ply = poly(j);  % polynomial degree
        
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






