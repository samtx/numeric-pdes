% Numeric PDEs
% Programming Assignment #6
% Sam Friedman
% 10/11/2016

% ref: http://www.math.chalmers.se/~mohammad/teaching/PDEbok/Draft_I+II.pdf

P = 2;  %  interpolating polynomial degree
NN = [5, 10, 20]; % number of finite elements
xnum = 20;  % number of additional points per finite element

% constants
S = 100;    % lb/in
q = 200;    % lb/in2
D = 8.8e7;  % lb*in
L = 50;     % in
Z = q*(L^4)/(2*D);

t0 = 0; tf = 1;   % domain
u0 = 0; uf = 0;   % boundary conditions

allU = {};
allE = {}; 
allT = {};
allUe = {};
allTprime = {};

% approx u(x) ~ U(x)
% Prob 1
data_dir = 'data/';
hw = 6;
prob = 1;
part = 1;
prbsfx = [data_dir,'hw',num2str(hw),'pb',num2str(prob),'pt',num2str(part)];

gauss_pts = 4;  % use 4 quadrature points for all polynomial basis functions
[wi, xi] = gauss_pts_1D(gauss_pts);

[phi0, phi1] = basis_functions_1D(P);  % get basis functions for mass and stiffness matrices

for k = 1:length(NN)
    N = NN(k);  % number of finite elements
    
    h = (tf-t0)/N;  % set step size
    t = [t0:h/P:tf]'; %#ok<NBRAK> % list of t points at all nodes %#ok<NBRAK>
    soln = table(t(2:end-1));  % init soln data table with FD points except boundary points
    
    % ------------------- Create Local matrices ------------------
    
    % init global matrices
    glbA0 = zeros(N*P+1);
    glbA1 = zeros(N*P+1);
    
    % create local mass (u) and stiffness (u') matrices
    [A0, A1] = local_matrices_1D(phi0, phi1, wi, xi);
    
    % scale local matrices for step size
    A0 = A0 * h/2;
    A1 = A1 * 2/h;
    
    % ------------------- Create Global matrices ------------------
    FEidx = zeros(N,P+1);  % index each element to their respective nodes
    % loop over the finite elements and create global matrix.
    for n = 1:N
        % add local matrix to global matrix
        idx1 = n*P-(P-1);% begin row/column = n*P-1
        idx2 = (n+1)*P-(P-1); % end row/column = (n+1)*P-1
        glbA0(idx1:idx2,idx1:idx2) = glbA0(idx1:idx2,idx1:idx2) + A0;
        glbA1(idx1:idx2,idx1:idx2) = glbA1(idx1:idx2,idx1:idx2) + A1;
        
        FEidx(n,:) = idx1:idx2;  % which nodes are in each element
        % each row # is the element ID
        % each column in the row is the node ID
    end
    
    % ----------    Create force vector -----------------
    % build f(t) vector
    f = Z*t.*(1-t);
    f = glbA0 * f;
    
    % ----------------- apply boundary conditions -----------------
    tmpA0 = glbA0; tmpA1 = glbA1;
    glbA0(1,:)   = [];   glbA1(1,:) = [];   f(1) = [];   % remove first row
    glbA0(end,:) = []; glbA1(end,:) = []; f(end) = [];  % remove last row
    
    % apply constant to mass matrix
    gamma = S*(L^2)/D;
    glbA0  = glbA0 * gamma;
    
    % add mass and stiffness matrices
    glbA = glbA0 + glbA1;
    
    % -------------- solve AU = f matrix equation --------------
    uapx = glbA\f;  % u approximate
    uapx = uapx(2:end-1);  % remove boundary entries
    
    % exact solution
    u_exact = @(t) Z/gamma*(-t.^2 + t - 2/gamma + ...
        2/(gamma * sinh(sqrt(gamma))) * ...
        (sinh(sqrt(gamma)*t) + sinh(sqrt(gamma)*(1-t))));
    uext = u_exact(t(2:end-1));
    
    u_prime = @(t) Z/gamma*( ...
        ( 2*csch(sqrt(gamma))*( sqrt(gamma)*cosh(t*sqrt(gamma))-sqrt(gamma)*cosh(sqrt(gamma)*(1-t))))/gamma - 2*t + 1);
    
    
    % ----------------- compute error -------------------
    e = uext - uapx;
    
    % save norms data to matrices
    L2 = norm(e,2);
    Linf = norm(e,Inf);
    H1 = sqrt(e'*glbA1(:,2:end-1)*e + e'*glbA0(:,2:end-1)*e);
    
    %             glbA0 = tmpA0; glbA1 = tmpA1;  % replace temporary global matrices
    
    % take abs value and add eps to error for graphing purposes
    e = abs(e) + eps;
    
    % save solution data to table
    solnvars = {'uext','uapx','e'};
    tmp = table(uext,uapx,e);
    tmp.Properties.VariableNames = solnvars;
    soln = [soln, tmp]; % concatenate tmp table with solutions table
    
    % write solutions table to data file for each n value
    writetable(soln,[prbsfx,'ply',num2str(P),'n',num2str(n),'soln','.dat'],'Delimiter','\t');
    
    % Find errors within each element
    %     snum = P + 1;  % number of s pts per element
    
    % s in [-1, 1]
    
    % add the boundary conditions back into uapx, uext, and e vectors
    t;   % values of t that the nodes are computed at
    uapx = [0;uapx;0];
    uext = [0;uext;0];
    e = [0;e;0];
    
    % xnum = 5;  % 5 points between nodes
    U = [];
    T = [];
    Uprime = [];
    Tprime = [];
    Ugauss = [];
    Tgauss = [];
    % U = zeros(length(uapx)+xnum*(length(uapx)-1),1);  % values at all points, including previously computed nodes
    % T = zeros(length(uapx)+xnum*(length(uapx)-1),1);  % t value at all computed points
    for i = 1:length(uapx)-1
        u1 = uapx(i);   % first node value
        u2 = uapx(i+1); % second node value
        t1 = t(i);   % domain at first node
        t2 = t(i+1);  % domain at second node
        dx = (t2-t1)/(xnum+1);  % distance between extra points between the two nodes
        xx = t1+dx:dx:t2-dx;   % list of new domain points to interpolate between
        idx = (i-1)*xnum + 1;  % starting index of extra points
        y = zeros(xnum,1);
%         yprime = zeros(xnum,1);
        for j = 1:xnum
            x = xx(j); % domain at extra point that is being evaluated
            s = (2*(x-t1)/(t2-t1));  % transformed point into s space of [-1, 1]
            y(j) = u1*(1-s)/2 + u2*(1+s)/2;  % interpolated value of that point
%             yprime(j) = u1*(-s)/2 + u2*s/2;
        end
        % add y values and original node values into master list
        U = [U;u1;y;u2];
        T = [T;t1;xx';t2];
%         Uprime = [Uprime;yprime];
%         Tprime = [Tprime;xx'];
    end
    
    Ue = u_exact(T);
    
    E = abs(Ue-U) + eps;

    allT{k} = T;
    allU{k} = U;
    allUe{k} = Ue;
    allE{k} = E;
    

    
    
    % Take central finite difference for U and Ugauss to get derivative
    Uprime = zeros(length(T)-1,1);
    Tprime = T(2:end-1);
    for i = 2:length(Tprime)
        Uprime(i) = (U(i+1)-U(i-1))/(2*(T(i+1)-T(i)));
    end
    %remove first entry
    Uprime(1) = [];
    Ueprime = u_prime(Tprime);
    Eprime = abs(Ueprime-Uprime) + eps;
    
    allTprime{k} = Tprime;
    allUprime{k} = Uprime;
    allUeprime{k} = Ueprime;
    allEprime{k} = Eprime;
    allUgauss{k} = Ugauss;
%     allTgauss{k} = Tgauss;
%     allEgauss{k} = Egauss;
    
    % print norms to screen
    fprintf('Norms for P=%d,  N=%d   \n',[P,N]);
    fprintf(' max(e) = %.8e  \n',norm(E,Inf));
    fprintf('  L2(e) = %.8e  \n',norm(E,2));
    fprintf(' L2(ep) = %.8e  \n',norm(Eprime,2));
    
end

% subplot(1,2,1);
% plot(allT{3},allUe{3},'-',allT{1},allU{1},'--',allT{2},allU{2},'--',allT{3},allU{3},'--','Linewidth',1);
% legend('U Exact','N=5','N=10','N=20');
% subplot(1,2,2);
plot(allT{1},allE{1},allT{2},allE{2},allT{3},allE{3},'Linewidth',1);
legend('N=5','N=10','N=20');
title(['Error with ',num2str(xnum),' extra points']);
% return
% subplot(1,2,1);
% plot(allTprime{3},allUeprime{3},'-',allTprime{1},allUprime{1},'--',allTprime{2},allUprime{2},'--',allTprime{3},allUprime{3},'--','Linewidth',1);
% legend('U Exact','N=5','N=10','N=20');
% subplot(1,2,2);
% plot(allTprime{1},allEprime{1},allTprime{2},allEprime{2},allTprime{3},allEprime{3},'Linewidth',1);
% legend('N=5','N=10','N=20');
% title(['Derivative Error with ',num2str(xnum),' extra points']);

% subplot(1,2,1);
% plot(allTprime{3},allUeprime{3},'-',allTprime{1},allUprime{1},'--',allTprime{2},allUprime{2},'--',allTprime{3},allUprime{3},'--','Linewidth',1);
% legend('U Exact','N=5','N=10','N=20');
% subplot(1,2,2);
% plot(allTprime{1},allEprime{1},allTprime{2},allEprime{2},allTprime{3},allEprime{3},'Linewidth',1);
% legend('N=5','N=10','N=20');
% title(['Derivative Error with ',num2str(xnum),' extra points']);

return
% should be an even number
% P = 2;
% if P == 2
%     ds = 1;
%     spts = [-1, 0, 1];
%     dx = 2/(xnum+1);
%     %         dx = 1/6;  % xnum = 10
%     %         xpts = dx * [-xnum/2:1:xnum/2];
%     %         xpts = dx*[-5,-4,-3,-2,-1,1,2,3,4,5];
% elseif P == 3
%     ds = 2/3;
%     spts = -1:ds:1;  % t values at u nodes
%     dx = 2/(xnum+1);
%     %         xpts = dx*[-5,-4,-3,-2,-1,1,2,3,4,5];
%
% else
%     ds = 2/P;
%     spts = -1:ds:1;
%     dx = 2/(xnum+1);
%
% end
% % xpts = dx * [-xnum/2:-1, 1:xnum/2]
% xpts = -1+dx:dx:1-dx;
%
%

% -1+dx:dx:1-dx;
%
% allx = zeros(xnum,1);  % number of total additional points in each FE
% allxval = zeros(xnum*N,1);  % total num of additional points in domain
% %loop over each finite element
% for n = 1:N
%     fe = FEidx(n,:);  % get the nodes in the particular finite element
%     %loop over each extra x point
%     for xidx = 1:length(xval)
%         x = xpts(xidx);
%         %loop over each basis function
%         for j = 1:length(phi0)  % number of polynomial basis functions
%             xpts(xidx) = xpts(xidx) + phi0{j}(x)*uapx(fe(j));
%         end
%     end
%     idx1 = (n-1)*xnum + 1;  % begin xpts in allxval
%     idx2 = n*xnum;          % end xpts in allxval
%     allxval(idx1:idx2) = allxval(idx1:idx2) + xpts;
% end
% xpts

%  ------------------------    Plot data   --------------------
xx = t(2:end-1);
subplot(1,2,1);
plot(t(2:end-1),uext,'-',xx,uapx,'--','Linewidth',1);
legend('U Exact','U Approx');
subplot(1,2,2);
plot(xx,e,'Linewidth',1);
legend('Error');

%
% xval = zeros(length(xpts),1);  % number of total additional points in each FE
% allxval = zeros(length(xpts)*N,1);  % total num of additional points in domain
% %loop over each finite element
% for n = 1:N
%     fe = FEidx(n,:);  % get the nodes in the particular finite element
%     %loop over each extra x point
%     for xidx = 1:length(xpts)
%         x = xpts(xidx);
%         %loop over each basis function
%         for j = 1:length(phi0)  % number of polynomial basis functions
%             xval(xidx) = xval(xidx) + phi0{j}(x)*uapx(fe(j));
%         end
%     end
%     idx1 = (n-1)*xnum + 1;  % begin xpts in allxval
%     idx2 = n*xnum;          % end xpts in allxval
%     allxval(idx1:idx2) = allxval(idx1:idx2)  + xval;
% end
% xval
% export data to .dat file

return


% -----------------  write norms data to table   -----------------
nrms = zeros(size(L2,1),size(L2,2)*3);
nrmsvars = {'N'};
for j = 1:size(L2,2)
    
    nrms(:,(j-1)*3+1) = L2(:,j);
    nrms(:,(j-1)*3+2) = Linf(:,j);
    nrms(:,(j-1)*3+3) = H1(:,j);
    
    S = sList(j);
    nrmssfix = ['S',num2str(S)];
    nrmsvars = [nrmsvars, {['L2', nrmssfix], ['Linf', nrmssfix], ['H1', nrmssfix]}];
    
end
nrms = [nn', nrms];
nrmtbl = array2table(nrms);
nrmtbl.Properties.VariableNames = nrmsvars;
writetable(nrmtbl,[prbsfx,'ply',num2str(P),'nrms','.dat'],'Delimiter','\t');


%  ------------------------    Plot data   --------------------
xx = t(2:end-1);
subplot(1,2,1);
plot(xx,uext,'-',xx,uapx,'--','Linewidth',1);
legend('U Exact','U Approx');
subplot(1,2,2);
plot(xx,e,'Linewidth',1);
legend('Error');




