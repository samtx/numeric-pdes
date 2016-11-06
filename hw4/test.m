

% p = 3;
% n = 10;
% c = 1;
% origfname = 'p3c.1';

nnode = 3;  % nodes per element
q = 5;  % constant


fname = 'mesh/p1n40';

[n, verts, elems, bounds] =  get_mesh(fname);

nV = n(1);  % number of vertices of all triangle elements
nE = n(2);  % number of 3-noded triangle elements
nN = nV;    % total number of nodes
% draw_mesh(verts, elems);

% init global triplet indecies
gI0 = zeros(nE*(nnode^2),1);
gJ0 = zeros(nE*(nnode^2),1);
gK0 = zeros(nE*(nnode^2),1);
gI1 = zeros(nE*(nnode^2),1);
gJ1 = zeros(nE*(nnode^2),1);
gK1 = zeros(nE*(nnode^2),1);
gF  = zeros(nN,1);
idx1 = 1;  % start index for first element

% force vector
f = @(x,y) (5+10*pi^2)*cos(3*pi*x).*cos(pi*y);

tic;
for elemID = 1:size(elems,1)
    
    % choose element by ID
    e = elems(elemID,:);  % element vertex indicies
    v = verts(e,:);  % element vertex positions in domain
    
    % get element mass matrix triplets
    [eI0, eJ0, eK0] = get_elem_mass_matrix(e,v);
    
    % get element stiffness matrix triplets
    [eI1, eJ1, eK1] = get_elem_stiffness_matrix(e,v);
    
    % get element force vector
    %     eF = get_elem_force(e,v,f);
    x = v(:,1); y = v(:,2);  % get x,y coordinates of element vertices
    eFidx = e;
    %     eF = f(x,y);
    eF = (5+10*pi^2)*cos(3*pi*x).*cos(pi*y);
    
    % add element triplets to global triplets vectors
    idx2 = idx1 + nnode^2 - 1; % end index
    gI0(idx1:idx2) = eI0;
    gJ0(idx1:idx2) = eJ0;
    gK0(idx1:idx2) = eK0;
    gI1(idx1:idx2) = eI1;
    gJ1(idx1:idx2) = eJ1;
    gK1(idx1:idx2) = eK1;
    idx1 = idx2 + 1;  % start index for next element
    gF(e) = gF(e) + eF;
    
end

% create global sparse matrices
gA0 = sparse(gI0, gJ0, gK0, nN, nN);
gA1 = sparse(gI1, gJ1, gK1, nN, nN);
gF = sparse(gF);

% Solve equation
A = gA1 + q*gA0;
b = gA0*gF;
uh = A\b;
toc;

% graph FE approx solution
[vi, ~, Z] = find(uh);
X = verts(vi,1); Y = verts(vi,2);
% X = verts(:,1); Y = verts(:,2);
% Z = uh;
tr = triangulation(elems,X,Y,Z);
subplot(1,2,1)
trisurf(tr);
title('Approx'); xlabel('X'); ylabel('Y'); zlabel('u_h');

% Exact solution
uext = cos(3*pi*X).*cos(pi*Y);
tr = triangulation(elems,X,Y,uext);
subplot(1,2,2)
trisurf(tr);
title('Exact'); xlabel('X'); ylabel('Y'); zlabel('U');

