

function [uapx, time, X, Y] = fem2d_elliptic(n, verts, elems, bounds, f, g) 

nnode = 3;  % nodes per element
q = 1;  % constant


nV = n(1);  % number of vertices of all triangle elements
nE = n(2);  % number of 3-noded triangle elements
nB = n(3);  % number of boundary nodes
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
    eF = f(x,y);
    
    % add element triplets to global triplets vectors
    idx2 = idx1 + nnode^2 - 1; % end index
    gI0(idx1:idx2) = eI0;
    gJ0(idx1:idx2) = eJ0;
    gK0(idx1:idx2) = eK0;
    gI1(idx1:idx2) = eI1;
    gJ1(idx1:idx2) = eJ1;
    gK1(idx1:idx2) = eK1;
    idx1 = idx2 + 1;  % start index for next element
    gF(e) = eF;
    
end

% create global sparse matrices
gA0 = sparse(gI0, gJ0, gK0, nN, nN);
gA1 = sparse(gI1, gJ1, gK1, nN, nN);
gF = sparse(gF);

% apply boundary conditions

% loop over the boundary nodes, compute the mass matrix for 1D, 2 node,
% linear elements. multiply that by g, then add that into the appropriate
% nodes in the global 2D mass matrix
gIg = zeros(nB*2,1);
gJg = zeros(nB*2,1);
gKg = zeros(nB*2,1);

idx1 = 1;
G = [g, g];  % boundary condition

for i = 1:size(bounds,1)
    idx2 = idx1 + 1;
    b1 = bounds(i,1); b2 = bounds(i,2);  % boundary node index
    x1 = verts(b1,1); y1 = verts(b1,2);  % xy coordinates of bound node 1
    x2 = verts(b2,1); y2 = verts(b2,2);  % xy coordinates of bound node 2
    h = sqrt((x2-x1)^2+(y2-y1)^2);  % distance between nodes
    eA0 = h/6 * [2, 1;
                 1, 2];
    gIg(idx1:idx2) = [b1; b2];
    gJg(idx1:idx2) = [1; 1];
    gKg(idx1:idx2) =  eA0 * G'; % element g vector
    idx1 = idx2 + 1;
end
gG = sparse(gIg, gJg, gKg, nN, 1);
% Solve equation
A = gA1 + q*gA0;
b = gA0*gF + gG;
uh = A\b;
time = toc;

% get x-y coordinates of nodes
[vi, ~, uapx] = find(uh);
X = verts(vi,1); Y = verts(vi,2);  

            
end


