

% p = 3;
% n = 10;
% c = 1;
% origfname = 'p3c.1';

nnode = 3;  % nodes per element

fname = 'mesh/p1n10';

[n, verts, elems, bounds] =  get_mesh(fname);

nV = n(1);  % number of vertices of all triangle elements
nE = n(2);  % number of 3-noded triangle elements
nN = nV;    % total number of nodes
% draw_mesh(verts, elems);

e = elems(k,:);  % element vertex indicies
v = verts(e,:);  % element vertex positions in domain

% init global triplet indecies
gI0 = zeros(nE*(nnode^2),1);
gJ0 = zeros(nE*(nnode^2),1);
gK0 = zeros(nE*(nnode^2),1);
gI1 = zeros(nE*(nnode^2),1);
gJ1 = zeros(nE*(nnode^2),1);
gK1 = zeros(nE*(nnode^2),1);
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
    
    % add element triplets to global triplets vectors
    idx2 = idx1 + nnode^2 - 1; % end index
    gI0(idx1:idx2) = eI0;
    gJ0(idx1:idx2) = eJ0;
    gK0(idx1:idx2) = eK0;
    gI1(idx1:idx2) = eI1;
    gJ1(idx1:idx2) = eJ1;
    gK1(idx1:idx2) = eK1;
    idx1 = idx2 + 1;  % start index for next element    
    
end

% create global sparse matrices
gA0 = sparse(gI0, gJ0, gK0, nN, nN);
gA1 = sparse(gI1, gJ1, gK1, nN, nN);
toc;