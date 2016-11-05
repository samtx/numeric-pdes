function [n, verts, elems, bounds] = get_mesh(fname)

% p = 3;
% n = 10;
% c = 1;
% origfname = 'p3c.1';

% % fname = 'mesh/p2n10';
% 
% % fname = ['p',num2str(p),'n',num2str(n),'c',num2str(c)];
% % make them text files so that dlmread can work
% if 0
%     movefile([origfname,'.node'],[fname,'.node.txt']);
%     movefile([origfname,'.ele'],[fname,'.ele.txt']);
%     movefile([origfname,'.poly'],[fname,'.poly.txt']);
% end

% Get vertices
lines = dlmread([fname, '.node.txt']);
nvert = lines(1,1);  % number of vertices
verts = lines(2:end,2:3); % vertex positions
bounds = lines(2:end,4);  % boundary markers

% Get elements
lines = dlmread([fname, '.ele.txt']);
nelem = lines(1,1);     % number of elements
elems = lines(2:end,2:4);  % element markers

n = [nvert, nelem];

% Draw triangulation
% trimesh(elems,verts(:,1),verts(:,2))
end