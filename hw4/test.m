

% p = 3;
% n = 10;
% c = 1;
% origfname = 'p3c.1';


fname = 'mesh/p2n40';

[n, verts, elems, bounds] =  get_mesh(fname);

draw_mesh(verts, elems);