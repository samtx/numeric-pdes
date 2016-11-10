

% p = 3;
% n = 10;
% c = 1;
% origfname = 'p3c.1';

nnode = 3;  % nodes per element
q = 5;  % constant

prob = 3;
g = 1;
fname = 'mesh/p3n40';

[n, verts, elems, bounds] =  get_mesh(fname);
% draw_mesh(verts,elems);

if prob == 1
    
    [uapx, uext, time, L2, H1, X, Y] = hw4pb1(n, verts, elems);
    
    [plot1, plot2] = plot_soln(elems, X, Y, uapx, uext);
    
elseif prob == 2
    
    [uapx, time, X, Y] = hw4pb2(n, verts, elems, bounds, g);
    
    [plot1] = plot_soln2(elems, X, Y, uapx);
    
elseif prob == 3
    
    f = @(x,y) x.*y;
        
    [uapx, time, X, Y] = fem2d_elliptic(n, verts, elems, bounds, f, g);
    
    [plot1] = plot_soln2(elems, X, Y, uapx);
    
    
end





