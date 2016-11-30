

% p = 3;
% n = 10;
% c = 1;
% origfname = 'p3c.1';
clf
nnode = 3;  % nodes per element
q = 5;  % constant

prob = 1;
N = 40;
fname = ['mesh/p',num2str(prob),'n',num2str(N)];

[n, verts, elems, bounds] =  get_mesh(fname);
% draw_mesh(verts,elems);

if prob == 1
       
    [uapx, uext, time, L2, H1, X, Y] = hw4pb1(n, verts, elems);
    
    [plot1, plot2] = plot_soln(elems, X, Y, uapx, uext);
    figure;
    plot_soln2(elems, X, Y, abs(uext-uapx));
    title('Error');
    
   % Print norms to screen
   fprintf(' Norms for N = %i\n',N);
   fprintf('    L2 Norm          H1 Norm   \n');
   fprintf('--------------  ---------------\n');
   fprintf('%12.8e  %12.8e\n',[L2, H1]);
    
elseif prob == 2
    
    f = @(x,y) 1;
    g = 1;
    
    [uapx, time, X, Y] = fem2d_elliptic(n, verts, elems, bounds, f, g);
    
    [plot1] = plot_soln2(elems, X, Y, uapx);
    
elseif prob == 3
    
    f = @(x,y) x.*y;
    g = 1;
    
    [uapx, time, X, Y] = fem2d_elliptic(n, verts, elems, bounds, f, g);
    
    [plot1] = plot_soln2(elems, X, Y, uapx);
    
    
end





