function [plot1, plot2] = plot_soln(elems,X,Y,uapx,uext)

% graph FE approx and exact solutions
tr = triangulation(elems,X,Y,uapx);   % Approx Solution
subplot(1,2,1)
plot1 = trisurf(tr);
axis([0 1 0 1 -1 1]);
title('Approx'); xlabel('X'); ylabel('Y'); zlabel('u_h');
tr = triangulation(elems,X,Y,uext);  % Exact solution
subplot(1,2,2)
plot2 = trisurf(tr);
title('Exact'); xlabel('X'); ylabel('Y'); zlabel('U');

end