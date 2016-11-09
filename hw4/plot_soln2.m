function [plot1] = plot_soln2(elems,X,Y,uapx)

% graph FE approx solution
tr = triangulation(elems,X,Y,uapx);   % Approx Solution
plot1 = trisurf(tr);
axis([0 1 0 1 -1 1]);
title('Approx'); xlabel('X'); ylabel('Y'); zlabel('u_h');

end