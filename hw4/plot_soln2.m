function [plot1] = plot_soln2(elems,X,Y,uapx)

% graph FE approx solution
tr = triangulation(elems,X,Y,uapx);   % Approx Solution
plot1 = trisurf(tr);
view(0,90)
% axis([0 1 0 2]);
title('Approx'); xlabel('X'); ylabel('Y'); zlabel('u_h');
colorbar

end