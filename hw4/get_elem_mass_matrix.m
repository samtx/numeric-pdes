% Create element mass matrix
function [I, J, K] = get_elem_mass_matrix(e,v)

x = v(:,1); y = v(:,2);  % get x,y coordinates of element vertices
area = abs(.5*det([x';y';1,1,1]));  % area of element triangle

% create element mass matrix
eA0 = [2, 1, 1;
       1, 2, 1;
       1, 1, 2];
eA0 = area/12*eA0;  % multiply by area

% save nodal values in index vectors for sparse global matrix
nrow = size(eA0,1); ncol = size(eA0,2);
nn = nrow * ncol;
I = zeros(nn,1); J = zeros(nn,1); K = zeros(nn,1);
idx = 0;  % idx of triplet
for i = 1:nrow
    for j = 1:ncol
        idx = idx + 1;
        I(idx) = e(i);
        J(idx) = e(j);
        K(idx) = eA0(i,j);
    end
end

 