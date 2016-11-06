% Create element stiffness matrix
function [I, J, K] = get_elem_stiffness_matrix(e,v)


% basis functions lambda... barycentric coordinates

x = v(:,1); y = v(:,2);  % get x,y coordinates of element vertices
area = abs(.5*det([x';y';1,1,1]));  % area of element triangle
nab = area*2;  % triangle symbol... nabla
cyc = [1,2,3,1,2,3];
a = zeros(3,1); b = zeros(3,1); c = zeros(3,1);
for i = [1,2,3]
    j = cyc(i+1);
    k = cyc(i+2);
    a(i) = x(j)*y(k)-x(k)*y(j);
    b(i) = y(j) - y(k);
    c(i) = x(k) - x(j);
%     lmb(i) = (a + b + c)/area;
end

% create element stiffness matrix from b,c
eA1 = [b(1)*b(1)+c(1)*c(1), b(1)*b(2)+c(1)*c(2), b(1)*b(3)+c(1)*c(3);
                         0, b(2)*b(2)+c(2)*c(2), b(2)*b(3)+c(2)*c(3);
                         0,                   0, b(3)*b(3)+c(3)*c(3)];
eA1 = eA1 + triu(eA1,1)'; % copy upper triangular part to lower triangular part for symmetry
eA1 = 1/(2*nab)*eA1;  % multiply by area of triangle


% save nodal values in index vectors for sparse global matrix
nrow = size(eA1,1); ncol = size(eA1,2);
nn = nrow * ncol;
I = zeros(nn,1); J = zeros(nn,1); K = zeros(nn,1);
idx = 0;  % idx of triplet
for i = 1:nrow
    for j = 1:ncol
        idx = idx + 1;
        I(idx) = e(i);
        J(idx) = e(j);
        K(idx) = eA1(i,j);
    end
end

end

% for k = 1:size(elems,2)
% e = elems(k,:);  % element vertex indicies
% v = verts(e,:);
% end
 