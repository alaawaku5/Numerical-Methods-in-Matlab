function A = assemblyLaplace(coordinates,elements)

% This function builds the P1-FEM matrix for the 2D Poisson problem.
% Even though the code appears to be of linear cost, one observes a
% quadratic growth of the runtime.

nC = size(coordinates,1);

A = sparse(nC,nC);
for i = 1:size(elements,1)
   nodes = elements(i,:);
   B = [1 1 1 ; coordinates(nodes,:)'];
   grad = B \ [0 0 ; 1 0 ; 0 1];
   A(nodes,nodes) = A(nodes,nodes) + det(B)*(grad*grad')/2;
end
