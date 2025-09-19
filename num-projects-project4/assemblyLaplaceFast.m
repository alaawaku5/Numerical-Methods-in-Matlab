function A = assemblyLaplaceFast(coordinates,elements)

% We aim to build the P1-FEM matrix for the 2D Poisson problem
% in (almost) linear complexity.

nC = size(coordinates,1);

% Instead of updating the matrix, in the element loop below,
% we will build vectors of these updates. Later we can then
% use the full power of the sparse command to build the P1-FEM
% matrix.

I = zeros(9*nC,1);
J = zeros(9*nC,1);
A = zeros(9*nC,1);

% Now comes the typical FEM loop over the elements

for i = 1:size(elements,1)
   nodes = elements(i,:);
   B = [1 1 1 ; coordinates(nodes,:)'];
   grad = B \ [0 0 ; 1 0 ; 0 1];
 
   % From here on you should modify the code to compute the vectors
   % of the matrix updates.
   %
   % old code: A(nodes,nodes) = A(nodes,nodes) + det(B)*grad*grad'/2;
    
   Alocal = det(B)*(grad*grad')/2;
   I(1+9*(i-1):9*i) = ???
   J(1+9*(i-1):9*i) = ???
   A(1+9*(i-1):9*i) = ???
end

% Now, we can build the matrix by means of sparse:
            
A = sparse(I,J,A,nC,nC);
