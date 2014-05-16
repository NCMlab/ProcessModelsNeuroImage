function [vector,value]=eig_habeck(matrix)
 
% function [vector,value]=eig_habeck(matrix)
%
% Modification of Matlab eigen vector decomposition,
% descending order according to size of eigen value, and don't
% output to the screen!
 
 
dim=size(matrix);
 
if(dim(1)~=dim(2))
  display([' Error: matrix is not quadratic!'])
end
[v,lambda]=eig(matrix);
 
llambda=diag(lambda);
 
value=llambda(dim(1):-1:1);
vector=v(:,dim(1):-1:1);
