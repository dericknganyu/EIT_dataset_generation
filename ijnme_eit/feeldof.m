function [index]=feeldof(nd,nnel,ndof)
% Purpose: compute system dofs with each element
% Variables:
% index: system dof vector with element nd
% nd: element # whose system dofs
% nnel: # nodes per element
% ndof: # dofs per node
edof=nnel*ndof;
k=0;
for i=1:nnel
    start=(nd(i)-1)*ndof;
    for j=1:ndof
        k=k+1;
        index(k)=start+j;
    end
end
%--------------------------------------------------------------------------
