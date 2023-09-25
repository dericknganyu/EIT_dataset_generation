function KK = stiffmat(fem,sigma)
%-------------------------------------------------------------------------%
% forming the global stiffness matrix                                     %
% JIN Bangti (kimbtsing@yahoo.com.cn) Sept. 21, 2010                      %
%-------------------------------------------------------------------------%

% extracting finite element parameters
nnel   = fem.nnel;
ndof   = fem.ndof;
nel    = fem.nel;
gcoord = fem.gcoord;
nodes  = fem.nodes;

KK=sparse(fem.sdof,fem.sdof);
if ischar(sigma)
    x = gcoord(:,1);
    y = gcoord(:,2);
    sigma = eval(sigma);
end
for iel=1:nel
    nd(1)=nodes(iel,1);nd(2)=nodes(iel,2);nd(3)=nodes(iel,3);
    index=feeldof(nd,nnel,ndof);
    x1=gcoord(nd(1),1); y1=gcoord(nd(1),2);
    x2=gcoord(nd(2),1); y2=gcoord(nd(2),2);
    x3=gcoord(nd(3),1); y3=gcoord(nd(3),2);
    [k]=festifftriang(x1,y1,x2,y2,x3,y3,sigma(index));
    KK(index,index) = KK(index,index) + k;
end