function MM = massmat(fem,acoef)
%-------------------------------------------------------------------------%
% forming the global stiffness matrix                                     %
% JIN Bangti (kimbtsing@yahoo.com.cn) Feb. 2, 2009                        %
%-------------------------------------------------------------------------%

% extrac fem parameters
nnel   = fem.nnel;
ndof   = fem.ndof;
nel    = fem.nel;
gcoord = fem.gcoord;
nodes  = fem.nodes;

% computation of the mass matrix
MM = sparse(fem.sdof,fem.sdof);
if nargin<2
    acoef = ones(fem.sdof,1);
end
for iel=1:nel
    nd(1)=nodes(iel,1);nd(2)=nodes(iel,2);nd(3)=nodes(iel,3);
    index=feeldof(nd,nnel,ndof);
    x1=gcoord(nd(1),1); y1=gcoord(nd(1),2);
    x2=gcoord(nd(2),1); y2=gcoord(nd(2),2);
    x3=gcoord(nd(3),1); y3=gcoord(nd(3),2);
    m = femasstriang(x1,y1,x2,y2,x3,y3,acoef(index));
    MM(index,index) = MM(index,index) + m;
end