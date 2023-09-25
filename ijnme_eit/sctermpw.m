function ffs = sctermpw(fem,gradfcn)
%---------------------------------------------------------------------%
% generate the source term given piecewise sources                    %
% JIN Bangti (kimbtsing@yahoo.com.cn) Feb. 4, 2009                    %
%---------------------------------------------------------------------%

% extracting fe parameters
nnel   = fem.nnel;
ndof   = fem.ndof;
nel    = fem.nel;
gcoord = fem.gcoord;
nodes  = fem.nodes;

ffs=zeros(fem.sdof,1);
for iel=1:nel
    nd(1) = nodes(iel,1);
    nd(2) = nodes(iel,2);
    nd(3) = nodes(iel,3);
    ind   = feeldof(nd,nnel,ndof);
    x1=gcoord(nd(1),1); y1=gcoord(nd(1),2);
    x2=gcoord(nd(2),1); y2=gcoord(nd(2),2);
    x3=gcoord(nd(3),1); y3=gcoord(nd(3),2);
    [f]  = feforctermtriang(x1,y1,x2,y2,x3,y3,gradfcn(iel)*ones(3,1));
    ffs(ind) = ffs(ind) + f;
end
