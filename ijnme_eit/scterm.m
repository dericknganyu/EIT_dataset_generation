function ffo =  scterm(fem,fcoef)
%---------------------------------------------------------------------%
% generate the source term given analytical or numerical sources      %
% JIN Bangti (kimbtsing@yahoo.com.cn) Feb. 1, 2009                    %
%---------------------------------------------------------------------%

gcoord= fem.gcoord;
nodes = fem.nodes;

ffo=zeros(fem.sdof,1);
if ischar(fcoef)
    x = gcoord(:,1);
    y = gcoord(:,2);
    fcoef = eval(fcoef);
    clear x y
end
for iel=1:fem.nel
    nd(1) = nodes(iel,1);
    nd(2) = nodes(iel,2);
    nd(3) = nodes(iel,3);
    ind = feeldof(nd,fem.nnel,fem.ndof);
    x1 = gcoord(nd(1),1); 
    y1 = gcoord(nd(1),2);
    x2 = gcoord(nd(2),1);
    y2 = gcoord(nd(2),2);
    x3 = gcoord(nd(3),1);
    y3 = gcoord(nd(3),2);
    f  = feforctermtriang(x1,y1,x2,y2,x3,y3,fcoef(ind));
    ffo(ind) = ffo(ind) + f;
end
