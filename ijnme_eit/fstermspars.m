function ff =  fstermspars(fem,inode,t,u)
%---------------------------------------------------------------------%
% generating source term computing Jacobian, observing the sparsity   %
% and piecewise constant nature of the gradient                       %
% JIN Bangti (kimbtsing@yahoo.com.cn) Feb. 6, 2009                    %
%---------------------------------------------------------------------%

% extracting parameters
nnel = fem.nnel;
ndof = fem.ndof;
nel  = fem.nel;
gcoord = fem.gcoord;
nodes  = fem.nodes;

ff=zeros(fem.sdof,1);

for i = 1:3
    ind = find(t(i,:)==inode);
    sigma = zeros(3,1); 
    sigma(i) =1;
    
    for iel=1:length(ind)
        nd(1)=nodes(ind(iel),1);
        nd(2)=nodes(ind(iel),2);
        nd(3)=nodes(ind(iel),3);
        index=feeldof(nd,nnel,ndof);
        x1=gcoord(nd(1),1); y1=gcoord(nd(1),2);
        x2=gcoord(nd(2),1); y2=gcoord(nd(2),2);
        x3=gcoord(nd(3),1); y3=gcoord(nd(3),2);
        K = festifftriang(x1,y1,x2,y2,x3,y3,sigma);
        f = K*u(index);
        ff(index) = ff(index) + f;
    end
end
ff = -ff;