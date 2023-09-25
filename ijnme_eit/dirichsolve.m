function [u,ut,un]=dirichsolve(KK,dbcval,ffo,fem,sigma)
%------------------------------------------------------------------------%
% solving the Dirichlet problem, output: Dirichlet, Neumann values and   %
% solution in the domain                                                 %
% JIN Bangti (kimbtsing@yahoo.com.cn), Feb. 1, 2009                      %
%------------------------------------------------------------------------%

e      = fem.e;
gcoord = fem.gcoord;
dbcdof = fem.bcdof;
KKo    = KK;
sdof   = fem.sdof;
n      = length(dbcdof);
for i = 1:n
    ii = dbcdof(i);
    for j = 1:sdof
        KKo(ii,j) = 0;
    end
    KKo(ii,ii) = 1;
end
if ischar(dbcval)
    x = gcoord(dbcdof,1);
    y = gcoord(dbcdof,2);
    dbcval = eval(dbcval);
end
ffo(dbcdof) = dbcval;
u = KKo\ffo;
ut = dbcval;
if nargin==4
    un = get_neumann(u,fem);
else
    un = get_neumann(u,fem,sigma);
end