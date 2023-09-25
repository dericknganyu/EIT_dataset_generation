function [MM]=massmatbdy(fem)
%------------------------------------------------------------------------%
% forming the boundary mass matrix (analogous to Robin b.c.)             %
% JIN Bangti (kimbtsing@yahoo.com.cn) Feb 3, 2009                        %
%------------------------------------------------------------------------%

% quadrature constants
kexi4  = [-0.86113631;-0.33998104;0.33998104;0.86113631];
omega4 = [0.34785485;0.65214515;0.65214515;0.34785485];

% extracting fe parameters
e      = fem.e;
gcoord = fem.gcoord;
Ned    = length(fem.bcdof);

MM     = sparse(Ned,Ned);
bcdof1 = e(1,:);
bcdof2 = e(2,:);
ind    = zeros(Ned,1);
for i=1:Ned
    ind(i)=find(e(2,i)==e(1,:));    
end
for i=1:Ned
    ii = bcdof1(i); 
    jj = bcdof2(i);
    x1 = gcoord(ii,1);
    y1 = gcoord(ii,2);
    x2 = gcoord(jj,1);
    y2 = gcoord(jj,2);
    xoffs = (x2-x1)/2; 
    yoffs = (y2-y1)/2;
    len   = norm([xoffs yoffs]);
    for j=1:4
        MM(i,i)           = MM(i,i)+omega4(j)*len*((kexi4(j)+1)/2)^2;
        MM(i,ind(i))      = MM(i,ind(i))+omega4(j)*len*((kexi4(j)+1)/2)*(1-(kexi4(j)+1)/2);
        MM(ind(i),i)      = MM(ind(i),i)+omega4(j)*len*((kexi4(j)+1)/2)*(1-(kexi4(j)+1)/2);
        MM(ind(i),ind(i)) = MM(ind(i),ind(i))+omega4(j)*len*(1-(kexi4(j)+1)/2)^2;
    end
end