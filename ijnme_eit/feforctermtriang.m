function [f]=feforctermtriang(x1,y1,x2,y2,x3,y3,fcoef)
% Purpose: forming the forcing term
A=det([1 x1 y1;1 x2 y2;1 x3 y3]);
kexi4=[1/3 1/3 1/3;3/5 1/5 1/5;1/5 3/5 1/5;1/5 1/5 3/5];
omega4=[-27/48 25/48 25/48 25/48];
f=zeros(3,1);
for i=1:4
    x=kexi4(i,1)*x1+kexi4(i,2)*x2+kexi4(i,3)*x3;
    y=kexi4(i,1)*y1+kexi4(i,2)*y2+kexi4(i,3)*y3;
    % using interpolation value
    b = kexi4(i,1)*fcoef(1)+kexi4(i,2)*fcoef(2)+kexi4(i,3)*fcoef(3);
    h1=(x2*y3-x3*y2+(y2-y3)*x+(x3-x2)*y)/A;
    h2=(x3*y1-x1*y3+(y3-y1)*x+(x1-x3)*y)/A;
    h3=(x1*y2-x2*y1+(y1-y2)*x+(x2-x1)*y)/A;
    f(1)=f(1)+omega4(i)*b*h1;
    f(2)=f(2)+omega4(i)*b*h2;
    f(3)=f(3)+omega4(i)*b*h3;
end
f=f*A/2;
%--------------------------------------------------------------------------