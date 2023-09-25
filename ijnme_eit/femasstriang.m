function [m]=femasstriang(x1,y1,x2,y2,x3,y3,acoef)
%-------------------------------------------------------------------------%
% computing element mass matrix for triangular elements                   %
% JIN Bangti (kimbtsing@yahoo.com.cn) Feb. 2, 2009                        %
%-------------------------------------------------------------------------%

A=det([1 x1 y1;1 x2 y2;1 x3 y3]);
kexi4=[1/3 1/3 1/3;3/5 1/5 1/5;1/5 3/5 1/5;1/5 1/5 3/5];
omega4=[-27/48 25/48 25/48 25/48];
m=zeros(3,3);
for i=1:4
    x=kexi4(i,1)*x1+kexi4(i,2)*x2+kexi4(i,3)*x3;
    y=kexi4(i,1)*y1+kexi4(i,2)*y2+kexi4(i,3)*y3;
    a = kexi4(i,1)*acoef(1)+kexi4(i,2)*acoef(2)+kexi4(i,3)*acoef(3);
   
    h1=((x2*y3-x3*y2)+(y2-y3)*x+(x3-x2)*y)/A;
    h2=((x3*y1-x1*y3)+(y3-y1)*x+(x1-x3)*y)/A;
    h3=((x1*y2-x2*y1)+(y1-y2)*x+(x2-x1)*y)/A;
   
    m(1,1)=m(1,1)+omega4(i)*a*h1*h1;
    m(1,2)=m(1,2)+omega4(i)*a*h1*h2;
    m(1,3)=m(1,3)+omega4(i)*a*h1*h3;
    m(2,2)=m(2,2)+omega4(i)*a*h2*h2;
    m(2,3)=m(2,3)+omega4(i)*a*h2*h3;
    m(3,3)=m(3,3)+omega4(i)*a*h3*h3;
end
m(2,1)=m(1,2);
m(3,1)=m(1,3);
m(3,2)=m(2,3);
m=m*A/2;