function [k]=festifftriang(x1,y1,x2,y2,x3,y3,acoef)
%-------------------------------------------------------------------------%
% forming element stiff matrix                                            %
% JIN Bangti (kimbtsing@yahoo.com.cn) Feb. 1, 2009                        %
%-------------------------------------------------------------------------%

A=det([1 x1 y1;1 x2 y2;1 x3 y3]);
kexi4=[1/3 1/3 1/3;3/5 1/5 1/5;1/5 3/5 1/5;1/5 1/5 3/5];
omega4=[-27/48 25/48 25/48 25/48];
k=zeros(3,3);
for i=1:4
    x=kexi4(i,1)*x1+kexi4(i,2)*x2+kexi4(i,3)*x3;
    y=kexi4(i,1)*y1+kexi4(i,2)*y2+kexi4(i,3)*y3;    
    % a=eval(acoef);
    % computing the conductivity using interpolation value
    a = kexi4(i,1)*acoef(1)+kexi4(i,2)*acoef(2)+kexi4(i,3)*acoef(3);
    
    k(1,1)=k(1,1)+omega4(i)*((y2-y3)*(y2-y3)+(x3-x2)*(x3-x2))*a/(A^2);
    k(1,2)=k(1,2)+omega4(i)*((y2-y3)*(y3-y1)+(x3-x2)*(x1-x3))*a/(A^2);
    k(1,3)=k(1,3)+omega4(i)*((y2-y3)*(y1-y2)+(x3-x2)*(x2-x1))*a/(A^2);
    k(2,2)=k(2,2)+omega4(i)*((y3-y1)*(y3-y1)+(x1-x3)*(x1-x3))*a/(A^2);
    k(2,3)=k(2,3)+omega4(i)*((y3-y1)*(y1-y2)+(x1-x3)*(x2-x1))*a/(A^2);
    k(3,3)=k(3,3)+omega4(i)*((y1-y2)*(y1-y2)+(x2-x1)*(x2-x1))*a/(A^2);
end
k(2,1)=k(1,2);k(3,1)=k(1,3);k(3,2)=k(2,3); k=k*A/2;
%--------------------------------------------------------------------------
