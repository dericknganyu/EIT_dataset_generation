function [ffo]=fenbc(ff,gcoord,nbcdof,fluxcoef)
% Purpose: treating the Neumann condition
n=size(nbcdof,1)-1;
ffo=ff;
% four-point Gauss quadrature abscissas and weights
kexi4 = [-0.86113631;-0.33998104;0.33998104;0.86113631];
omega4 = [0.34785485;0.65214515;0.65214515;0.34785485];
for i=1:n
    ii=nbcdof(i);
    jj=nbcdof(i+1);
    x1=gcoord(ii,1); y1=gcoord(ii,2);
    x2=gcoord(jj,1); y2=gcoord(jj,2);
    xmid=(x1+x2)/2; ymid=(y1+y2)/2;   % mid point of the line segment
    xoffs=(x2-x1)/2; yoffs=(y2-y1)/2; % offset from the mid point
    len=norm([xoffs yoffs]);
    for j=1:4
        x=xmid+kexi4(j)*xoffs;
        y=ymid+kexi4(j)*yoffs;
        flux=fluxcoef(i)*(1-kexi4(j))/2+fluxcoef(i+1)*(1-(1-kexi4(j))/2);
        ffo(ii)=ffo(ii)+omega4(j)*len*flux*(1-kexi4(j))/2;
        ffo(jj)=ffo(jj)+omega4(j)*len*flux*(1-(1-kexi4(j))/2);
    end
end
%--------------------------------------------------------------------------