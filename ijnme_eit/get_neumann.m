function dudn = get_neumann(u,fem,sigma)
%-------------------------------------------------------------------------%
% calculate the Neumann data from solution u on the mesh (p,e,t) by finite%
% difference method                                                       %
% JIN Bangti (kimbtsing@yahoo.com.cn), Feb. 1, 2009                       %
%-------------------------------------------------------------------------%

p = fem.p;
e = fem.e;
t = fem.t; 

nbcdof = e(1,:);
n      = length(nbcdof);
dudn   = zeros(n,1);
xy     = p(1:2,nbcdof);
dun    = u(nbcdof);
for i=1:n
    
    a1=e(1,i);    
    a2=e(2,find(e(1,:)==a1));
    a3=e(1,find(e(2,:)==a1));
    % we have to search for the triangle which has only one point on the
    % boundary for the finite difference approximation of the derivative
    
    for k=1:3
        ix=find(t(k,:)==a1);
        for l=1:length(ix)
            if t(mod(k,3)+1,ix(l))~=a2 && t(mod(k,3)+1,ix(l))~=a3 && t(mod(k+1,3)+1,ix(l))~=a2 && t(mod(k+1,3)+1,ix(l))~=a3
                six=ix(l);
                sk=k;
                break;
            end
        end
    end
    
    d1=p(:,a1);
    d2=p(:,t(mod(sk,3)+1,six));
    d3=p(:,t(mod(sk+1,3)+1,six));
    % now that we know the three points of the triangle we determine where
    % the normal vector goes thru the triangle to gather parameters for a
    % bilinear interpolation
    
%     phi1=acos( ((d2-d1)'*d1)/(norm(d2-d1)*norm(d1)));
%     phi2=acos( ((d3-d2)'*d1)/(norm(d3-d2)*norm(d1)));
%     l1=norm(d2-d1).*sin(phi1)./sin(phi2);
    l1=norm(d2-d1).*sqrt(1-(((d2-d1)'*d1)./(norm(d2-d1)*norm(d1)))^2)./sqrt(1-(((d3-d2)'*d1)/(norm(d3-d2)*norm(d1)))^2);
    
    
%     phi1=acos( ((d3-d1)'*d1)/(norm(d3-d1)*norm(d1)));
%     phi2=pi-phi2;
%     l2=norm(d3-d1).*sin(phi1)./sin(phi2);
    l2=norm(d3-d1).*sqrt(1-(((d3-d1)'*d1)./(norm(d3-d1)*norm(d1)))^2)./sqrt(1-(((d2-d3)'*d1)/(norm(d2-d3)*norm(d1)))^2);

    d4=d2+l1./(l1+l2).*(d3-d2);
    l3=norm(d4-d1);
    
    % now we know enough to do a bilinear interpolation and compute the
    % normal derivative with a difference quotient.
    if nargin==2
        dudn(i)=(dun(i)-((l1.*u(t(mod(sk+1,3)+1,six))+l2.*u(t(mod(sk,3)+1,six)))./(l1+l2)))./l3;
    else
        dudn(i)=(dun(i)-((l1.*u(t(mod(sk+1,3)+1,six))+l2.*u(t(mod(sk,3)+1,six)))./(l1+l2)))./l3;
        tmp = (l1.*sigma(t(mod(sk+1,3)+1,six))+l2.*sigma(t(mod(sk,3)+1,six)))./(l1+l2)
        dudn(i) = dudn(i).*tmp;
    end
    
end