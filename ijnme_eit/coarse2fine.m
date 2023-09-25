function fluxf = coarse2fine(flux,fem,femf)
%------------------------------------------------------------------------%
% interpolate the flux on coarse mesh to fine mesh for computing exact   %
% measurements                                                           %
% JIN Bangti (kimbtsing@yahoo.com.cn) February 22, 2009                  %
%------------------------------------------------------------------------%

% extract fem parameters
p = fem.p;
e = fem.e;
p1= femf.p;
e1= femf.e;


if length(e(1,:))==length(e1(1,:))
    fluxf = flux;
    return;
end

ptc = p(1:2,e(1,:));
for i=1:length(e1(1,:))
    cpt = p1(1:2,e1(1,i));  % coordinate of current node
    dist = zeros(length(e(1,:)),1);
    for j=1:length(e(1,:))
        dist(j) = norm(ptc(1:2,j)-cpt);
    end
    [Y,I] = sort(dist);
    if (Y(1)==0)
        fluxf(i) = flux(I(1));
    else
        % linear interpolation from nearest two-point
        fluxf(i) = (flux(I(1))*Y(2)+flux(I(2))*Y(1))/(Y(1)+Y(2));
    end
end