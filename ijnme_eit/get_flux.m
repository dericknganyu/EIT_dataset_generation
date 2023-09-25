function [flux] = get_flux(p,e,t,N)
%------------------------------------------------------------------------%
% generate pulsive Neumann b.c. for computing optimal-current pattern    %
% JIN Bangti (kimbtsing@yahoo.com.cn) Feb. 3, 2009                       %
%------------------------------------------------------------------------%

dof = e(1,:);
Ned = length(dof);
pts = p(1:2,dof)';

flux = zeros(Ned,N);
warning off
%  the first one
ind1 = find(pts(:,1)>0);
ind2 = find(pts(:,2)>0);
ind3 = find(atan(pts(:,1)./pts(:,2))>0);
ind4 = find(atan(pts(:,1)./pts(:,2))<=pi/16);
ind3 = intersect(ind3,ind4);
ind = intersect(ind1,intersect(ind2,ind3));
flux(ind,1) = 1;
ind1 = find(pts(:,1)<0);
ind2 = find(pts(:,2)<0);
ind3 = find(atan(pts(:,1)./pts(:,2))>0);
ind4 = find(atan(pts(:,1)./pts(:,2))<=pi/16);
ind3 = intersect(ind3,ind4);
ind = intersect(ind1,intersect(ind2,ind3));
flux(ind,1) = -1;

% the second
ind1 = find(pts(:,1)>0);
ind2 = find(pts(:,2)>0);
ind3 = find(atan(pts(:,1)./pts(:,2))>6*pi/16);
ind4 = find(atan(pts(:,1)./pts(:,2))<=7*pi/16);
ind3 = intersect(ind3,ind4);
ind = intersect(ind1,intersect(ind2,ind3));
flux(ind,2) = 1;
ind1 = find(pts(:,1)<0);
ind2 = find(pts(:,2)<0);
ind3 = find(atan(pts(:,1)./pts(:,2))>6*pi/16);
ind4 = find(atan(pts(:,1)./pts(:,2))<=7*pi/16);
ind3 = intersect(ind3,ind4);
ind = intersect(ind1,intersect(ind2,ind3));
flux(ind,2) = -1;

%  the third one
ind1 = find(pts(:,1)>0);
ind2 = find(pts(:,2)<0);
ind3 = find(atan(pts(:,1)./pts(:,2))+pi>9*pi/16);
ind4 = find(atan(pts(:,1)./pts(:,2))+pi<=10*pi/16);
ind3 = intersect(ind3,ind4);
ind = intersect(ind1,intersect(ind2,ind3));
flux(ind,3) = 1;
ind1 = find(pts(:,1)<0);
ind2 = find(pts(:,2)>0);
ind3 = find(atan(pts(:,1)./pts(:,2))+pi>9*pi/16);
ind4 = find(atan(pts(:,1)./pts(:,2))+pi<=10*pi/16);
ind3 = intersect(ind3,ind4);
ind = intersect(ind1,intersect(ind2,ind3));
flux(ind,3) = -1;

% the forth
ind1 = find(pts(:,1)<0);
ind2 = find(pts(:,2)>0);
ind3 = find(atan(pts(:,1)./pts(:,2))+pi>11*pi/16);
ind4 = find(atan(pts(:,1)./pts(:,2))+pi<=12*pi/16);
ind3 = intersect(ind3,ind4);
ind = intersect(ind1,intersect(ind2,ind3));
flux(ind,4) = 1;
ind1 = find(pts(:,1)>0);
ind2 = find(pts(:,2)<0);
ind3 = find(atan(pts(:,1)./pts(:,2))+pi>11*pi/16);
ind4 = find(atan(pts(:,1)./pts(:,2))+pi<=12*pi/16);
ind3 = intersect(ind3,ind4);
ind = intersect(ind1,intersect(ind2,ind3));
flux(ind,4) = -1;

% the fifth
ind1 = find(pts(:,1)<0);
ind2 = find(pts(:,2)>0);
ind3 = find(atan(pts(:,1)./pts(:,2))+pi>14*pi/16);
ind4 = find(atan(pts(:,1)./pts(:,2))+pi<=15*pi/16);
ind3 = intersect(ind3,ind4);
ind = intersect(ind1,intersect(ind2,ind3));
flux(ind,5) = 1;
ind1 = find(pts(:,1)>0);
ind2 = find(pts(:,2)<0);
ind3 = find(atan(pts(:,1)./pts(:,2))+pi>14*pi/16);
ind4 = find(atan(pts(:,1)./pts(:,2))+pi<=15*pi/16);
ind3 = intersect(ind3,ind4);
ind = intersect(ind1,intersect(ind2,ind3));
flux(ind,5) = -1;
