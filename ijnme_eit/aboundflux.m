function [flux] = aboundflux(fem,type,N)
    %------------------------------------------------------------------------%
    % generate pulsive Neumann b.c. for computing optimal-current pattern    %
    % JIN Bangti (kimbtsing@yahoo.com.cn) Feb. 3, 2009                       %
    %------------------------------------------------------------------------%
    
    p   = fem.p;
    %e   = fem.e;
    dof = fem.bcdof;
    Ned = length(dof);
    pts = p(1:2,dof)';
    
    if (nargin<=1)
        type = 'sin';
        N = 1;
    elseif (nargin<=2)
        N = 1;
    end
    
    switch lower(type)
        
        case 'pulse'  % impulsive flux
            flux = zeros(Ned,N);
            for j=1:N
                ind1 = find(pts(:,1)>0); 
                ind2 = find(pts(:,2)>0); 
                ind3 = find(atan(pts(:,1)./pts(:,2))>(j-1)*pi/16);
                ind4 = find(atan(pts(:,1)./pts(:,2))<=j*pi/16);
                ind3 = intersect(ind3,ind4);
                ind = intersect(ind1,intersect(ind2,ind3));
                flux(ind,j) = 1;
                ind1 = find(pts(:,1)<0);
                ind2 = find(pts(:,2)<0);
                ind3 = find(atan(pts(:,1)./pts(:,2))>(j-1)*pi/16);
                ind4 = find(atan(pts(:,1)./pts(:,2))<=j*pi/16);
                ind3 = intersect(ind3,ind4);
                ind = intersect(ind1,intersect(ind2,ind3));
                flux(ind,j) = -1;
            end
            
        case 'sin'   % sinusoidal flux
            flux = zeros(Ned,N);
            for j=1:N
                if j <=16
                    ind = find(pts(:,1)>0);
                    flux(ind,j) = sin(j*(atan(pts(ind,2)./pts(ind,1))));
                    ind = find(pts(:,1)<0);
                    flux(ind,j) = sin(j*(atan(pts(ind,2)./pts(ind,1))+pi));
                elseif j>16
                    ind = find(pts(:,1)>0);
                    flux(ind,j) = cos((j-16)*(atan(pts(ind,2)./pts(ind,1))));
                    ind = find(pts(:,1)<0);
                    flux(ind,j) = cos((j-16)*(atan(pts(ind,2)./pts(ind,1))+pi));
                end       
            end
            flux = flux*(1/sqrt(pi)); 
            
        case 'piecewise'   % piecewise constant flux
            flux = ones(Ned,1);
            ind = find((pts(:,2)-0).^2<1e-10);
            flux(ind) = 0;
            ind = find(pts(:,2)<0);
            flux(ind) = -flux(ind);
            
        otherwise
            error('Undefined input current pattern!')
    end