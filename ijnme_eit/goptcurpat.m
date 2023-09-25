function [flux,fluxeigs] = goptcurpat(KKt,KK1,fem,maxitcur,modcurrpatt,Neig,sigma)
%------------------------------------------------------------------------%
% computing Neig optimal current patterns simultaneously, with the true  %
% conductivity encoded in KKt, and its estimate in KK1                   %
% with new-definition of inner product                                   %
% JIN Bangti (kimbtsing@yahoo.com.cn), March 3, 2009                     %
% last modified on July 7, 2009, calculate also the eigenvalues          %
%------------------------------------------------------------------------%

%------------------------- FEM Parameter ---------------------------%
e = fem.e;
p = fem.p;
t = fem.t;
color = {'k','r','b','g','y','k','k','k','k','k'};

%----------- stiffness matrix & source term ------------------------%
MMi = massmat(fem);       % domain mass matrix
MMb = massmatbdy(fem);    % bdy mass matrix
ffo = zeros(fem.sdof,1);  % homogeneous source
xy = p(1:2,e(1,:));       % coordinate of bdy nodes
Ned = length(fem.bcdof);  % # of bdy nodes
flux = randn(Ned,Neig);
eig_legend = [];
for i=1:Neig
    flux(:,i) = flux(:,i) - mean(flux(:,i));
    if i<=5
        eig_legend = [eig_legend, ['eig # ' num2str(i)]];
    end
end
if Neig == 1
    flux = boundflux(fem,'sin');
end
[Q,R]=qr(flux);
flux = Q(:,1:Neig);
figure(12), 
plotbound(xy,flux,color),legend('eig #1','eig #2','eig #3')

for iter=1:maxitcur
    
    switch modcurrpatt
        
        case 'isacson'   % perform Isacson iteration     
            for j=1:Neig
                [uN_t,uNt_t,uNn_t] = neumansolve(KKt,flux(:,j),ffo,fem);
                [uN1,uNt1,uNn1]    = neumansolve(KK1,flux(:,j),ffo,fem);
                flux(:,j) = uNt_t - uNt1;
            end
            flux = orthtrans(flux,'L2',KKt,fem);
        
        case 'XmY1'
            % from H^{-1/2} to H^1
            for j = 1:Neig
                [uN_t,uNt_t,uNn_t] = neumansolve(KKt,flux(:,j),ffo,fem);  % F_N^\sigma^\ast
                [uD1,uDt1,uDn1]    = dirichsolve(KK1,uNt_t,ffo,fem);
                temp = uDn1;
                [uN_t,uNt_t,uNn_t] = neumansolve(KKt,temp,ffo,fem);
                [uD1,uDt1,uDn1]    = dirichsolve(KK1,uNt_t,ffo,fem);
                flux(:,j) = flux(:,j) - 2*temp + uDn1;
                flux(:,j) = flux(:,j) - mean(flux(:,j));
            end
            flux = orthtrans(flux,'Hm',KK1,fem);
            
        case 'X0Y1'
            % from L^2 to H^1
            for j = 1:Neig
                [uN1,uNt1,uNn1]    = neumansolve(KK1,flux(:,j),ffo,fem);
                [uN_t,uNt_t,uNn_t] = neumansolve(KKt,flux(:,j),ffo,fem);
                utdiff = uNt1 - uNt_t;
                [uD1,uDt1,uDn1]    = dirichsolve(KK1,utdiff,ffo,fem);
                [uN_t,uNt_t,uNn_t] = neumansolve(KKt,uDn1,ffo,fem);
                flux(:,j) = utdiff - uNt_t;
            end
            flux = orthtrans(flux,'L2',KK1,fem);
            
        case 'X0Y0'
            % from L^2 to L^2(\Omega)
            for j = 1:Neig
                % solve for w
                [uN1,uNt1,uNn1]    = neumansolve(KK1,flux(:,j),ffo,fem);
                [uD1,uDt1,uDn1]    = dirichsolve(KK1,uNt1,ffo,fem);
                [uN_t,uNt_t,uNn_t] = neumansolve(KKt,uDn1,ffo,fem);
                [uD1,uDt1,uDn1]    = dirichsolve(KK1,uNt1-uNt_t,ffo,fem);
                cnst = -sum(MMi*uD1)/(2*pi);
                ff = scterm(fem,uD1);
                [w,wt,wn] = neumansolve(KK1,cnst*ones(Ned,1),ff,fem);
                [uD1,uDt1,uDn1] = dirichsolve(KK1,wt,ffo,fem);
                [uN_t,uNt_t,uNn_t] = neumansolve(KKt,uDn1,ffo,fem);
                flux(:,j) = wt - uNt_t;
            end
            flux = orthtrans(flux,'L2',KK1,fem);
          
        case 'XmY0'
            % from H^{-1/2} to L^2(\Omega)
            for j=1:Neig
                % solve for w
                [uN_t,uNt_t,uNn_t] = neumansolve(KKt,flux(:,j),ffo,fem);
                [uD1,uDt1,uDn1]    = dirichsolve(KK1,uNt_t,ffo,fem);
                [uN1,uNt1,uNn1]    = neumansolve(KK1,flux(:,j)-uDn1,ffo,fem);  % for the source term
                cnst = -sum(MMi*uN1)/(2*pi);
                ff = scterm(fem,uN1);  % computing the source term
                [w,wt,wn] = neumansolve(KK1,cnst*ones(Ned,1),ff,fem);
                [uD1,uDt1,uDn1]    = dirichsolve(KK1,wt,ffo,fem);
                temp = uDn1;
                [uN_t,uNt_t,uNn_t] = neumansolve(KKt,uDn1,ffo,fem);
                [uD1,uDt1,uDn1]    = dirichsolve(KK1,uNt_t,ffo,fem);
                flux(:,j) = temp - uDn1;
            end
            flux = orthtrans(flux,'Hm',KK1,fem);
            
        otherwise
            display('Undefined model for optimal current patterns!')
    end

    display([' ....... inner iter No. = ' num2str(iter)])
end

% computing assoc. eigenvalue with lambda_j = \|A^\astA f_j\|/\|f_j\|, 
% the Rayleigh quotient, if required
if nargout == 2
    
    fluxeigs = [];
    
    switch modcurrpatt
        
        case 'isacson'   % perform Isacson iteration     
            for j=1:Neig
                [uN_t,uNt_t,uNn_t] = neumansolve(KKt,flux(:,j),ffo,fem);
                [uN1,uNt1,uNn1]    = neumansolve(KK1,flux(:,j),ffo,fem);
                flxtemp = uNt_t - uNt1;
                fluxeigs = [fluxeigs orthtransnorm(flxtemp,'L2',KKt,fem)];
            end
        
        case 'XmY1'
            % from H^{-1/2} to H^1
            for j = 1:Neig
                [uN_t,uNt_t,uNn_t] = neumansolve(KKt,flux(:,j),ffo,fem);  % F_N^\sigma^\ast
                [uD1,uDt1,uDn1]    = dirichsolve(KK1,uNt_t,ffo,fem);
                temp = uDn1;
                [uN_t,uNt_t,uNn_t] = neumansolve(KKt,temp,ffo,fem);
                [uD1,uDt1,uDn1]    = dirichsolve(KK1,uNt_t,ffo,fem);
                flxtemp = flux(:,j) - 2*temp + uDn1;
                flxtemp = flxtemp - mean(flxtemp);
                fluxeigs = [fluxeigs orthtransnorm(flxtemp,'Hm',KK1,fem)];
            end
            fluxeigs = sqrt(fluxeigs);
            
        case 'X0Y1'
            % from L^2 to H^1
            for j = 1:Neig
                [uN1,uNt1,uNn1]    = neumansolve(KK1,flux(:,j),ffo,fem);
                [uN_t,uNt_t,uNn_t] = neumansolve(KKt,flux(:,j),ffo,fem);
                utdiff = uNt1 - uNt_t;
                [uD1,uDt1,uDn1]    = dirichsolve(KK1,utdiff,ffo,fem);
                [uN_t,uNt_t,uNn_t] = neumansolve(KKt,uDn1,ffo,fem);
                flxtemp = utdiff - uNt_t;
                fluxeigs = [fluxeigs orthtransnorm(flxtemp,'L2',KK1,fem)];
            end
            fluxeigs = sqrt(fluxeigs);
                        
        case 'X0Y0'
            % from L^2 to L^2(\Omega)
            for j = 1:Neig
                % solve for w
                [uN1,uNt1,uNn1]    = neumansolve(KK1,flux(:,j),ffo,fem);
                [uD1,uDt1,uDn1]    = dirichsolve(KK1,uNt1,ffo,fem);
                [uN_t,uNt_t,uNn_t] = neumansolve(KKt,uDn1,ffo,fem);
                [uD1,uDt1,uDn1]    = dirichsolve(KK1,uNt1-uNt_t,ffo,fem);
                cnst = -sum(MMi*uD1)/(2*pi);
                ff = scterm(fem,uD1);
                [w,wt,wn] = neumansolve(KK1,cnst*ones(Ned,1),ff,fem);
                [uD1,uDt1,uDn1] = dirichsolve(KK1,wt,ffo,fem);
                [uN_t,uNt_t,uNn_t] = neumansolve(KKt,uDn1,ffo,fem);
                flxtemp = wt - uNt_t;
                fluxeigs = [fluxeigs orthtransnorm(flxtemp,'L2',KK1,fem)];
            end
            fluxeigs = sqrt(fluxeigs);
            
        case 'XmY0'
            % from H^{-1/2} to L^2(\Omega)
            for j=1:Neig
                % solve for w
                [uN_t,uNt_t,uNn_t] = neumansolve(KKt,flux(:,j),ffo,fem);
                [uD1,uDt1,uDn1]    = dirichsolve(KK1,uNt_t,ffo,fem);
                [uN1,uNt1,uNn1]    = neumansolve(KK1,flux(:,j)-uDn1,ffo,fem);  % for the source term
                cnst = -sum(MMi*uN1)/(2*pi);
                ff = scterm(fem,uN1);  % computing the source term
                [w,wt,wn] = neumansolve(KK1,cnst*ones(Ned,1),ff,fem);
                [uD1,uDt1,uDn1]    = dirichsolve(KK1,wt,ffo,fem);
                temp = uDn1;
                [uN_t,uNt_t,uNn_t] = neumansolve(KKt,uDn1,ffo,fem);
                [uD1,uDt1,uDn1]    = dirichsolve(KK1,uNt_t,ffo,fem);
                flxtemp = temp - uDn1;
                fluxeigs = [fluxeigs orthtransnorm(flxtemp,'Hm',KK1,fem)];
            end
            fluxeigs = sqrt(fluxeigs);
            
        otherwise
            display('Undefined model for optimal current patterns!')
    end
end