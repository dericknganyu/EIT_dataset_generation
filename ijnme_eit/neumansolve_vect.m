function [u,ut,un]=neumansolve_vect(KK,nbcval,ffo,fem)
    %-----------------------------------------------------------------------%
    % Neumann solver: incorporating Neumann data                            %
    % output: Dirichlet & Neumann data with solution in the domain          %
    % JIN Bangti (kimbtsing@yahoo.com.cn) Feb. 21, 2010                     %
    %-----------------------------------------------------------------------%
    
    % four-point Gauss quadrature abscissas and weights
    kexi4 = [-0.86113631;-0.33998104;0.33998104;0.86113631];
    omega4 = [0.34785485;0.65214515;0.65214515;0.34785485];
    
    e = fem.e;
    gcoord = fem.gcoord;
    
    % incorporating Neumann boundary condition
    dof = fem.bcdof;
    n = size(e,2);  % # of edges
    if ischar(nbcval)
        x = gcoord(dof,1);
        y = gcoord(dof,2);
        fluxcoef = eval(nbcval);
    else
        fluxcoef = nbcval;
    end
 
    % Calculate ind using logical indexing    
    %disp(size(e))
    for i=1:n
        ind(i)=find(e(2,i)==e(1,:));    
    end

    ii   = e(1,:); 
    jj   = e(2,:);
    x1   = gcoord(ii,1);
    y1   = gcoord(ii,2);
    x2   = gcoord(jj,1);
    y2   = gcoord(jj,2);
    xmid = (x1+x2)/2;
    ymid = (y1+y2)/2;   % mid point
    xoffs= (x2-x1)/2;  
    yoffs= (y2-y1)/2;   % offset
    len  = norm([xoffs yoffs]);
    for j=1:4
        x = xmid+kexi4(j)*xoffs;
        y = ymid+kexi4(j)*yoffs;
        flux    = fluxcoef*(1-kexi4(j))/2 + fluxcoef(ind)*(1-(1-kexi4(j))/2);
        ffo(ii) = ffo(ii)+omega4(j)*len*flux*(1-kexi4(j))/2;
        ffo(jj) = ffo(jj)+omega4(j)*len*flux*(1-(1-kexi4(j))/2);
    end


    
    % vanishing mean on bdy
    vec  = zeros(1,fem.sdof);
    vec(dof) = 1/length(dof); 
    KKo = [KK;vec];
    ffo = [ffo;0];
    u = KKo\ffo;
    % u = u-um;
    ut = u(dof);
    un = fluxcoef;