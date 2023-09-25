%-------------------------------------------------------------------------%
% sparsity regularization for domain-based data-fitting with L1-H1 regu.  %
% by primal-dual active set strategy (fixed point algorithm)              %
% JIN Bangti (kimbtsing@yahoo.com.cn) April 15, 2010                      %
%-------------------------------------------------------------------------%
% clear, close all
% 
% %%
% %-------------------------- mesh generation ---------------------------%
% dl=circle;
% [p,e,t]=initmesh(dl);
% [p,e,t]=refinemesh(dl,p,e,t);
% [p1,e1,t1]=refinemesh(dl,p,e,t);
% [p1,e1,t1]=refinemesh(dl,p1,e1,t1);
% % [p,e,t]=refinemesh(dl,p,e,t); figure, pdemesh(p,e,t), axis square, axis
% % off
% 
% % fem parameters for coarse mesh
% nnel=3;              % # nodes per element
% ndof=1;              % # dof per node
% gcoord = p(1:2,:)';  % coordinates of nodes 
% nnode = size(p,2);   % # of system nodes
% sdof=nnode*ndof;     % system dofs
% nodes = t(1:3,:)';   % triangle element
% nel = size(nodes,1); % # of elements
% bcdof = e(1,:);
% % fem parameters for fine mesh
% gcoordf = p1(1:2,:)';  % coordinates of nodes 
% nnodef = size(p1,2);   % # of system nodes
% sdoff = nnodef*ndof;     % system dofs
% nodesf = t1(1:3,:)';   % triangle element
% nelf = size(nodesf,1); % # of elements
% bcdoff = e(1,:);
% 
% % conductivity and source term
% sigma0 = 'exp(x-x)';     % homogeneous background
% % background + inhomogeneities
% sigmat = '1 + 5*heavis(0.09-x.^2-(y+0.65).^2) + 5*heavis(0.09-x.^2-(y-0.65).^2) + 5*heavis(0.09-(x+0.65).^2-(y).^2) + 5*heavis(0.09-(x-0.65).^2-(y).^2)';
% % sigmat = '1 +
% % 5*(heavis(0.09-x.^2-(y-.65).^2)+heavis(0.09-x.^2-(y+.65).^2))'; sigmat =
% % '1+5*heavis(0.09-(x-0.5).^2-(y-0.45).^2)'; sigmat =
% % '1+5*heavis(0.09-(x-0.1).^2-(y).^2)';
% x = gcoord(:,1); y = gcoord(:,2);
% sigmatt = eval(sigmat);
% 
% %----------- stiffness matrix & source term ------------------------%
% KK0 = stiffmat(sdof,nnel,ndof,nel,gcoord,nodes,sigma0);
% KKt = stiffmat(sdof,nnel,ndof,nel,gcoord,nodes,sigmat);    % true conduct.
% KKs = stiffmat(sdof,nnel,ndof,nel,gcoord,nodes,'exp(x-x)');
% MM = massmat(sdof,nnel,ndof,nel,gcoord,nodes);             % used for computing domain integrals
% % KKs = KKs + MM;        % required by smoothing step of gradient, i.e.
% % Sobolev gradient
% ffo = zeros(sdof,1);   % homogeneous source
% xy = p(1:2,e(1,:));    % coordinate of bdy nodes
% Ned = length(e(1,:));  % # of bdy nodes
% KKtf = stiffmat(sdoff,nnel,ndof,nelf,gcoordf,nodesf,sigmat);
% 
% %%
% %------------------------ data generation -------------------------%
% % computing the optimal input current
% Neig = 5;           % # of optimal current
% maxitcur = 50;     % maximum # of iterations for computing optimal currents
% modcurrpatt =  'isacson'; %'isacson'; % 'XmY0';
% [flux,wt] = goptcurpat(KKt,KK0,p,e,t,maxitcur,modcurrpatt,Neig);
% figure(5), plotbound(xy,flux)
% wt 
% 
% temp = zeros(size(flux));
% for j = 1 : Neig
%     fluxf = coarse2fine(flux(:,j),p,e,t,p1,e1,t1);    % interpolation to fine mesh
%     [uNt,uNtt,uNnt] = neumansolve(KKtf,fluxf,zeros(sdoff,1),gcoordf,p1,e1,t1);
%     temp(:,j) = uNt(e(1,:));    % electrical potential measurements
% end
% std = max(max(abs(temp)));
% randn('state',0)
% NL = 1e-2; % NL = 5e-2;
% % temp = temp + NL*std*randn(size(temp));  % generate noisy data
% 
% %----------------------- inverse iteration ------------------------%
% sigma = ones(sdof,1);      % background conductivity
% dsig = zeros(sdof,1);      % initial guess of inhomogeneities
% % sigmatmp = '1+5*heavis(0.09-(x+0.2).^2-(y).^2)';
% % sigmatmp = eval(sigmatmp);
% % dsig = sigmatt-1;
% % dsig = dsig + 5*rand(sdof,1);
% maxiter = 600;             % maximum # of iterations for sparse updating for each optimal current
% eta = 3e-4;   % regu. param. for L^1
% beta= 1e-3;   % regu. param. for H^1
%  
% %%
% % % pre-calcuation of needed variables
% w = eta*ones(sdof,1);
% [gradfcn,objfcn,regval] = grad_fcn(Neig,MM,sigma,dsig,p,e,t,flux,temp,wt,w);
% lamda = zeros(sdof,1);  % variable for complementarity condition
% c = 1e9;                % damping parameter for ssn
% dt = 2e-1;              % time-step size
% lamda0 = lamda;
% dsig0  = dsig;

% %%
% for iter = 1:maxiter
%     
%     % determine active and inactive sets
%     pp = abs(lamda0 + c*dsig0);
%     IndI = find(pp<=1);
%     IndA = find(pp>1);
%     size(IndA)
%     
%     % update u on inactive set
%     dsig = zeros(sdof,1);
%     
%     % smoothed gradient
%     ffs = sctermpw(sdof,nnel,ndof,nel,gcoord,nodes,gradfcn);
%     [gradfcn] = dirichsolve(KKs+MM,zeros(Ned,1),ffs,gcoord,p,e,t);
%     
%     % update u on inactive set
%     tmp = max(1./(pp-1),0);
%     MMt = massmat(sdof,nnel,ndof,nel,gcoord,nodes,tmp);
%     A = MM(IndA,IndA) + dt*beta*KKs(IndA,IndA) + dt*eta*c*MMt(IndA,IndA);
%     % A = MM(IndA,IndA) + dt*beta*KKs(IndA,IndA) + dt*eta*c*diag(1./(pp(IndA)-1))*MM(IndA,IndA);
%     % ffs = sctermpw(sdof,nnel,ndof,nel,gcoord,nodes,gradfcn);
%     b = MM*dsig0 - dt*gradfcn;
%     b = b(IndA);
%     dsig(IndA) = A\b;
%     
%     % update lambda
%     lamda(IndA) = c*dsig(IndA)./(pp(IndA)-1);
%     res = (dsig-dsig0)/dt - beta*KKs*dsig + ffs;
%     lamda(IndI) = - res(IndI)/eta;
%     
%     [gradfcn,objfcn,regval] = grad_fcn(Neig,MM,sigma,dsig,p,e,t,flux,temp,wt,w);
% 
%     % record & display results
%     err = (dsig+sigma-sigmatt)'*MM*(dsig+sigma-sigmatt);    % squared error
%     
%     figure(3), clf
%     subplot(2,1,1), pdeplot(p,e,t,'xydata',dsig), drawnow
%     subplot(2,1,2), pdecont(p,t,dsig,20), drawnow % colorbar
%     lamda0 = lamda;
%     dsig0  = dsig;
%     display(sprintf('iter=%0.5g',iter))
%     if rem(iter,20)==0 & iter<=300
%         eta = eta/0.8;
%         beta = beta*0.8;
%     end
% end


%%
dsig = zeros(sdof,1);
dt = 5e-1;
eta = 3e-4;   % regu. param. for L^1
beta= 1e-3;   % regu. param. for H^1
for iter = 1:maxiter
    
    % smoothed gradient
    ffs = sctermpw(sdof,nnel,ndof,nel,gcoord,nodes,gradfcn);
    % [gradfcn] = dirichsolve(KKs+MM,zeros(Ned,1),ffs,gcoord,p,e,t);
    gradfcn = ffs;
    
    % update u 
    MMt = massmat(sdof,nnel,ndof,nel,gcoord,nodes,1./max(abs(dsig),1e-10));
    % A = speye(549)/dt+beta*KKs+eta*spdiags(1./max(abs(dsig),1e-10),0,549,549);
    A = MM/dt+beta*KKs+eta*MMt;
    b = MM*dsig0/dt-gradfcn;
    dsig = neumansolve(A,zeros(Ned,1),b,gcoord,p,e,t);
    [gradfcn,objfcn,regval] = grad_fcn(Neig,MM,sigma,dsig,p,e,t,flux,temp,wt,w);

    % record & display results
    err = (dsig+sigma-sigmatt)'*MM*(dsig+sigma-sigmatt);
    figure(3), clf
    % subplot(2,1,1), 
    pdeplot(p,e,t,'xydata',dsig), drawnow, colormap('hsv')
    % subplot(2,1,2), pdecont(p,t,dsig,20), drawnow % colorbar
    lamda0 = lamda;
    dsig0  = dsig;
    display(sprintf('iter=%0.5g',iter))
    if rem(iter,20)==0
        eta = eta/0.8;
        beta = beta*0.8;
    end
end