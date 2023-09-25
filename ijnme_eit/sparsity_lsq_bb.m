%-------------------------------------------------------------------------%
% sparsity reconstruction with least-squares fidelity via soft-shrinkage  %
% Created on Sept. 18 by Bangti JIN (btjin@math.tamu.edu)                 %
%-------------------------------------------------------------------------%
clear, close all

%-------------------------- mesh generation ---------------------------%
dl=circle;
[p,e,t]=initmesh(dl);
[p,e,t]=refinemesh(dl,p,e,t);
[p1,e1,t1]=refinemesh(dl,p,e,t);
[p1,e1,t1]=refinemesh(dl,p1,e1,t1);
% [p1,e1,t1] = refinemesh(dl,p1,e1,t1);
% figure, pdemesh(p,e,t), axis square, axis off
size(p)
size(e)
size(t)
% structure fem: parameters for coarse mesh
fem.nnel=3;                      % # nodes per element
fem.ndof=1;                      % # dof per node
fem.gcoord = p(1:2,:)';          % coordinates of nodes 
fem.nnode = size(p,2);           % # of system nodes
fem.sdof=fem.nnode*fem.ndof;     % system dofs
fem.nodes = t(1:3,:)';           % triangle element
fem.nel = size(fem.nodes,1);     % # of elements
fem.bcdof = e(1,:);
fem.p = p;
fem.e = e;
fem.t = t;
% structure femf: parameters for fine mesh
femf.nnel = 3;                   % # nodes per element
femf.ndof = 1;                   % # dof per node
femf.gcoord = p1(1:2,:)';        % coordinates of nodes 
femf.nnode = size(p1,2);         % # of system nodes
femf.sdof = femf.nnode*femf.ndof;% system dofs
femf.nodes = t1(1:3,:)';         % triangle element
femf.nel = size(femf.nodes,1);   % # of elements
femf.bcdof = e1(1,:);
femf.p = p1;
femf.e = e1;
femf.t = t1;

%-------------------------- example setup ----------------------------%
exam = 6;
oc   = 1;
sig0 = 'exp(x-x)';   % homogeneous background
switch exam
    case 1
        sigmat = '1+5*heavis(0.09-(x-0.5).^2-(y-0.45).^2)';
    case 2
        sigmat = '1 + 5*heavis(0.09-x.^2-(y+0.6).^2) + 5*heavis(0.09-x.^2-(y-0.6).^2) + 5*heavis(0.04-(x+0.76).^2-(y).^2) + 5*heavis(0.04-(x-0.75).^2-(y).^2)';
    case 3
        sigmat = '1 + 5*(heavis(0.64-x.^2-(y-.15).^2)-heavis(0.36-x.^2-(y-.15).^2))';
    case 4
        sigmat = '1+2*heavis(0.09-x.^2-(y+.6).^2) - .6*heavis(0.09-x.^2-(y-.6).^2)';
    case 5
        sigmat = '1+.5*heavis(y) + 5*heavis(0.09-x.^2-(y+0.6).^2) + 5*heavis(0.09-x.^2-(y-0.6).^2)';
    case 6
        sigmat = '1+5*heavis(0.09-(x-0.5).^2-(y-0.45).^2)';
    otherwise
        error('Undefined example!')
end
x = fem.gcoord(:,1);
y = fem.gcoord(:,2);
sigmatt = eval(sigmat);
if exam == 6
    rand('state',0)
    sigmatt = sigmatt + (2*rand(size(sigmatt))-1)*0.2;
end

%------------------------ input current ------------------------%
if oc == 0 %IGNORE
    % load the flux from stored data
    Neig = 5;
    switch exam
        case 1
            load flux_1circ.mat;
        case 2
            load flux_4circ.mat;
        case 3
            load flux_torus.mat;
        case 4
            load flux_pmcirc.mat
    end
    flux = flux(:,1:Neig);
    wt   = wt(1:Neig);
elseif oc == 1
    % sine
    Neig = 5;
    flux = boundflux(fem,'sin',Neig);
    wt   = ones(Neig,1);
else %IGNORE
    % compute the flux on the fly
    KK0  = stiffmat(fem,sig0);     
    KKt  = stiffmat(fem,sigmat);
    Neig = 10;                  % # of optimal current
    maxitcur    = 200;         % maximum # of iters
    modcurrpatt =  'isacson';  %'isacson'; % 'XmY0';
    randn('state',0)
    [flux,wt]   = goptcurpat(KKt,KK0,fem,maxitcur,modcurrpatt,Neig);
    % figure(5), plotbound(xy,flux)    
end

% simulated data
if ((exam>=1) & (exam<=5))
    KKtf= stiffmat(femf,sigmat);% stiff matrix - fine mesh
    temp = zeros(size(flux));
    for j = 1 : Neig
        fluxf = coarse2fine(flux(:,j),fem,femf);    % interpolation to fine mesh
        [uNt,uNtt,uNnt] = neumansolve(KKtf,fluxf,zeros(femf.sdof,1),femf);
        temp(:,j) = uNt(e(1,:));                    % potential measurements
    end
    
elseif exam==6
    KKt = stiffmat(fem,sigmatt);% stiff matrix - fine mesh 
    temp = zeros(size(flux));
    for j = 1 : Neig
        [uNt,uNtt,uNnt] = neumansolve(KKt,flux(:,j),zeros(fem.sdof,1),fem);
        temp(:,j) = uNt(e(1,:));                    % potential measurements
    end
end
std = max(max(abs(temp)));
randn('state',0)
NL = 3e-2;
temp = temp + NL*std*randn(size(temp));
    
%% sparse reconstruction
%----------------------- inverse iteration ------------------------%
KKs = stiffmat(fem,'exp(x-x)');
MM  = massmat(fem);
KKs = KKs + MM;             % Sobolev grad
MMb = massmatbdy(fem);
invpm.Neig = Neig;
invpm.MM   = MM;
invpm.KKs  = KKs;
invpm.MMb  = MMb;
invpm.sig0 = ones(fem.sdof,1);     % bgd cond.
if exam == 5
    x = fem.gcoord(:,1); 
    y = fem.gcoord(:,2);
    invpm.sig0 = eval('1+.5*heavis(y)');
end
invpm.flux = flux;
invpm.temp = temp;
invpm.wt   = wt;
invpm.maxit= 500;

% % pre-calcuation of needed variables
tic
dim = 10;
etaens = logspace(-5,-2.5,dim); 
for k = 6
    invpm.eta = etaens(k);
    [dsig,hist] = recon_sparsity(invpm,fem,0,sigmatt);
    hf = figure(5);
    display(sprintf('save reconstruction ...  %d \n',k))
    pdeplot(fem.p,fem.e,fem.t,'xydata',invpm.sig0+dsig), drawnow, colormap('jet')
    pfadplot='recon/sparse/';
    filename = [pfadplot 'sparse-exp' num2str(exam) '-NL' num2str(100*NL) '-alpha' num2str(k)];
    axis image, axis off, colormap('jet')
    set(colorbar,'FontSize',16), drawnow
    % print(5, '-dpng' ,filename)
    % saveas(5,filename,'fig')
end
toc
    
%%
% reconstruction results with standard method
etaens = logspace(-4,-1,dim);
for k = 5
    sig0  = ones(fem.sdof,1);
    if exam == 5
        x = fem.gcoord(:,1); 
        y = fem.gcoord(:,2);
        sig0 = eval('1+.5*heavis(y)');
    end
    hist  = [];
    alpha = etaens(k);
    dsig0  = zeros(fem.sdof,1);
    tic
    for iter = 1:2
        % linearized inverse model
        [J,FN0d] = jacobian(fem,sig0,flux);
        b = temp(:) - FN0d;
        
        % Tikhonov regularized solution
        n = size(J,2);
        ind = setdiff([1:fem.sdof],e(1,:));
        x = (J'*J+alpha*KKs(ind,ind))\(J'*b-alpha*KKs(ind,ind)*dsig0(ind));
        dsig0 = zeros(fem.sdof,1);
        dsig0(ind) = x;
        figure(1), pdeplot(p,e,t,'xydata',sig0+dsig0), colormap('jet'), drawnow,
        sig0 = sig0 + dsig0;
    end
    toc
    err = sqrt((sig0-sigmatt)'*invpm.MM*(sig0-sigmatt))
    % save the figures
    hf = figure(10);
    pfadplot='recon/smooth/';
    filename = [pfadplot 'smooth-exp' num2str(exam) '-NL' num2str(100*NL) '-alpha' num2str(k)];
    pdeplot(fem.p,fem.e,fem.t,'xydata',sig0)
    axis image, axis off, colormap('jet')
    set(colorbar,'FontSize',16), drawnow
    %print(10, '-dpng' ,filename)
    %saveas(10,filename,'fig')
end