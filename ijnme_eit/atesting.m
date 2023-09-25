%-------------------------------------------------------------------------%
% sparsity reconstruction with least-squares fidelity via soft-shrinkage  %
% Created on Sept. 18 by Bangti JIN (btjin@math.tamu.edu)                 %
%-------------------------------------------------------------------------%
clear, close all
%load dataset
texture = false;

if ~texture
    PATH_read = '/pvfs2/Derick/EIT/Mine/data';%'/pvfs2/Derick/EIT/Mine/data'; %'/localdata/Derick/EIT/Mine/data';
    name = '100_samples__max_Inclusions_3__2023-09-25-10-58-00';
% else
%     PATH_read = '/pvfs2/Derick/EIT/Mine/data_texture'; %'/localdata/Derick/EIT/Mine/data';
%     name = '1_samples__max_Inclusions_3__2023-08-14-16-46-00';
end

inclusions = load(sprintf('%s/%s/inclusions',PATH_read,name));

%samp = 6;
n_current = 16;
hval = 0.3;
max_samp = 100;

%size(aa)
%[batchSize, Nmax, n_nodes] = size(inputConductivity);

mypde = createpde();
geometryFromEdges(mypde,@circleg);

mypde_fwd = createpde();
geometryFromEdges(mypde_fwd,@circleg);

hmax = hval;%0.03;% 0.013;%38%0.01;%4;%0.038;%0.005;
hmin = hval;%0.03;% 0.013;%0.01;%4;
mesh = generateMesh(mypde,"Hmax",hmax,"Hmin",hmin,"GeometricOrder","linear");

mesh_fwd = generateMesh(mypde_fwd,"Hmax",0.03,"Hmin",0.03,"GeometricOrder","linear");        
Nodes_fwd = mesh_fwd.Nodes;
x_fwd = Nodes_fwd(1,:);
y_fwd = Nodes_fwd(2,:);



%-------------------------- mesh generation ---------------------------%
[p,e,t] = meshToPet(mypde.Mesh);

%size(p)
%size(e)
%size(t)
%pause(20)
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
%-------------------------- example setup ----------------------------%
x = fem.gcoord(:,1);
y = fem.gcoord(:,2);

sig0 = exp(x-x);


PATH_save = '/pvfs2/Derick/EIT/Mine/sparsity';
plot_dir  = sprintf('%s/current_%d/hval_%f-nodes_%d/plots/',PATH_save,n_current,hval,length(p(1,:)));
cond_dir  = sprintf('%s/current_%d/hval_%f-nodes_%d/result/',PATH_save,n_current,hval,length(p(1,:)));

if ~exist(plot_dir, 'dir'), mkdir(plot_dir), end
if ~exist(cond_dir, 'dir'), mkdir(cond_dir), end

for samp = 1:max_samp
    fprintf('Working on sample %d of %d \n',samp, max_samp)
    aa    = inclusions.AA(samp,:);
    alpha = inclusions.ALPHA(samp,:);
    bb    = inclusions.BB(samp,:);
    cond  = inclusions.COND(samp,:);
    hh    = inclusions.HH(samp,:);
    kk    = inclusions.KK(samp,:);

    mask = aa ~= 0;
    aa    = aa(mask);
    alpha = alpha(mask); 
    bb    = bb(mask);
    cond  = cond(mask); 
    hh    = hh(mask); 
    kk    = kk(mask); 
    
    if length(cond) == 3
        sigmatt = cond(1)*heavis(1-((((x-hh(1)).*cos(alpha(1)) + (y-kk(1)).*sin(alpha(1))).^2)./(aa(1)^2))-((((x-hh(1)).*sin(alpha(1)) - (y-kk(1)).*cos(alpha(1))).^2)./(bb(1)^2)))  +  cond(2)*heavis(1-((((x-hh(2)).*cos(alpha(2)) + (y-kk(2)).*sin(alpha(2))).^2)./(aa(2)^2))-((((x-hh(2)).*sin(alpha(2)) - (y-kk(2)).*cos(alpha(2))).^2)./(bb(2)^2)))   +   cond(3)*heavis(1-((((x-hh(3)).*cos(alpha(3)) + (y-kk(3)).*sin(alpha(3))).^2)./(aa(3)^2))-((((x-hh(3)).*sin(alpha(3)) - (y-kk(3)).*cos(alpha(3))).^2)./(bb(3)^2)));
    elseif length(cond) == 2
        sigmatt = cond(1)*heavis(1-((((x-hh(1)).*cos(alpha(1)) + (y-kk(1)).*sin(alpha(1))).^2)./(aa(1)^2))-((((x-hh(1)).*sin(alpha(1)) - (y-kk(1)).*cos(alpha(1))).^2)./(bb(1)^2)))  +  cond(2)*heavis(1-((((x-hh(2)).*cos(alpha(2)) + (y-kk(2)).*sin(alpha(2))).^2)./(aa(2)^2))-((((x-hh(2)).*sin(alpha(2)) - (y-kk(2)).*cos(alpha(2))).^2)./(bb(2)^2)));
    elseif length(cond) == 1
        sigmatt = cond(1)*heavis(1-((((x-hh(1)).*cos(alpha(1)) + (y-kk(1)).*sin(alpha(1))).^2)./(aa(1)^2))-((((x-hh(1)).*sin(alpha(1)) - (y-kk(1)).*cos(alpha(1))).^2)./(bb(1)^2)));
    end

    if length(cond) == 3
        sigmatt_fwd = cond(1)*heavis(1-((((x_fwd-hh(1)).*cos(alpha(1)) + (y_fwd-kk(1)).*sin(alpha(1))).^2)./(aa(1)^2))-((((x_fwd-hh(1)).*sin(alpha(1)) - (y_fwd-kk(1)).*cos(alpha(1))).^2)./(bb(1)^2)))  +  cond(2)*heavis(1-((((x_fwd-hh(2)).*cos(alpha(2)) + (y_fwd-kk(2)).*sin(alpha(2))).^2)./(aa(2)^2))-((((x_fwd-hh(2)).*sin(alpha(2)) - (y_fwd-kk(2)).*cos(alpha(2))).^2)./(bb(2)^2)))   +   cond(3)*heavis(1-((((x_fwd-hh(3)).*cos(alpha(3)) + (y_fwd-kk(3)).*sin(alpha(3))).^2)./(aa(3)^2))-((((x_fwd-hh(3)).*sin(alpha(3)) - (y_fwd-kk(3)).*cos(alpha(3))).^2)./(bb(3)^2)));
    elseif length(cond) == 2
        sigmatt_fwd = cond(1)*heavis(1-((((x_fwd-hh(1)).*cos(alpha(1)) + (y_fwd-kk(1)).*sin(alpha(1))).^2)./(aa(1)^2))-((((x_fwd-hh(1)).*sin(alpha(1)) - (y_fwd-kk(1)).*cos(alpha(1))).^2)./(bb(1)^2)))  +  cond(2)*heavis(1-((((x_fwd-hh(2)).*cos(alpha(2)) + (y_fwd-kk(2)).*sin(alpha(2))).^2)./(aa(2)^2))-((((x_fwd-hh(2)).*sin(alpha(2)) - (y_fwd-kk(2)).*cos(alpha(2))).^2)./(bb(2)^2)));
    elseif length(cond) == 1
        sigmatt_fwd = cond(1)*heavis(1-((((x_fwd-hh(1)).*cos(alpha(1)) + (y_fwd-kk(1)).*sin(alpha(1))).^2)./(aa(1)^2))-((((x_fwd-hh(1)).*sin(alpha(1)) - (y_fwd-kk(1)).*cos(alpha(1))).^2)./(bb(1)^2)));
    end
    sigmatt(sigmatt==0) = 1;
    sigmatt_fwd(sigmatt_fwd==0) = 1;

    ffff = figure(2);
    pixel = [925*3 598];
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, pixel/150])
    
    subplot(1,3,1)
    pdeplot(mesh_fwd,XYData=sigmatt_fwd)%, mesh='on')
    colormap jet
    title("Truth (fwd mesh)");
    xlabel("x")
    ylabel("y")
    axis off
    axis image%equal

    subplot(1,3,2)
    pdeplot(mesh,XYData=sigmatt)%, mesh='on')
    colormap jet
    title("Truth (inverse mesh)");
    xlabel("x")
    ylabel("y")
    axis off
    axis image%equal
    %drawnow
    %pause(1) 

    %------------------------ input current ------------------------%
    Neig = n_current;
    flux = aboundflux(fem,'sin',Neig);
    wt   = ones(Neig,1);

    % simulated data
        

    KKt = stiffmat(fem,sigmatt);% stiff matrix - fine mesh 
    temp = zeros(size(flux));
    for j = 1 : Neig
        [uNt,uNtt,uNnt] = neumansolve(KKt,flux(:,j),zeros(fem.sdof,1),fem);
        temp(:,j) = uNt(e(1,:));                    % potential measurements
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
        %hf = figure(5);
        subplot(1,3,3)
        fprintf('save reconstruction ...  %d \n\n',k)
        pdeplot(fem.p,fem.e,fem.t,'xydata',invpm.sig0+dsig), colormap('jet')%, drawnow
        title("Sparsity");
        axis image, axis off, colormap('jet')
    end
    toc

    figsname =  [plot_dir 'sample_' num2str(samp) '-NL' num2str(100*NL) '-alpha' num2str(k) '.png'];
    filename =  [cond_dir 'sample_' num2str(samp) '-NL' num2str(100*NL) '-alpha' num2str(k)];

    saveas(ffff, figsname) 
    out = invpm.sig0+dsig;
    save(filename, 'out') 
   
end

