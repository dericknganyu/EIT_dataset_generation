%Script to produce a set of random ellipses 
clearvars, 
%close all
%clc
%%
% seed = 1;
% rng(seed)
for iii=1:100
    texture = false;
    max_numInc = 3;
    batchSize = 100;%32;%4%32; 
    Nmin = 1;  Nmax = 32;

    %%
    mypde = createpde();
    geometryFromEdges(mypde,@circleg);

    hmax = 0.03;% 0.013;%38%0.01;%4;%0.038;%0.005;
    hmin = 0.03;% 0.013;%0.01;%4;
    mesh = generateMesh(mypde,"Hmax",hmax,"Hmin",hmin);

    %size(mesh.Nodes)
    % figure(1),clf
    % pdemesh(mypde)
    % axis equal
    % axis off
    % title("Mesh");

    %save data/mesh mesh
    Nodes = mesh.Nodes;

    x1 = Nodes(1,:);
    x2 = Nodes(2,:);
    [radius, theta] = polarcoord(x1, x2);
    idx_circum = find(radius==1);

    angl_circum = theta(idx_circum);
    [angl_circum,sortIdx] = sort(angl_circum,'ascend');

    %%
    name = sprintf('%d_samples__max_Inclusions_%d__',batchSize,max_numInc) + string(datetime('now','Format',"yyyy-MM-dd-HH-mm-ss"));

    nx = numel(x1);
    nx_circum = numel(angl_circum);
    nn = Nmax-Nmin+1; %basis function number 

    inputConductivity = zeros(batchSize,nx);
    outputVoltage = zeros(batchSize,nn,nx);
    outputBoundvoltage = zeros(batchSize,nn,nx_circum);
    outputBoundcurrent = zeros(batchSize,nn,nx_circum);
    HH    = zeros(batchSize,max_numInc);
    KK    = zeros(batchSize,max_numInc);
    AA    = zeros(batchSize,max_numInc);
    BB    = zeros(batchSize,max_numInc);
    ALPHA = zeros(batchSize,max_numInc);
    COND  = zeros(batchSize,max_numInc);


    count=0;
    batchRun = 1;
    %warning('')  % Clear last warning message
    %%

    %tic
    while batchRun <=batchSize
    %for batchRun=1:batchSize

        [numInc, backCond, cond, condOut, h, k, a, b, alpha] = gen_conductivity(mesh, max_numInc, texture);

        % figure(2),clf
        % pdeplot(mesh,XYData=condOut)%, mesh='on')
        % colormap jet
        % title("Conductivity");
        % xlabel("x")
        % ylabel("y")
        % axis off
        % axis equal
        % drawnow
        %pause() 
        tic
        %%
        for N=Nmin:Nmax
            fprintf('run %d: batch %d/%d sample %d/%d\n',iii, batchRun, batchSize, N, nn);
            [u, warnmsg]=fem_eit_fwd_v2(numInc, backCond, cond, h, k, a, b, alpha, N, mypde, texture);%, idx_circum, theta);%, mesh);
            if not(isempty(warnmsg)),break,end

            u_circum = u(idx_circum);
            u_circum = u_circum(sortIdx);
            bound_mean = mean(u_circum);
            u_circum = u_circum - bound_mean;
            u = u-bound_mean;
            
            outputVoltage(batchRun, N,:) = u;
            outputBoundvoltage(batchRun, N, :) = u_circum;
            
            % figure(4)
            % pdeplot(mesh,XYData=u)
            % colormap jet
            % title("Numerical Solution");
            % xlabel("x")
            % ylabel("y")  
            % axis off
            % axis equal


            % figure(5)
            
            % % Top plot
            
            if N<=16
                i_circum = (1/sqrt(pi))*sin(N*angl_circum);
                %legend_name = {"$\frac{1}{\sqrt{\pi}} \sin($"+N+"$\cdot \theta)$"};
            else 
                i_circum = (1/sqrt(pi))*cos((N-16)*angl_circum);
                %legend_name = {"$\frac{1}{\sqrt{\pi}} \cos($"+num2str(N-16)+"$\cdot \theta)$"};
            end        
            
            outputBoundcurrent(batchRun, N, :) = i_circum;
            % subplot(3,1,1)
            % plot(angl_circum, i_circum, 'color', 'red') 
            % legend(legend_name,'Interpreter','latex') 
            
            % %title("Injected current, with N ="+N)
            % ylabel('Injected current')
            % set(gca,'XTick',0:pi/2:2*pi) 
            % set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})

            % % Bottom plot
            % subplot(3,1,2)
            % %ax2 = nexttile;
            % plot(angl_circum, u_circum)
            % title('Resulting voltage')
            % ylabel('voltage')
            % set(gca,'XTick',0:pi/2:2*pi) 
            % set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
            
            % subplot(3,1,3)
            % %ax3 = nexttile;
            % diff = transpose(i_circum)-u_circum;
            % plot(angl_circum,diff)
            
            % ylabel('voltage-current')
            % set(gca,'XTick',0:pi/2:2*pi)  
            % set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
            
            
        end
        
        toc
        %%
        inputConductivity(batchRun, :) = condOut;
    
        HH(batchRun, 1:size(h,2))    = h;
        KK(batchRun, 1:size(h,2))    = k;
        AA(batchRun, 1:size(h,2))    = a;
        BB(batchRun, 1:size(h,2))    = b;
        ALPHA(batchRun, 1:size(h,2)) = alpha;
        COND(batchRun, 1:size(h,2))  = cond;
        
        if isempty(warnmsg)
            batchRun = batchRun+1;
        end
    
    end

    PATH = '/pvfs2/Derick/EIT/Mine/data'; %'/localdata/Derick/EIT/Mine/data';
    if ~exist(sprintf('%s/%s',PATH, name), 'dir')
        mkdir(PATH, name)
    end
    x1 = single(x1);
    x2 = single(x2);
    radius = single(radius);
    theta = single(theta);
    inputConductivity = single(inputConductivity);
    outputVoltage = single(outputVoltage);
    angl_circum = single(angl_circum);
    outputBoundcurrent = single(outputBoundcurrent);
    outputBoundvoltage = single(outputBoundvoltage);
    
    HH = single(HH);
    KK = single(KK);
    AA = single(AA);
    BB = single(BB);
    ALPHA = single(ALPHA);
    COND = single(COND);

    save(sprintf('%s/%s/dataset_domain',PATH, name), 'x1', 'x2', 'radius', 'theta', 'inputConductivity'); %, 'outputVoltage');
    save(sprintf('%s/%s/dataset_bound',PATH, name), 'angl_circum', 'outputBoundcurrent', 'outputBoundvoltage');
    save(sprintf('%s/%s/mesh'         ,PATH, name), 'mesh')
    save(sprintf('%s/%s/inclusions',PATH, name), 'HH', 'KK', 'AA', 'BB', 'ALPHA', 'COND');

    %sprintf('data/%s/dataset_domain',name)

    %{
    angl_circum;
    difference = angl_circum(2:end) - angl_circum(1:end-1);
    uniq = unique(difference)
    %}
    clearvars;                                                                 


end

function [rad, theta]=polarcoord(xx, yy)
    rad = sqrt(xx.^2 + yy.^2);
    %theta = zeros(size(xx));
    theta = atan(yy./xx);
    theta(xx <= 0) = theta(xx <= 0) + pi;  
    theta(xx >0 & yy<0) = theta(xx >0 & yy<0) + 2*pi;
    
    %theta(theta>pi & theta<=1.5*pi) = theta(theta>pi & theta<=1.5*pi) - 2*pi;
    %theta = rad2deg(theta);
end                                              