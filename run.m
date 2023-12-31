%Script to produce a set of random ellipses 
%clearvars, 
%close all
%clc
%%
% seed = 1;
% rng(seed)

% We generate and save 'batchSize' samples (with 'Nmax-Nmin+1' currents), 'numRuns' times.
% Each sample has not more than 'max_numInc' inclusions, with or without
% textures, depending on the values of textures.

function [] = run(numRuns, texture, max_numInc, batchSize, Nmax, Nmin)
    seed = 2;
    rng(seed)
    % Check if the input arguments are provided, and if not, set default values
    if nargin < 1; numRuns    = 100 ; end
    if nargin < 2; texture    = true; end
    if nargin < 3; max_numInc = 3   ; end
    if nargin < 4; batchSize  = 100 ; end
    if nargin < 5; Nmax       = 32  ; end
    if nargin < 6; Nmin       = 1   ; end
    
    if ischar(numRuns   ); numRuns    = str2double(numRuns   ); end
    if ischar(max_numInc); max_numInc = str2double(max_numInc); end  
    if ischar(batchSize ); batchSize  = str2double(batchSize ); end
    if ischar(Nmax      ); Nmax       = str2double(Nmax      ); end   
    if ischar(Nmin      ); Nmin       = str2double(Nmin      ); end

    if strcmpi(texture, 'true') || strcmpi(texture, 'false') 
        texture = strcmpi(texture, 'true');
    end
    if ~islogical(texture)
        numRuns = 1; max_numInc = 1; batchSize = 1;
    end

    for iii=1:numRuns
        %texture = true;
        %max_numInc = 3;
        %batchSize = 100;%4%32; 
        %Nmin = 1;  
        %Nmax = 32;

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
        KX    = zeros(batchSize,max_numInc);
        KY    = zeros(batchSize,max_numInc);

        %count=0;
        batchRun = 1;
        %warning('')  % Clear last warning message
        %%

        %tic
        while batchRun <=batchSize
        %for batchRun=1:batchSize

        [numInc, backCond, cond, condOut, h, k, a, b, alpha, kx, ky] = gen_conductivity(mesh, max_numInc, texture);
        % warning('off','all')   
        % fprintf('\nnumber of inclusions = %d,',numInc) 
        % disp(cond)

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
                fprintf('run %d of %d: batch %d/%d sample %d/%d\n',iii, numRuns, batchRun, batchSize, N, nn);
                [u, warnmsg]=fem_eit_fwd(numInc, backCond, cond, h, k, a, b, alpha, N, mypde, texture);%, idx_circum, theta);%, mesh);
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
            %%
            if isempty(warnmsg)
                inputConductivity(batchRun, :) = condOut;    
                HH(batchRun, 1:numInc)    = h;
                KK(batchRun, 1:numInc)    = k;
                AA(batchRun, 1:numInc)    = a;
                BB(batchRun, 1:numInc)    = b;
                ALPHA(batchRun, 1:numInc) = alpha;
                COND(batchRun, 1:numInc)  = cond;
                KX(batchRun, 1:numInc)    = kx;
                KY(batchRun, 1:numInc)    = ky;
                % fprintf('\nnumber of inclusions = %d,',numInc) 
                % disp(cond)
                % fprintf('\n saved value: ') 
                % disp(COND(batchRun, :))
                % fprintf('\n unique value: ') 
                % disp(unique(condOut))

                fprintf('\n') 
                toc

                batchRun = batchRun+1;
            end
            
        end
        if texture==true
            PATH = '/pvfs2/Derick/EIT/Mine/data_texture'; %'/localdata/Derick/EIT/Mine/data';
        elseif texture==false
            PATH = '/pvfs2/Derick/EIT/Mine/data'; %'/localdata/Derick/EIT/Mine/data';
        else
            PATH = '/pvfs2/Derick/EIT/Mine/data_constant'; %'/localdata/Derick/EIT/Mine/data';
        end
        
        if ~exist(sprintf('%s/%s',PATH, name), 'dir')
            mkdir(PATH, name)
        end
        x1 = single(x1);
        x2 = single(x2);
        radius = single(radius);
        theta = single(theta);
        inputConductivity = single(inputConductivity);
        %outputVoltage = single(outputVoltage); #Not needed
        angl_circum = single(angl_circum);
        outputBoundcurrent = single(outputBoundcurrent);
        outputBoundvoltage = single(outputBoundvoltage);
        
        HH = single(HH);
        KK = single(KK);
        AA = single(AA);
        BB = single(BB);
        ALPHA = single(ALPHA);
        COND = single(COND);
        KX   = single(KX);
        KY   = single(KY);

        save(sprintf('%s/%s/dataset_domain',PATH, name), 'x1', 'x2', 'radius', 'theta', 'inputConductivity'); %, 'outputVoltage');
        save(sprintf('%s/%s/dataset_bound',PATH, name), 'angl_circum', 'outputBoundcurrent', 'outputBoundvoltage');
        %save(sprintf('%s/%s/mesh'         ,PATH, name), 'mesh')
        save(sprintf('%s/%s/inclusions',PATH, name), 'HH', 'KK', 'AA', 'BB', 'ALPHA', 'COND', 'KX', 'KY'); %initially forgot to save KX, KY in the dataset saved in google drive
        if texture==true
            save(sprintf('%s/%s/textures',PATH, name), 'KX', 'KY');
        end

        %sprintf('data/%s/dataset_domain',name)

        %{
        angl_circum;
        difference = angl_circum(2:end) - angl_circum(1:end-1);
        uniq = unique(difference)
        %}
        %clearvars;                                                                 


    end
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