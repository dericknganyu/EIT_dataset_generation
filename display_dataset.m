%Script to read dataset 

%load dataset
%dataset_bound__
name = '100_samples__max_Inclusions_3__2023-07-29-14-39-04';

PATH = '/pvfs2/Derick/EIT/Mine/';

if ~exist(sprintf('%sdata/%s/files',PATH, name), 'dir')
    mkdir(sprintf('%sdata/%s/',PATH, name), 'files')
else 
    fprintf('already exists\n');
    return
end

%load data\mesh.mat mesh;
load(sprintf('%sdata/%s/mesh.mat',PATH, name));

[p,e,t] = meshToPet(mesh);
save(sprintf('%sdata/%s/mesh_pet',PATH, name), 'p', 'e', 't');

data_domain = load(sprintf('%sdata/%s/dataset_domain',PATH, name));%,'inputConductivity');
inputConductivity  = data_domain.inputConductivity;
%data_domain = load(sprintf('%sdata/%s/dataset_domain',PATH, name),'outputVoltage')
outputVoltage      = data_domain.outputVoltage;

data_bound  = load(sprintf('%sdata/%s/dataset_bound',PATH, name));
outputBoundvoltage = data_bound.outputBoundvoltage;
outputBoundcurrent = data_bound.outputBoundcurrent;
angl_circum        = data_bound.angl_circum;

[batchSize, Nmax, n_nodes] = size(outputVoltage);


%get(0,'ScreenSize')
for batchRun=1:batchSize
    
    condOut = inputConductivity(batchRun,:);
    set(0,'DefaultFigureVisible','off'); %prevent matlab from showing figure
    f= figure('Name',num2str(batchRun),'Position',[-300 1120  2100 7200]);%figure(batchRun);

    subplot(Nmax+1,3,2);
    pdeplot(mesh,'XYData',condOut);
    colormap jet
    title("Conductivity",'FontWeight','Normal');
    xlabel("x")
    ylabel("y")
    axis off
    %drawnow

    
    for N=1:Nmax
        fprintf('batch %d/%d sample %d/%d\n',batchRun, batchSize, N, Nmax)
        u = outputVoltage(batchRun, N, :);
        u = squeeze(u);
        
        u_circum = outputBoundvoltage(batchRun, N, :);
        u_circum = squeeze(u_circum);
        
        i_circum = outputBoundcurrent(batchRun, N, :);
        i_circum = squeeze(i_circum);
        
        subplot(Nmax+1,3,3*N+1);
        plot(angl_circum, i_circum);
        if N<=16
            title_name = {"Injected (boundary) current: $\frac{1}{\sqrt{\pi}} \sin($"+N+"$\cdot \theta)$"};
        else 
            title_name = {"Injected (boundary) current: $\frac{1}{\sqrt{\pi}} \cos($"+num2str(N-16)+"$\cdot \theta)$"};
        end        
        title(title_name,'Interpreter','latex')
        ylabel('Boundary current')
        set(gca,'XTick',0:pi/2:2*pi) 
        set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
        
        
 
        subplot(Nmax+1,3,3*N+2);
        pdeplot(mesh,'XYData',u);
        colormap jet
        title("Numerical Solution (Domain voltage)",'FontWeight','Normal');
        xlabel("x")
        ylabel("y")  
        axis off
        
        subplot(Nmax+1,3,3*N+3);
        plot(angl_circum, u_circum);
        if N<=16
            title_name = {"Resulting boundary voltage for Injected current: $\frac{1}{\sqrt{\pi}} \sin($"+N+"$\cdot \theta)$"};
        else 
            title_name = {"Resulting boundary voltage for Injected current: $\frac{1}{\sqrt{\pi}} \cos($"+num2str(N-16)+"$\cdot \theta)$"};
        end        
        title(title_name,'Interpreter','latex')
        ylabel('Boundary voltage')
        set(gca,'XTick',0:pi/2:2*pi) 
        set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
        
    end
    
    set(f,'Units','inches');
    screenposition = get(f,'Position');
    set(f,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperOrientation','portrait',...
        'PaperType','a1',...
        'PaperSize',[screenposition(3:4)]);

    fprintf('Now saving to >> data/%s/files/batch%d.png\n',name,batchRun)
    tic
    filename = sprintf('%sdata/%s/files/batch%d.png',PATH,name,batchRun);
    saveas(f,filename)

    %print(f,'-dpdf',filename, '-opengl')
    toc

end



clearvars; 