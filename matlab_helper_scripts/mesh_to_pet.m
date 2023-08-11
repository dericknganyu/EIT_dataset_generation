%Script to read dataset 

%load dataset
%dataset_bound__
name = '100_samples__max_Inclusions_3__2023-08-05-03-34-06';

PATH = '/pvfs2/Derick/EIT/Mine/';

% if ~exist(sprintf('%sdata/%s/files',PATH, name), 'dir')
%     mkdir(sprintf('%sdata/%s/',PATH, name), 'files')
% else 
%     fprintf('already exists\n');
%     return
% end

%load data\mesh.mat mesh;
load(sprintf('%sdata/%s/mesh.mat',PATH, name));

[p,e,t] = meshToPet(mesh);
save(sprintf('%sdata/%s/mesh_pet',PATH, name), 'p', 'e', 't');


%--------------------------------------------------------
mypde = createpde();
geometryFromEdges(mypde,@circleg);

hmax = 0.03;% 0.013;%38%0.01;%4;%0.038;%0.005;
hmin = 0.03;% 0.013;%0.01;%4;
mesh = generateMesh(mypde,"Hmax",hmax,"Hmin",hmin);

vtkFileName = sprintf('%sdata/%s/mesh.vtk',PATH, name);
%exportGeometry(mesh, vtkFileName);

pdetool('export', vtkFileName, 'vtk', 'data', mesh);