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