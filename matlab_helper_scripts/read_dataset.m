% # How to read dataset in the sub-directories of the dataset directory
% # Usage #
% path = '../dataset';

% data_boundary = read_dataset('bound.mat');
% data_domain   = read_dataset('domain.mat');
% mesh          = read_dataset('mesh.mat');

function [data] = read_dataset(path, ending)
    dir_list = dir(path);
    for i = 1:numel(dir_list)
        file = dir_list(i).name;
        if ~dir_list(i).isdir && endsWith(file, ending)
            file_path = fullfile(path, file);
            % Process the file
            data = load(file_path);
        end
    end
end

