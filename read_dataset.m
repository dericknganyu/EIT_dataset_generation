% #Usage#
% path = '../dataset';

% data_boundary = read('bound.mat');
% data_domain   = read('bound.mat');
% mesh          = read('mesh.mat');

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

