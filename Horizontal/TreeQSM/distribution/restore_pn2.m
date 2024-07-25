clear; clc;

% Base folders
baseTreeqsm = 'treeqsm';
basePn2 = 'pn2';
baseRestore = 'restore';

% Get all subfolders in the baseTreeqsm folder
subfolders = dir(fullfile(baseTreeqsm, '**', ''));
subfolders = subfolders([subfolders.isdir] & ~ismember({subfolders.name}, {'.', '..'}));

% Create restore folders if they don't exist
for folder = subfolders'
    restore_folder = fullfile(baseRestore, folder.name);
    if ~exist(restore_folder, 'dir')
        mkdir(restore_folder);
    end
end

% Process each subfolder
for folder = subfolders'
    treeqsm_folder = fullfile(baseTreeqsm, folder.name);
    pn2_folder = fullfile(basePn2, folder.name);
    restore_folder = fullfile(baseRestore, folder.name);

    % Get prediction files
    pred_files = dir(fullfile(pn2_folder, '*_pred.obj'));
    num_files = length(pred_files);

    for i = 1:num_files
        % Define file paths
        original_file = fullfile(treeqsm_folder, [num2str(i-1), '.txt']);
        pred_file = fullfile(pn2_folder, [num2str(i-1), '_pred.obj']);

        % Check if original file exists
        if ~isfile(original_file)
            disp(['Original file ', original_file, ' does not exist. Skipping...']);
            continue;
        end

        % Check if prediction file exists
        if ~isfile(pred_file)
            disp(['Prediction file ', pred_file, ' does not exist. Skipping...']);
            continue;
        end

        % Load original point cloud
        pc = load(original_file);
        pc_coords = pc(:, 1:3);
        centroid = mean(pc_coords, 1);
        pc_centered = pc_coords - centroid;
        m = max(sqrt(sum(pc_centered .^ 2, 2)));

        % Load prediction data
        pred_data = load_obj_file(pred_file);
        pred_coords = pred_data(:, 1:3);
        pred_classes = pred_data(:, 4:6);

        % Check dimensions
        if size(pred_coords, 2) ~= size(centroid, 2)
            error(['Dimension mismatch between pred_coords and centroid for file ', pred_file]);
        end

        % Restore coordinates
        restored_coords = pred_coords * m + repmat(centroid, size(pred_coords, 1), 1);
        restored_data = [restored_coords, pred_classes];
        restored_data = unique(restored_data, 'rows');

        % Save restored data
        restore_file = fullfile(restore_folder, [num2str(i-1), '_restore.txt']);
        save(restore_file, 'restored_data', '-ascii');
    end
end

function points = load_obj_file(filename)
    % Load OBJ file and extract points
    fid = fopen(filename, 'r');
    points = [];
    while ~feof(fid)
        line = fgetl(fid);
        if startsWith(line, 'v ')
            coords = sscanf(line(3:end), '%f %f %f %f %f %f');
            points = [points; coords'];
        end
    end
    fclose(fid);
end
