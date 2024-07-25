clear; clc;

% Path to the restore folder
basePath = 'restore';

% Get a list of all subfolders in the restore folder
subfolders = dir(fullfile(basePath, '**', '')); % Include all subfolders

% Iterate through each subfolder
for folderIdx = 1:length(subfolders)
    if subfolders(folderIdx).isdir && ~ismember(subfolders(folderIdx).name, {'.', '..'})
        % Get the path to the current subfolder
        subfolderPath = fullfile(subfolders(folderIdx).folder, subfolders(folderIdx).name);
        
        % Get all .txt files in the current subfolder
        files = dir(fullfile(subfolderPath, '*.txt'));

        % Process each file in the current subfolder
        for fileIdx = 1:length(files)
            % Load the data from the file
            data = load(fullfile(subfolderPath, files(fileIdx).name));
            
            % Check if the data has enough columns
            if size(data, 2) < 6
                warning('File %s does not have enough columns, skipping.', files(fileIdx).name);
                continue;
            end
            
            % Separate coordinates and classes
            coords = data(:, 1:3);
            classes = data(:, 4:6);
            
            % Identify boll and branch indices
            bollIdx = all(classes == [0.65, 0.95, 0.05], 2); % cotton bolls
            branchIdx = all(classes == [0.35, 0.05, 0.35], 2); % branches
            
            % Check if all points belong to one class
            if all(bollIdx)
                % All points are cotton bolls
                bollFileName = fullfile(subfolderPath, sprintf('boll_%d.txt', fileIdx));
                dlmwrite(bollFileName, coords, 'delimiter', '\t', 'precision', '%.6f');
            elseif all(branchIdx)
                % All points are branches
                branchFileName = fullfile(subfolderPath, sprintf('branch_%d.txt', fileIdx));
                dlmwrite(branchFileName, coords, 'delimiter', '\t', 'precision', '%.6f');
            else
                % Separate and save boll data if it exists
                bollData = coords(bollIdx, :);
                if ~isempty(bollData)
                    bollFileName = fullfile(subfolderPath, sprintf('boll_%d.txt', fileIdx));
                    dlmwrite(bollFileName, bollData, 'delimiter', '\t', 'precision', '%.6f');
                end
                
                % Separate and save branch data if it exists
                branchData = coords(branchIdx, :);
                if ~isempty(branchData)
                    branchFileName = fullfile(subfolderPath, sprintf('branch_%d.txt', fileIdx));
                    dlmwrite(branchFileName, branchData, 'delimiter', '\t', 'precision', '%.6f');
                end
            end
        end
    end
end

disp('All files have been processed and restored successfully.');
