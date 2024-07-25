clear; clc;

% Main treeqsm folder
treeqsm_folder = 'treeqsm';

% Get all branch files from all subfolders
original_files = dir(fullfile(treeqsm_folder, '**', 'branch_*.txt'));

% Extract file numbers
file_numbers = arrayfun(@(f) sscanf(f.name, 'branch_%d.txt'), original_files);
[~, sorted_idx] = sort(file_numbers);

% Rename files
for i = 1:length(sorted_idx)
    old_name = fullfile(original_files(sorted_idx(i)).folder, original_files(sorted_idx(i)).name);
    new_name = fullfile(original_files(sorted_idx(i)).folder, [num2str(i-1), '.txt']);
    movefile(old_name, new_name);
end

disp('Files have been renamed successfully.');
