clear;close all;clc;
%%--1.-------------????????????????????????branchBollsCell?--------------------%%
%folderId = '1f';
folderId = input('Enter folder identifier: ', 's'); %input plant name
folderPath = sprintf('./restore/%s', folderId); %%read files
disFolderPath = sprintf('./dis/%s', folderId); %%save boll level
trunk = importdata(sprintf('treeqsm/%s/trunk.txt', folderId)); %%reand trunk

% ???????????
Files = dir(fullfile(folderPath, 'boll*'));
len_txt=length(Files);
name=sort_nat({Files.name})'; %%????


branchBollsCell = cell(1, numel(Files));
% ????????
for i = 1:numel(Files)
    filePath = fullfile(folderPath, Files(i).name);
    data = load(filePath);
    
    minDist = 0.008;
    ptCloud_in = pointCloud(data(:, 1:3));
    [labels, numClusters] = pcsegdist(ptCloud_in, minDist);
    
    branchBolls = cell(1, numClusters);
    
    for j = 1:numClusters
        clusterIndices = find(labels == j);
        clusterData = data(clusterIndices, :);
        branchBolls{j} = clusterData;
    end
    
    branchBollsCell{i} = branchBolls;
end

%%----2.-----------??branchBollsCell??????????100????--------------------%%
newBranchBollsCell = {};

for i = 1:numel(branchBollsCell)
    branchBolls = branchBollsCell{i};
    newBranchBolls = {};
    
    for j = 1:numel(branchBolls)
        bollData = branchBolls{j};
        
        % Check if the cluster has 100 or more points
        if size(bollData, 1) >= 100
            newBranchBolls{end+1} = bollData;
        end
    end
    
    if ~isempty(newBranchBolls)
        newBranchBollsCell{end+1} = newBranchBolls;
    end
end

% ????? branchBollsCell ???????
branchBollsCell = newBranchBollsCell;



%%----3.-----------??????branchBollsCell?????????????????????????--------------------%%
% ?????????
branchCentroids = cell(size(branchBollsCell));
for i = 1:numel(branchBollsCell)
    branchBolls = branchBollsCell{i};
    centroids = zeros(numel(branchBolls), 3);
    for j = 1:numel(branchBolls)
        bollData = branchBolls{j};
        centroids(j, :) = mean(bollData(:, 1:3));
    end
    branchCentroids{i} = centroids;
end
% ???????????????
% ????????????
distancesToTrunk = cell(size(branchCentroids));
for i = 1:numel(branchCentroids)
    centroids = branchCentroids{i};
    distances = pdist2(centroids, trunk);
    [minDistance, sortedIndices1] = min(distances, [], 2);
    distancesToTrunk{i} = minDistance;
end

% ????????????---????????sortedBranchBollsCell
sortedBranchBollsCell = cell(size(branchBollsCell));
for i = 1:numel(branchBollsCell)
    branchBolls = branchBollsCell{i};
    distances = distancesToTrunk{i};
    [~, sortedIndices] = sort(distances);
    sortedBranchBolls = cell(size(branchBolls));
    for j = 1:numel(branchBolls)
        sortedBranchBolls{j} = branchBolls{sortedIndices(j)};
    end
    sortedBranchBollsCell{i} = sortedBranchBolls;
end

%%-----4.----------????????--------------------%%
% % ??????? % ?? ?? ?? ?? ?? ??
colorMap = [
    1 0 0;  
    0 1 0;
    0 0 1; 
    1 1 0; 
    1 0 1; 
    0 1 1; 
    0.0500 0.0500 0.0500;
    0.0500 0.3500 0.0500;
    0.0500 0.3500 0.9500;
    0.0500 0.6500 0.3500;
    0.0500 0.6500 0.9500;
    0.0500 0.9500 0.0500;
    0.0500    0.9500    0.3500;
    0.0500    0.9500    0.6500;
    0.3500    0.0500    0.0500;
    0.3500    0.0500    0.3500;
    0.3500    0.0500    0.6500;
    0.3500    0.0500    0.9500;
    0.3500    0.3500    0.6500;
    0.3500    0.6500    0.9500;
    0.3500    0.9500    0.3500;
    0.3500    0.9500    0.6500;
    0.3500    0.9500    0.9500;
    0.6500    0.3500    0.3500;
    0.6500    0.6500    0.9500;
    0.6500    0.9500    0.0500;
    0.6500    0.9500    0.3500;
    0.6500    0.9500    0.9500;
    0.9500    0.0500    0.0500;
    0.9500    0.3500    0.0500;
    0.9500    0.6500    0.0500;
    0.9500    0.6500    0.6500;
    0.9500    0.9500    0.6500;
    0.9500    0.9500    0.9500;];

% ??????
%disFolderPath = './dis/2b';
if ~exist(disFolderPath, 'dir')
    mkdir(disFolderPath);
end

% ??????????????
allBollsMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

% ???????????????
for i = 1:numel(sortedBranchBollsCell)
    branchBolls = sortedBranchBollsCell{i};
    
    % ????????
    for j = 1:numel(branchBolls)
        clusterData = branchBolls{j};
        colorKey = sprintf('row_%d', j);
        
        if isKey(allBollsMap, colorKey)
            allBollsMap(colorKey) = [allBollsMap(colorKey); clusterData]; 
        else
            allBollsMap(colorKey) = clusterData;
        end
    end
end

% ???????????????
keys = allBollsMap.keys;
for k = 1:numel(keys)
    key = keys{k};
    dlmwrite(fullfile(disFolderPath, [key, '.txt']), allBollsMap(key), 'delimiter', '\t');
end



% ??????Figure??
figure;
% ???????????????
for i = 1:numel(sortedBranchBollsCell)
    branchBolls = sortedBranchBollsCell{i};    
    % ??????????????
    for j = 1:numel(branchBolls)
        clusterData = branchBolls{j};
        scatter3(clusterData(:, 1), clusterData(:, 2), clusterData(:, 3), 10, colorMap(j, :), 'filled');
        hold on;
    end
end
%%%%%%%%%%%%%%%%%%
numBolls = 100;  % Specify the number of bolls to extract
bollCell = cell(numBolls, numel(sortedBranchBollsCell));
for i = 1:numel(sortedBranchBollsCell)
    branchBolls = sortedBranchBollsCell{i};
    
    for j = 1:numBolls
        if numel(branchBolls) >= j
            boll = branchBolls{j};
            bollCell{j, i} = boll;
        end
    end
end
nonEmptyCounts = zeros(size(bollCell, 1), 1);
% Loop over each row of bollCell
for row = 1:size(bollCell, 1)
    % Get the sub-cells in the current row
    rowBolls = bollCell(row, :);   
    % Remove empty cells from the rowBolls array
    rowBolls = rowBolls(~cellfun(@isempty, rowBolls));    
    % Count the number of non-empty sub-cells
    nonEmptyCounts(row) = numel(rowBolls);    
    % Save the non-empty sub-cells in a text file
    if ~isempty(rowBolls)
        %fileName = sprintf('./dis/2b/row%d.txt', row);
        fileName = sprintf('%s/row%d.txt', disFolderPath, row);
        cellData = cell2mat(rowBolls(:));
        %dlmwrite(fileName, cellData, 'precision', '%10.6f');
    end
end
nonZeroCounts = nonEmptyCounts(nonEmptyCounts ~= 0); %%????????

%%%%????
for i = 1:numel(sortedBranchBollsCell)
    branchBolls = sortedBranchBollsCell{i};    
    % Plot each cluster of cotton bolls
    for j = 1:numel(branchBolls)
        clusterData = branchBolls{j};
        scatter3(clusterData(:, 1), clusterData(:, 2), clusterData(:, 3), 10, colorMap(j, :), 'filled');
        hold on;
    end
end
% Create the legend labels
legendLabels = cell(size(nonZeroCounts));
for i = 1:numel(nonZeroCounts)
    legendLabels{i} = sprintf('boll%d: %d', i, nonZeroCounts(i));
end

% ????????????????????????????
scatterHandles = gobjects(numel(branchBollsCell), 1);
for i = 1:numel(sortedBranchBollsCell)
    branchBolls = sortedBranchBollsCell{i};
    
    for j = 1:numel(branchBolls)
        clusterData = branchBolls{j};
        scatter3(clusterData(:, 1), clusterData(:, 2), clusterData(:, 3), 10, colorMap(j, :), 'filled');
        hold on;
    end
    
    % ????????????????
    scatterHandles(i) = scatter3(NaN, NaN, NaN, 10, colorMap(i, :), 'filled');
end

% ????????????????????
legend(scatterHandles, legendLabels);

%%-----5.----------???????????????????--------------------%%
hold on;
% Specify the folder path
%folderPath = './restore/2b';
% List all files starting with "branch"
Files = dir(fullfile(folderPath, 'branch*.txt'));
% Initialize variables
combinedData = [];
branchPoints = {};
% Loop through each file
for i = 1:numel(Files)
    % Get the file name
    fileName = Files(i).name;    
    % Build the file's complete path
    filePath = fullfile(folderPath, fileName);
    % Read the file's content
    data = importdata(filePath);    
    % Store the data in the cell array
    branchPoints{i} = data;
    % Concatenate the data to create a combined point cloud
    combinedData = [combinedData; data];
end

% Create a point cloud from the combined data
pointCloudData = pointCloud(combinedData);
% Visualize the point cloud
xyz=pointCloudData.Location;
br_tr=[xyz;trunk];
plot3(br_tr(:,1),br_tr(:,2),br_tr(:,3),'k.');
xlabel('X(m)');
ylabel('Y(m)');
zlabel('Z(m)');
set(gca, 'Color', 'w');
set(gca, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
set(gcf, 'Color', 'w');
grid on;
axis equal;




 



