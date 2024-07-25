close all; clear; clc;

% Define custom color map
colorMap = [
    1 0 0;         % Red
    0 1 0;         % Green
    0 0 1;         % Blue
    1 1 0;         % Yellow
    1 0 1;         % Magenta
    0 1 1;         % Cyan
    1 0.5 0;       % Orange
    0.5 0 0.5;     % Purple
    0.5 0.5 0;     % Olive
    0.5 0 0;       % Dark Red
    0 0.5 0;       % Dark Green
    0 0 0.5;       % Dark Blue
    1 0.75 0.8;    % Light Pink
    0.5 0.5 0.5;   % Gray
    0.75 0.75 0;   % Light Yellow
    0 0.75 0.5;    % Teal
    1 0.25 0.5;    % Coral
    0.25 0.75 0.25; % Light Green
    0.75 0.25 0;   % Brown
    0.25 0.25 0.75; % Slate Blue
    0.75 0.75 0.75; % Light Gray
    0.25 0.75 0.75; % Light Cyan
    0.75 0.5 0.25; % Gold
    0.5 0 0.25;    % Dark Magenta
    0.75 0.25 0.5; % Salmon
    0.25 0.5 0.75; % Light Blue
    0.5 0.75 0.5;  % Mint
];

% Prompt for folder identifier
folderID = input('Enter the folder identifier: ', 's');

% Define file paths
trunkFile = fullfile('treeqsm', folderID, 'trunk.txt');
branchFolder = fullfile('restore', folderID);

% Load trunk data
trunkData = importdata(trunkFile);
trunkPointCloud = trunkData(:, 1:3);

% Load branches data
branchFiles = dir(fullfile(branchFolder, '*_restore.txt'));

figure;
hold on;
pcshow(trunkPointCloud, 'k', 'MarkerSize', 10); % Display trunk in black

for i = 1:numel(branchFiles)
    branchFile = fullfile(branchFolder, branchFiles(i).name);
    branchData = importdata(branchFile);
    branchPointCloud = branchData(:, 1:3);
    colorIndex = mod(i-1, size(colorMap, 1)) + 1;
    pcshow(branchPointCloud, colorMap(colorIndex, :), 'MarkerSize', 10); % Display each branch in a different color
end

%xlabel('X(m)');
%ylabel('Y(m)');
%zlabel('Z(m)');
%title('Branches and Trunk with Custom Colors');
set(gca, 'Color', 'w');
set(gca, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
set(gcf, 'Color', 'w');
%grid on;
%hold off;

% Hide the axis
axis off;
