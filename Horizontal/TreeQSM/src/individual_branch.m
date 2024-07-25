% Main script
path = 'segment_data/1b/';
if ~exist(path, 'dir')
    mkdir(path);
end

% Save trunk data to trunk.txt
filename = 'trunk.txt';
dlmwrite(fullfile(path, filename), trunk, 'precision', '%10.6f');

% Process and save branch segments
ChiSeg = segment2.ChildSegment;
sop = segment2.SegmentOfPoint;
Br_i = ChiSeg{1, 1}; % Get indexes of first level branches

for i = 1:length(Br_i)
    branchIndex = Br_i(i);
    branchFilename = fullfile(path, ['branch_', num2str(branchIndex), '.txt']);
    saveBranches(path, segment2, sop, P, branchIndex, branchFilename);
end