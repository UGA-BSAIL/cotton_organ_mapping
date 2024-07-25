function saveBranches(path, segment2, sop, P, parentBranchIndex, filename)
    ChiSeg = segment2.ChildSegment;
    branchIndices = ChiSeg{parentBranchIndex, 1};
    branchPoints = [];

    for i = 1:length(branchIndices)
        branchIndex = branchIndices(i);
        branchPoints = [branchPoints; find(sop == branchIndex)];

        subBranchIndices = ChiSeg{branchIndex, 1};

        if ~isempty(subBranchIndices)
            % Recursively process sub-branches
            saveBranches(path, segment2, sop, P, branchIndex, filename);
        end
    end

    % Accumulate points for the current branch
    branchPoints = [branchPoints; find(sop == parentBranchIndex)];

    % Save the combined points of the parent branch and its sub-branches
    if ~isempty(branchPoints)
        sop_pc = P(branchPoints, :); % Get point cloud data for the branches
        if length(sop_pc) >= 200 % Only save if there are 200 or more points
            fid = fopen(filename, 'a');
            lenp = size(sop_pc, 1);
            X = zeros(lenp, 1); % Create a column of zeros
            a = [sop_pc, X];
            fprintf(fid, [repmat('%10.6f\t', 1, size(a, 2)), '\n'], a'); % Write data to file
            fclose(fid);
        end
    end
end