clear; close all; clc;

% Define the folder identifier
folderID = input('Enter the folder identifier: ', 's');
% Base folder paths
folderPath = fullfile('restore', folderID);
trunkPath = fullfile('treeqsm', folderID, 'trunk.txt');


Files = dir(fullfile(folderPath, '*_restore.txt'));
name = sort_nat({Files.name})';
% Add the trunk points in black color
trunk = importdata(trunkPath);

branchBollsCell = cell(1, numel(Files));
branchCentroids = cell(1, numel(Files));
branch = cell(1, numel(Files));
for i_ = 1:numel(Files)
    filePath = fullfile(folderPath, name{i_});
    data = load(filePath);
    positions = data(:, 1:3);
    colors = data(:, 4:6);
    cotton_boll_points = positions(colors(:, 1) > 0.5 & colors(:, 2) > 0.5 & colors(:, 3) > 0.045, :);
    branch_points = positions(colors(:, 1) < 0.5 & colors(:, 2) < 0.5 & colors(:, 3) < 0.5, :);
        % Save the branch points and cotton boll points in 'branch' variable
    branch{i_} = [branch_points; cotton_boll_points];
    %%%-----1.---------??????????-------------------------
    V=branch_points(:,1:3);
    %%%---??????????????????????????
%%????????,V-???F-??????????????????
[F,volumeConcave] = boundary(V(:,1),V(:,2),V(:,3),1);% Pierce ??????? Qhull ??????? %s ????? 0 ? 1 ??????? s ??? 0 ??????? s ??? 1 ??????????????
%%????????????????5000???         %Boundary ???? Pierce ?????????????????s = 1 ?????????????????
[N,I,B,r] = random_points_on_mesh(V, F, 2000, 'Color', 'blue', 'MaxIter', 200);
% %writematrix([N,I],'my_sampling_p.txt')
%%%??mesh????????
% trisurf(F,V(:,1),V(:,2),V(:,3))
% axis equal;xlabel('X');ylabel('Y');zlabel('Z');
% figure;
% plot3(N(:,1),N(:,2),N(:,3),'b.');
% axis equal;xlabel('X');ylabel('Y');zlabel('Z');grid on;
%%%---

%%?????
% ptCloud = pointCloud(V);
% nSamplePoints = 1500;
% ptCloudSampled = pcdownsample(ptCloud, 'random', nSamplePoints/ptCloud.Count);
% data=ptCloudSampled.Location;


pts=N;
%pts=data;
DebugShow=1;

   [coeff, ~,score]= pca(pts);
   center=mean(pts,1);%?????
   PNum=size(pts,1); % ??????
   Dir1=coeff(:,1)';%???????
   %% %????????????axis_(1,:) ????  axis_(2,:)????  axis_(3,:)????? axis_(4,:)?????
   axis_=find_AxisByPrincipalDir_mt(pts,Dir1,center,false); %???????XYZ??????AXIS???
   newpts=Transfer_XYZ2AXIS_mt(pts,axis_);                
   X1=newpts(:,3);Y1=newpts(:,1);Z1=newpts(:,2);
   [minX,minXId]=min(X1); [maxX,maxXId]=max(X1);
   
   %%%%%construct ajmatrix %??????%%%%%%%%%%%
   D=zeros(length(newpts),length(newpts));
   sumd=0;
   k=0;
   m=length(newpts);
   for i=1:length(newpts)-1
     for j=i+1:length(newpts)  
       vi=newpts(i,:);
       vj=newpts(j,:);
       d=norm(vi-vj);
       sumd=sumd+d;
       k=k+1;
       D(i,j)=d;
       D(j,i)=d;
     end
   end
   
%????????k????????????????
W1=zeros(m,m);
k=32;
for i=1:m
A=D(i,:);
t=sort(A(:));%??????????????????????
 k_=k;
if(length(t)<k)
    k_=length(t);
end
[row,col]=find(A<=t(k_),k_);%?????K???????
for j=1:k_
c=col(1,j);
 W1(i,c)=D(i,c); %W1(i,c)=1;%?k???????
end
end
for i=1:m
    for j=1:m
        if W1(i,j)==0&&i~=j
            W1(i,j)=inf;
        end
    end
end

[leafLen,LeafLenPath]=mydijkstra(W1,minXId,maxXId); %?????branch.txt????????????
%leafLen????????LeafLenPath?????
% ???W1--?????W1?i?j)??i?j????????????
% minXId--?????, maxXId--????? 
%%LenPatPts???????
   % LenPath
   pt1=pts(1,:); %?????????????
   pt2=pts(2,:); %?????????????
   %%???????????????????????????
   dis1=sqrt(pt1(2)*pt1(2)+pt1(3)*pt1(3));
   dis2=sqrt(pt2(2)*pt2(2)+pt2(3)*pt2(3));
   O_Pt=[];%%%%%%%%%???????
   if(dis1>dis2)
     O_Pt=pt2;
     LeafLenPath=fliplr(LeafLenPath);
   else
     O_Pt=pt1;
   end
   LenPathPts=pts(LeafLenPath,:);%???????????,????????????????????????  
   [~,maxId]=max(LenPathPts(:,1));
   T_Pt=pts(LeafLenPath(maxId),:); %%???????
   Dir=(T_Pt-O_Pt)./norm(T_Pt-O_Pt);
%%%??????
%      if(DebugShow)
%           figure('Name','LeafTrait' ,'NumberTitle','off');set(gcf,'color','white');movegui('southwest'); 
%      scatter3(pts(:,1),pts(:,2),pts(:,3),5,[0 0 0], 'filled');
%      hold on;
%      scatter3(pts(LeafLenPath,1),pts(LeafLenPath,2),pts(LeafLenPath,3),100,[1 0 0], 'filled');%     hold on;
%      hold on;
%      axis off; axis equal; camorbit(0,0,'camera'); axis vis3d; view(-180,0);%view3d rot;   
%      end
% leafLen
% dlmwrite('shortpath.txt',LenPathPts)

    %% ??? cotton_boll_points ????????????2???????
    minDist = 0.008;
    ptCloud_in = pointCloud(cotton_boll_points(:, 1:3));
    [labels, numClusters] = pcsegdist(ptCloud_in, minDist);  %%numClusters?????????
    %%???????
    branchBolls = cell(1, numClusters);
    for j_ = 1:numClusters
        clusterIndices = find(labels == j_);
        clusterData = cotton_boll_points(clusterIndices, :); 
        [~, volume] = boundary(clusterData(:, 1), clusterData(:, 2), clusterData(:, 3), 1);
        volume_cm3 = volume * (100^3); % ??????????
        if volume_cm3 >= 2
            branchBolls{j_} = clusterData;
            centerOfMass = mean(clusterData(:, 1:3)); % ????
            %%?????LenPathPts?????
            distances = pdist2(centerOfMass, LenPathPts); 
            [minDist, ~] = min(distances(:));
            if minDist>=0.1
            branchCentroids{i_}{j_} = minDist; % ????
            
            end
        end
    end 
    branchBolls = branchBolls(~cellfun('isempty', branchBolls));   %%??????????
    branchBollsCell{i_} = branchBolls;  %%????????????????
   

end
%%%%%%%---------------????????------------
nonEmptyCells = sum(~cellfun(@isempty, branchCentroids));
disp(['Number of vegetative branches: ', num2str(nonEmptyCells)]);
disp(['Number of fruit branches: ', num2str(length(branchCentroids)-nonEmptyCells)]);

%%%%%%%%%%%%%-------------------??????????????????????????

% % Create a figure for visualization
figure('Name', 'Point Cloud Visualization', 'NumberTitle', 'off');
hold on;
% Create scatter plot objects for each category
VB_scatter = scatter3(0, 0, 0, 5, 'r', 'filled');
FB_scatter = scatter3(0, 0, 0, 5, 'g', 'filled');
trunk_scatter = scatter3(0, 0, 0, 5, 'k', 'filled');

for i_ = 1:numel(branchCentroids)
    if ~isempty(branchCentroids{i_})
        % Non-empty centroid, plot branchBollsCell points as VB (red)
        branchBolls = branch{i_};
        if ~isempty(branchBolls)
            VB_points = branchBolls;
            scatter3(VB_points(:, 1), VB_points(:, 2), VB_points(:, 3), 5, 'r', 'filled');
            set(VB_scatter, 'XData', VB_points(:, 1), 'YData', VB_points(:, 2), 'ZData', VB_points(:, 3));
        end
    else
        % Empty centroid, plot branchBollsCell points as FB (green)
        branchBolls = branch{i_};
        if ~isempty(branchBolls)
            FB_points = branchBolls;
            scatter3(FB_points(:, 1), FB_points(:, 2), FB_points(:, 3), 5, 'g', 'filled');
            set(FB_scatter, 'XData', FB_points(:, 1), 'YData', FB_points(:, 2), 'ZData', FB_points(:, 3)); 
        end
    end
end

% Add the trunk points in black color
%trunk = importdata('treeqsm/DES56/trunk.txt');
scatter3(trunk(:, 1), trunk(:, 2), trunk(:, 3), 5, 'k', 'filled');
set(trunk_scatter, 'XData', trunk(:, 1), 'YData', trunk(:, 2), 'ZData', trunk(:, 3));
% Adjust the range of the Z-axis tick marks
%zlim([min(trunk(:, 3)) max(trunk(:, 3))]);
 
axis equal;
xlabel('X(m)');
ylabel('Y(m)');
zlabel('Z(m)');
grid on;
%legend('Vegetative Branch (VB)', 'Fruit Branch (FB)', 'Main stem');
title(['VB: ', num2str(nonEmptyCells), '   FB: ', num2str(length(branchCentroids) - nonEmptyCells)]);
 

% ?scatter3????????
% ???????????????
allPoints = [];
for i_ = 1:numel(branchCentroids)
    branchBolls = branch{i_};
    if ~isempty(branchBolls)
        allPoints = [allPoints; branchBolls];
    end
end
%trunk = importdata('treeqsm/1b/trunk.txt');
allPoints = [allPoints; trunk];
minValues = min(allPoints);
maxValues = max(allPoints);
% ???????
axis([minValues(1) maxValues(1) minValues(2) maxValues(2) minValues(3) maxValues(3)]);

set(gca, 'Color', 'w');
set(gca, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
set(gcf, 'Color', 'w');
%grid on;
%hold off;

% Hide the axis
axis off;


















