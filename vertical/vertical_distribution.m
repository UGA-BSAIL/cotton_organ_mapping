% Cotton boll count and vertical spatial distribution analysis on the whole plant.
% The first part is to implement point cloud Euclidean distance clustering, complete the count, and display it on the image
close all;clear;clc;

%% 1. load point cloud and Euclidean cluster
PlantNumber = input('Enter the plant number (from 0 to n): ', 's');
IdName = [PlantNumber '_pred'];
filename = ['output_plant\' IdName '.obj']; %%prediction file stored in this folder
disp(['Generated filename: ' filename]);
ptCloudData= importdata(filename);
ptCloud=ptCloudData.data;
%writematrix(ptCloud, 'plant_points.txt');

% Define category information
branchClass = [0.35, 0.05, 0.35];
bollClass = [0.65, 0.95, 0.05];
% Set tolerance for comparison
tolerance = 1e-6;
% Find points belonging to the branch category
branchIndices = all(abs(ptCloud(:, 4:6) - branchClass) < tolerance, 2);
branch = ptCloud(branchIndices, 1:3);
% Find points belonging to the boll category
bollIndices = all(abs(ptCloud(:, 4:6) - bollClass) < tolerance, 2);
boll = ptCloud(bollIndices, 1:3);
% Save the results (if needed)
%writematrix(branch, 'branch_points.txt');
%writematrix(boll, 'boll_points.txt');

boll_in= pointCloud(boll(:,1:3));
%pcwrite(ptCloud, 'test.pcd', 'Encoding', 'ascii'); %Write the txt data the pcd file
%pcshow(boll_in); %show point cloud
% Euclidean distance
minDist =0.008;       % 0.8cm in our boll count 
% Euclidean cluster
[labels,numClusters] = pcsegdist(boll_in,minDist);

% Count number of points in each cluster
clusterCounts = histcounts(labels, 1:(numClusters + 1));
% Plot line graph of cluster sizes
%figure;
plot(1:numClusters, clusterCounts, '-o');
title('Number of points in each cluster','Color', 'k');
xlabel('Cluster Index');
ylabel('Number of Points');
grid on;

% show raw cluster results
% figure; 
% pcshow(boll_in.Location,labels);
 %colormap(hsv(numClusters));
% title('Euclidean clustering');
 %xlabel('X(m)');
% ylabel('Y(m)');
 %zlabel('Z(m)');

%%2. Remove clusters with less than 100 points for accurate counting
table=tabulate(labels); %%The first column is the clustering, the second column is the number of points in each cluster, 
                        %%and the third column is the proportion.
Remo_labels=find(table(:,2)<=100);%Remove
%labels(labels>=min(Remo_labels)&labels<=max(Remo_labels))=0;%去除的位置赋值为0。在这个区间的数值赋为0.
%labels(sub2ind(size(labels), Remo_labels))=0;
%[c, ia, ib] = intersect(Remo_labels,labels,'rows');
%%Generate a cell with the size of the category to be merged and save it in a cell.
len_Re=length(Remo_labels);
cell_Re= cell(len_Re,1);%Store the points that need to be deleted
for i=1:len_Re
    cell_Re{i}=find(labels(:,1)==Remo_labels(i));
end
%%Parse all indexes in the cell
T = vertcat(cell_Re{:,1});%%Store all index values ??in the label and assign all values ??of these indexes to 0. 
                          %%At this time, the cell is saved in the form of columns. If it is a row -- T = vertcat(cell{1,:});
labels(sub2ind(size(labels), T))=0;%The deleted categories are replaced by 0
labels(all(labels==0,2),:) = []; %Delete rows with all zeros -- in the labels index
pc_r=boll_in.Location;
pc_r(T,:) = []; %Delete rows with all zeros - in the original point cloud, T is the index of the point to be deleted
%%show new cluster
figure;
pcshow(pc_r,labels'); % Optimized clustered cotton bolls
colormap(hsv(numClusters));
%title('');
xlabel('x(m)');
ylabel('y(m)');
zlabel('z(m)');
%%Number of clusters after merging
new_numClusters=numClusters-length(Remo_labels);
%%load branch
branch=pointCloud(branch(:,1:3));
hold on;
pcshow(branch.Location,"k");
title(['Predicted bolls: ',num2str(new_numClusters)],'Color', 'k') %show the boll number in the title
%title(['Predicted bolls: ',num2str(new_numClusters),', Vs ground truth bolls: 47'],'Color', 'k')
%%%change background color
set(gcf, 'Color', 'w');
set(gca, 'Color', 'w');
set(gca, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
%%The above shows the results of Euclidean clustering, achieving full plant counting and display, 
%%showing the first image, and refined clustering
%%%--------------------------------------------------------------



%%In the second part, the cotton bolls are divided into three parts: upper, middle and lower according to plant height, and displayed in different colors
%%Get the centroid of each cotton boll--the centroid is layered by plant height
cell_R= cell(new_numClusters,1); %%Store the points of each retained cluster(the last cotton boll)
cell_cen= cell(new_numClusters,1);%%Store the centroid of each cluster, and the point cloud corresponds to the centroid one by one.
                                  %%That is, the centroid of cell_cen{1,1} corresponds to all points in the interval cell_R{1,1}
labels_R=unique(labels);
for i=1:new_numClusters
    cell_R_index=find(labels(:,1)==labels_R(i));
    cell_R_pc=pc_r(cell_R_index,:);
    cell_R{i}=cell_R_pc;
    [cell_R_idx,cell_R_center] = kmeans(cell_R_pc,1);
    cell_cen{i}=cell_R_center;
end
%%How to display these centroid points--according to the upper, middle and lower parts, 
%%and the corresponding cotton bolls, and all the points after adding the trunk
%Display the centroid of the cluster cell_cen
center_=cell2mat(cell_cen); %centroid
%figure;
% pcshow(center_,labels_R,'MarkerSize',160);%MarkerSize--the size of each point
% colormap(hsv(new_numClusters));
%Category of centroids based on plant height
%Plant height
m_b_pc=branch.Location;
m_b=[min(m_b_pc(:,3)), max(m_b_pc(:,3))]; %branch
m_t_pc=boll_in.Location;
m_t=[min(m_t_pc(:,3)), max(m_t_pc(:,3))]; %cotton boll
m_h=[min(min(m_b_pc(:,3)),min(m_t_pc(:,3))),max(max(m_b_pc(:,3)),max(m_t_pc(:,3)))]; %Height range of the entire plant
height_plant=m_h(2)-m_h(1);
%Evenly divided into three parts: upper, middle and lower m_h(1)+height_plant/3;m_h(1)+height_plant*2/3; Z axis
part_1=m_h(1)+height_plant/3;
part_2=m_h(1)+height_plant*2/3;
bottom=find(center_(:,3)>=m_h(1) & center_(:,3)<part_1);
bottom_pc=center_(bottom,:);
middle=find(center_(:,3)>=part_1 & center_(:,3)<part_2);
middle_pc=center_(middle,:);
top=find(center_(:,3)>=part_2 & center_(:,3)<m_h(2));
top_pc=center_(top,:);
grid on;
%%Only show the vertical distribution of the centroid points
%figure;hold on; 
%pcshow(bottom_pc,'^r','MarkerSize',30);
%pcshow(middle_pc,'og','MarkerSize',30);
%pcshow(top_pc,'pentagramb','MarkerSize',40);
%xlabel('x');
%ylabel('y');
%zlabel('z');
%legend('Bottom','Middle','Top','FontSize',8)
%title(['Top: ',num2str(length(top)),', ','Middle: ',num2str(length(middle)),', ','Bottom: ',num2str(length(bottom))]);
%set(legend, 'TextColor', 'k');
%set(legend, 'Color', 'w'); 
%grid on;
%%%change background color
%set(gcf, 'Color', 'w');
%set(gca, 'Color', 'w');
%set(gca, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
%set(legend, 'FontName', 'Arial');
%set(legend, 'FontSize', 8);

%%Show the vertical distribution of the original point cloud(cotton boll)corresponding to the centroid
figure;hold on; 
if ~isempty(bottom) %If there is no cotton boll at the bottom, this part will not be displayed
bottom_pc_o=cell_R(bottom,:);
T_bottom=vertcat(bottom_pc_o{:,1});
pcshow(T_bottom,'r');
else 
    T_bottom=0;
end
if ~isempty(middle) 
middle_pc_o=cell_R(middle,:);
T_middle=vertcat(middle_pc_o{:,1});
pcshow(T_middle,'g');
else 
    T_middle=0;
end
if ~isempty(top) 
top_pc_o=cell_R(top,:);
T_top=vertcat(top_pc_o{:,1});
pcshow(T_top,'b');
else
    T_top=0;
end
pcshow(branch.Location,"k");
xlabel('x(m)');
ylabel('y(m)');
zlabel('z(m)');
%legend('Bottom','Middle','Top','FontSize',8);%Sometimes there is no boll below the boll, only branches
title(['Top: ',num2str(length(top)),', ','Middle: ',num2str(length(middle)),', ','Bottom: ',num2str(length(bottom))],'Color', 'k');
grid on;
%%%change background color
set(gcf, 'Color', 'w');
set(gca, 'Color', 'w');
%set(legend, 'TextColor', 'k');
%set(legend, 'Color', 'w'); 
%set(legend, 'FontName', 'Arial');
%set(legend, 'FontSize', 8);
set(gca, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
%%%The above shows the cotton boll divided into three parts: upper, 
%middle and lower (plant height is divided into three equal parts), and displayed in different colors. 
%The second and third pictures show the centroid in three different colors and the cotton boll in three different colors.
%%--------------------------------------------------------------

%%The third part draws the vertical distribution curve of cotton bolls
%center_, the center of mass coordinates of each cotton boll
figure;hold on;
zz=m_h(1):0.1:m_h(2); % Count the number of cotton bolls every 5cm
ZZ=cell(length(zz),1); %Save the number of cotton bolls in each height range and plot it
for i= 1:length(zz)
    z_s=zz(1,i);
    z_b=z_s+0.05;
    ind_z=find(center_(:,3)>=z_s & center_(:,3)<=z_b);
    count_i_z=length(ind_z);
    ZZ(i,1)={count_i_z}; %Store the number of points in each part into the corresponding position of the cell in turn    
end   
%barh(cell2mat(ZZ));
%xlabel('Number of bolls');
%ylabel('Plant height direction (5cm/interval)');
%pp_=1:1:length(ZZ);
%plot(cell2mat(ZZ),pp_,'g--','LineWidth',1.5);
%view(90,-90);
ZZ_mat = cell2mat(ZZ);
barh(ZZ_mat);

% Calculate and display proportions
total_bolls = sum(ZZ_mat);
for i = 1:length(ZZ_mat)
    percentage = (ZZ_mat(i) / total_bolls) * 100;
    text(ZZ_mat(i) + 0.1, i, [num2str(percentage, '%.f'), '%'], 'VerticalAlignment', 'middle', 'FontName', 'Arial', 'FontSize', 10);
end

xlabel('Number of bolls');
ylabel('Plant height direction (10cm/interval)');
title('Vertical distribution of cotton bolls');
grid on;
















