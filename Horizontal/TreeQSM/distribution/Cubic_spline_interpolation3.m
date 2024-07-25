clear;close all;clc;
% 读取数据
data1 = importdata('71.txt');

ptCloud = pointCloud(data1(:,1:3));
nSamplePoints = 1500;
ptCloudSampled = pcdownsample(ptCloud, 'random', nSamplePoints/ptCloud.Count);
data=ptCloudSampled.Location;

% 构造插值点网格
[X, Y] = meshgrid(linspace(min(data(:,1)),max(data(:,1)),100), linspace(min(data(:,2)),max(data(:,2)),100));
%  计算每个插值点的加权平均值??
Vq = griddata(data(:,1),data(:,2),data(:,3),X,Y,'v4');
%构造三次样条插值函数?
F = griddedInterpolant(X',Y',Vq');
% 计算新网格上的插值点 100*100的网格
[XI,YI] = meshgrid(linspace(min(data(:,1)),max(data(:,1)),100),linspace(min(data(:,2)),max(data(:,2)),100));
VI = F(XI,YI);
%绘图?
mesh(XI,YI,VI);
%grid on;
% 将插值结果保存为文本文件
dlmwrite('interp_data.txt', [XI(:) YI(:) VI(:)], 'delimiter', ' ', 'precision', 6);
new_data=[XI(:) YI(:) VI(:)];
last_data=[new_data;data];%%原始点云与构建的点合起来
%dlmwrite('last_data.txt',last_data, 'delimiter', ' ', 'precision', 6);
figure;
plot3(new_data(:,1),new_data(:,2),new_data(:,3),'g.');
hold on;
plot3(data(:,1),data(:,2),data(:,3),'r.');
grid on;axis equal;
xlabel('X');ylabel('Y');zlabel('Z');







