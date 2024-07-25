clear;close all;clc;
% ��ȡ����
data1 = importdata('71.txt');

ptCloud = pointCloud(data1(:,1:3));
nSamplePoints = 1500;
ptCloudSampled = pcdownsample(ptCloud, 'random', nSamplePoints/ptCloud.Count);
data=ptCloudSampled.Location;

% �����ֵ������
[X, Y] = meshgrid(linspace(min(data(:,1)),max(data(:,1)),100), linspace(min(data(:,2)),max(data(:,2)),100));
%  ����ÿ����ֵ��ļ�Ȩƽ��ֵ??
Vq = griddata(data(:,1),data(:,2),data(:,3),X,Y,'v4');
%��������������ֵ����?
F = griddedInterpolant(X',Y',Vq');
% �����������ϵĲ�ֵ�� 100*100������
[XI,YI] = meshgrid(linspace(min(data(:,1)),max(data(:,1)),100),linspace(min(data(:,2)),max(data(:,2)),100));
VI = F(XI,YI);
%��ͼ?
mesh(XI,YI,VI);
%grid on;
% ����ֵ�������Ϊ�ı��ļ�
dlmwrite('interp_data.txt', [XI(:) YI(:) VI(:)], 'delimiter', ' ', 'precision', 6);
new_data=[XI(:) YI(:) VI(:)];
last_data=[new_data;data];%%ԭʼ�����빹���ĵ������
%dlmwrite('last_data.txt',last_data, 'delimiter', ' ', 'precision', 6);
figure;
plot3(new_data(:,1),new_data(:,2),new_data(:,3),'g.');
hold on;
plot3(data(:,1),data(:,2),data(:,3),'r.');
grid on;axis equal;
xlabel('X');ylabel('Y');zlabel('Z');







