clear 
close all
clc
load('CalibrationData.mat');
addpath('um')
d1 = dir('um\*.mat');
d2 = dir('um\*.jpg');
l = length(d1);
for i=1:l
    file_names(i).depth = d1(i).name;
    file_names(i).rgb = d2(i).name;
end
[obj]=peopletracking(file_names,Depth_cam,RGB_cam,R_d_to_rgb,T_d_to_rgb);
rmpath('um');
%%
close all
for i=1:length(obj)
    figure(i);
    plot(obj{i}(:,1),obj{i}(:,2));
end