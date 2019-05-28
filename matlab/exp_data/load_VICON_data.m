clear all;
close all;

XLS=readtable('./Trial_196_Vicon_marker_position_data.xlsx');
XLS=table2array(XLS);

% read data / convert to meter
frame = str2double(XLS(5:end,1));
t = (frame-frame(1))/200;
varnames={'LW1', 'LW2', 'LW3', 'RW1', 'RW2', 'RW3', 'T1', 'T2', 'A1', 'A2', 'LHW1', 'RHW1'};
for i=1:length(varnames)
    eval([varnames{i} ' = str2double(XLS(5:end,'  num2str(3*i) ':' num2str(3*i+2) '))''/1e3;']);
end

figure;
plot3(T1(1,:),T1(2,:),T1(3,:));
hold on;
plot3(T1(1,[1, end]),T1(2,[1, end]),T1(3,[1, end]));
title('raw data');
axis equal;
grid on;

% translate the inertial frame such that T1 starts from the origin
T1_init = T1(:,1);
for i=1:length(varnames)
    eval([varnames{i} ' = ' varnames{i} '- T1_init*ones(1,length(T1));']);
end

% rotate the inertial frame such that the line connecting the initial T1
% and the terminal T1 is along the x axis
% also, convert it into NED frame 

yaw = atan2( T1(2,end), T1(1,end));
e3=[0 0 1]';
R = diag([1 -1 -1])*expm(-yaw*hat(e3));
for i=1:length(varnames)
    for k=1:length(T1)
        eval([varnames{i} '(:,k) = R*' varnames{i} '(:,k);']);
    end
end

figure;
plot3(T1(1,:),T1(2,:),T1(3,:));
set(gca,'YDir','reverse','ZDir','reverse');
title('transformed data');
axis equal;
grid on;

save('VICON_data', 't', 'LW1', 'LW2', 'LW3', 'RW1', 'RW2', 'RW3', 'T1', 'T2', 'A1', 'A2', 'LHW1', 'RHW1');
