function load_exp_data
clear all;
close all;
addpath('../');
e2=[0 1 0]';
e3=[0 0 1]';
g=9.81;
m = 0.5077e-3;

XLS=readtable('Time_history_data_for_Dr.Lee.xlsx');
XLS=table2array(XLS);

%%
I_cycle=[1 21;
    21 41;
    41 61;
    61 80;
    80 99];

N_cycle=5;
x=cell(N_cycle,1);
v=cell(N_cycle,1);
t=cell(N_cycle,1);
pitch_th=cell(N_cycle,1);
pitch_ab=cell(N_cycle,1);
for i=1:N_cycle
    I_range=I_cycle(i,1):I_cycle(i,2);
    t{i} = XLS(I_range,1);    
    t{i} = t{i}-t{i}(1);
    x{i} = XLS(I_range,2:4);    
    v{i} = XLS(I_range,5:7);    
    
    pitch_th{i} = XLS(I_range,8);
    pitch_ab{i} = XLS(I_range,9);
    
    phi{i} = XLS(I_range,10);
    theta{i} = XLS(I_range,11);
    psi{i} = XLS(I_range,12);
    
    CL{i} = XLS(I_range,13);
    CT{i} = XLS(I_range,14);
    
    for j=1:3
        a{i}(:,j) = diff(v{i}(:,j))./diff(t{i});
    end
    
    for k=1:size(a{i},1)
        F_body{i}(k,:) = ( diag([1 1 -1]) * expm(pitch_th{i}(k)*pi/180*hat(e2))*a{i}(k,:)' - g*e3)'*m/2 ;
        a{i}(k,:) = (diag([1 1 -1]) * a{i}(k,:)')';
        v{i}(k,:) = (diag([1 1 -1]) * v{i}(k,:)')';
    end
    
end

h_x=figure;
for i=1:N_cycle
    plot(t{i}/max(t{i}),x{i});
    hold on;
end

h_v=figure;
for i=1:N_cycle
    for j=1:3
        subplot(3,1,j);
        plot(t{i}/max(t{i}),v{i}(:,j),'b.');
        hold on;
    end    
end
xlabel('$t$','interpreter','latex');
subplot(3,1,2);
ylabel('$v$','interpreter','latex');

h_a=figure;
for i=1:N_cycle
    for j=1:3
        subplot(3,1,j);
        plot(t{i}(2:end)/max(t{i}),a{i}(:,j),'b.');
        hold on;
    end    
end
xlabel('$t$','interpreter','latex');
subplot(3,1,2);
ylabel('$a$','interpreter','latex');

h_F_body=figure;
for i=1:N_cycle
    for j=1:3
        subplot(3,1,j);
        plot(t{i}(2:end)/max(t{i}),F_body{i}(:,j),'b.');
        hold on;
    end    
end
xlabel('$t$','interpreter','latex');
subplot(3,1,2);
ylabel('$F_B$','interpreter','latex');

load('fit_pitch_thorax');
h_pitch=figure;
N=501;
t_fit=linspace(0,1/F_theta_th.f,N);
for k=1:N
    [~,~,~,theta_th_fit(k)]=body_attitude(t_fit(k));
end

for i=1:N_cycle
    subplot(2,1,1);
    plot(t{i}/max(t{i}),pitch_th{i},'b.');  
    hold on;
    plot(t_fit/t_fit(end),theta_th_fit*180/pi,'r');
    ylabel('$\theta_{th}$','interpreter','latex');
    subplot(2,1,2);
    plot(t{i}/max(t{i}),pitch_ab{i},'b.');
    hold on;
    ylabel('$\theta_{ab}$','interpreter','latex');    
    
end
xlabel('$t$','interpreter','latex');

keyboard;
return;

h_E=figure;
for i=1:N_cycle
    subplot(3,1,1);
    plot(t{i}/max(t{i}),phi{i},'b.');
    hold on;
    ylabel('$\phi$','interpreter','latex');
    
    subplot(3,1,2);
    plot(t{i}/max(t{i}),theta{i},'b.');
    hold on;
    ylabel('$\theta$','interpreter','latex');
    
    subplot(3,1,3);
    plot(t{i}/max(t{i}),psi{i},'b.');
    hold on;
    ylabel('$\psi$','interpreter','latex');    
end
xlabel('$t$','interpreter','latex');

h_C=figure;
for i=1:N_cycle
    subplot(2,1,1);
    plot(t{i}/max(t{i}),CL{i},'b.');
    hold on;
    ylabel('$C_L$','interpreter','latex');
    subplot(2,1,2);
    plot(t{i}/max(t{i}),CT{i},'b.');
    hold on;
    ylabel('$C_T$','interpreter','latex');    
    
end
xlabel('$t$','interpreter','latex');


bool_print=0;

filename='exp';
if bool_print
    print(h_v,[filename '_v'],'-depsc');
    print(h_a,[filename '_a'],'-depsc');
    print(h_F_body,[filename '_F_B'],'-depsc');
    print(h_pitch, [filename '_pitch'],'-depsc');
    print(h_E, [filename '_E'],'-depsc');
    print(h_E, [filename '_C'],'-depsc');
    
    evalin('base',['!mv ' filename '*.eps ../../doc/Figs']);
end

% Get a list of all variables
allvars = whos;
% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics object
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
% Pass these variable names to save
save(filename, allvars(tosave).name)
%evalin('base',['load ' filename]);
end

