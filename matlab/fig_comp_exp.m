clear all;
close all;

set(0, 'DefaultLineLineWidth', 2);

load sim_QS_x;

h_x3=figure;
plot3(x(1,:),x(2,:),x(3,:));
set(gca,'YDir','reverse','ZDir','reverse');
axis equal;

h_x=figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t*WK.f,x(ii,:));
end

h_v=figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t*WK.f,x_dot(ii,:));
end

h_E=figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t*WK.f,Euler_R(ii,:)*180/pi);
end

h_theta=figure;
subplot(2,1,1);
plot(t*WK.f,theta_B*180/pi);
subplot(2,1,2);
plot(t*WK.f,theta_A*180/pi);

h_FB=figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t*WK.f,F_B(ii,:));
end

t_sim_QS=t;

%%
XLS=readtable('./exp_data/Time_history_data_for_Dr.Lee.xlsx');
XLS=table2array(XLS);

t = XLS(:,1)-XLS(1,1);    
x = XLS(:,2:4);    
v = XLS(:,5:7);    
x(:,3)=-x(:,3);
v(:,3)=-v(:,3);
a= diff(v)./diff(t);

theta_th = XLS(:,8);
theta_ab = XLS(:,9)-theta_th;

phi = XLS(:,15);
theta = - XLS(:,16);
psi = XLS(:,17);
E=[phi, theta, psi];

g=9.81;
e2=[0 1 0]';
e3=[0 0 1]';
F_B=[];
for k=1:length(a)
    F_B(k,:) = ( expm(theta_th(k)*pi/180*hat(e2))*a(k,:)' - g*e3)'*MONARCH.m ;
end
% 
        
figure(h_x3);
hold on;
plot3(x(:,1),x(:,2),x(:,3),'r');
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
zlabel('$x_3$','interpreter','latex');

figure(h_x);
for ii=1:3
    subplot(3,1,ii);
    hold on;
    plot(t*WK.f,x(:,ii),'r');
    patch_downstroke(h_x,t_sim_QS*WK.f,Euler_R_dot);
end
xlabel('$t/T$','interpreter','latex');
subplot(3,1,2);
ylabel('$x$','interpreter','latex');


figure(h_E);
E_ylabel={'$\phi$','$\theta$','$\psi$'};
for ii=1:3
    subplot(3,1,ii);
    hold on;
    plot(t*WK.f,E(:,ii),'r');
    patch_downstroke(h_E,t_sim_QS*WK.f,Euler_R_dot);
    ylabel(E_ylabel{ii},'interpreter','latex');
end
xlabel('$t/T$','interpreter','latex');

figure(h_theta);
subplot(2,1,1);
hold on;
plot(t*WK.f,theta_th,'r');
patch_downstroke(h_theta,t_sim_QS*WK.f,Euler_R_dot);
ylabel('$\theta_B$','interpreter','latex');
subplot(2,1,2);
hold on;
plot(t*WK.f,theta_ab,'r');
patch_downstroke(h_theta,t_sim_QS*WK.f,Euler_R_dot);
ylabel('$\theta_A$','interpreter','latex');


figure(h_v);
for ii=1:3
    subplot(3,1,ii);
    hold on;
    plot(t*WK.f,v(:,ii),'r');
    patch_downstroke(h_v,t_sim_QS*WK.f,Euler_R_dot);
end
xlabel('$t/T$','interpreter','latex');
subplot(3,1,2);
ylabel('$\dot x$','interpreter','latex');


figure(h_FB);
for ii=1:3
    subplot(3,1,ii);
    hold on;
    plot(t(2:end)*WK.f,F_B(:,ii),'r');
    patch_downstroke(h_FB,t_sim_QS*WK.f,Euler_R_dot);
end
xlabel('$t/T$','interpreter','latex');
subplot(3,1,2);
ylabel('$F_B$','interpreter','latex');


bool_print=false;
filename='exp_comp';
if bool_print
    figure(h_x3);pause(1);print([filename '_x3'],'-depsc2');
    figure(h_x);pause(1);print([filename '_x'],'-depsc2');
    figure(h_theta);pause(1);print([filename '_theta'],'-depsc2');
    figure(h_E);print([filename '_E'],'-depsc2');
    figure(h_v);print([filename '_v'],'-depsc2');
    figure(h_FB);print([filename '_FB'],'-depsc2');
end
!mv *.eps ../doc/Figs

