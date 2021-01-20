clear all;
close all;

set(0,'DefaultAxesFontName','times');
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultLineLineWidth',2);

XLIM=[0 5];

load sim_QS_x;

h_x3=figure;
plot3(x(1,:),x(2,:),x(3,:));
set(gca,'YDir','reverse','ZDir','reverse');
axis equal;

h_v3=figure;
plot(x_dot(1,:),x_dot(3,:));
% set(gca,'YDir','reverse');
hold on
axis equal;
scatter(x_dot(1,1),x_dot(3,1),50,'r','filled');
text(x_dot(1,1),x_dot(3,1), '$\leftarrow v(0)$', 'FontSize', 18, 'Interpreter', 'latex')
xlabel('$v_1$','interpreter','latex');
ylabel('$v_3$','interpreter','latex');

h_x=figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t*WK.f,x(ii,:));
    xlim(XLIM);
end

h_v=figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t*WK.f,x_dot(ii,:));
    xlim(XLIM);
end

h_E=figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t*WK.f,Euler_R(ii,:)*180/pi);
    xlim(XLIM);
end

h_theta=figure;
subplot(2,1,1);
plot(t*WK.f,theta_B*180/pi);
xlim(XLIM);
subplot(2,1,2);
plot(t*WK.f,theta_A*180/pi);
xlim(XLIM);

h_FB=figure;
for ii=1:3
    subplot(3,1,ii);
    plot(t*WK.f,F_B(ii,:));
    xlim(XLIM);
end

h_wingtip=figure;
wing_tip_R=squeeze(Q_R(:,2,:));
wing_tip_L=squeeze(-Q_L(:,2,:));
plot3(wing_tip_R(1,:),wing_tip_R(2,:),wing_tip_R(3,:),'b');
hold on;
plot3(wing_tip_L(1,:),wing_tip_L(2,:),wing_tip_L(3,:),'g');
set(gca,'YDir','reverse','ZDir','reverse');
axis equal;
grid on;
xlabel('$\mathbf{b}_x$','interpreter','latex');
ylabel('$\mathbf{b}_y$','interpreter','latex');
zlabel('$\mathbf{b}_z$','interpreter','latex');



t_sim_QS=t;

%%
load exp_data/fit_VICON_data.mat

x=(T1+T2)/2;
v=(diff(x')./(diff(t)*ones(1,3)))';
a=(diff(v')./(diff(t(2:end))*ones(1,3)))';

g=9.81;
e2=[0 1 0]';
e3=[0 0 1]';
F_B=[];
for k=1:length(a)
    F_B(:,k) = ( expm(theta_th(k)*pi/180*hat(e2))*a(:,k) - g*e3)'*MONARCH.m ;
end
% 
        
figure(h_x3);
hold on;
plot3(x(1,:),x(2,:),x(3,:),'r');
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
zlabel('$x_3$','interpreter','latex');

figure(h_x);
for ii=1:3
    subplot(3,1,ii);
    hold on;
    plot(t*WK.f,x(ii,:),'r');
    ylabel('$x_'+string(ii)+'$','interpreter','latex');
    patch_downstroke(h_x,t_sim_QS*WK.f,Euler_R_dot);
end
xlabel('$t/T$','interpreter','latex');


figure(h_E);
E_ylabel={'$\phi$','$\theta$','$\psi$'};
E=0.5*(E_R+E_L)*180/pi;
for ii=1:3
    subplot(3,1,ii);
    hold on;
    plot(t*WK.f,E(ii,:),'r');
    patch_downstroke(h_E,t_sim_QS*WK.f,Euler_R_dot);
    ylabel(E_ylabel{ii},'interpreter','latex');
end
xlabel('$t/T$','interpreter','latex');

figure(h_theta);
subplot(2,1,1);
hold on;
plot(t*WK.f,theta_th*180/pi,'r');
patch_downstroke(h_theta,t_sim_QS*WK.f,Euler_R_dot);
ylabel('$\theta_B$','interpreter','latex');
subplot(2,1,2);
hold on;
plot(t*WK.f,theta_ab*180/pi,'r');
patch_downstroke(h_theta,t_sim_QS*WK.f,Euler_R_dot);
ylabel('$\theta_A$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');


figure(h_v);
for ii=1:3
    subplot(3,1,ii);
    hold on;
    plot(t(2:end)*WK.f,v(ii,:),'r');
    ylabel('$\dot{x}_'+string(ii)+'$','interpreter','latex');
    patch_downstroke(h_v,t_sim_QS*WK.f,Euler_R_dot);
end
xlabel('$t/T$','interpreter','latex');


figure(h_FB);
for ii=1:3
    subplot(3,1,ii);
    hold on;
    plot(t(3:end)*WK.f,F_B(ii,:),'r');
    ylabel('${F_B}_'+string(ii)+'$','interpreter','latex');
    patch_downstroke(h_FB,t_sim_QS*WK.f,Euler_R_dot);
end
xlabel('$t/T$','interpreter','latex');

h_E = figure;
my_ylabel={'$\phi$','$\theta$','$\psi$'};
for ii=1:3
    subplot(3,1,ii);
    plot(t*WK.f,E_R(ii,:)*180/pi,'r','LineWidth',1.5);
    hold on;
    plot(t*WK.f,E_L(ii,:)*180/pi,'b','LineWidth',1.5);
    ylabel(my_ylabel{ii},'interpreter','latex');
    plot(t*WK.f,0.5*(E_R(ii,:)+E_L(ii,:))*180/pi,'k--','LineWidth',1.5);
    axis('tight');
    if ii==1
        legend({'right wing','left wing','average'});
    end
end
xlabel('$t/T$','interpreter','latex');

bool_print=true;
filename='comp_VICON';
if bool_print
    figure(h_x3);pause(1);print([filename '_x3'],'-depsc2');
    figure(h_v3);pause(1);print([filename '_v3'],'-depsc2');
    figure(h_x);pause(1);print([filename '_x'],'-depsc2');
    figure(h_theta);pause(1);print([filename '_theta'],'-depsc2');
    figure(h_E);pause(1);print([filename '_E'],'-depsc2');
    figure(h_v);pause(1);print([filename '_v'],'-depsc2');
    figure(h_FB);pause(1);print([filename '_FB'],'-depsc2');
    figure(h_E);pause(1);print('fit_VICON_E','-depsc2');
end
% !mv *.eps ../doc/Figs

