clear all;
close all;

load sim_QS_x_wo_ab;

h_x3=figure;
plot(x(1,:),x(3,:),'g');
hold on;

h_x=figure;
for ii=1:2
    subplot(2,1,ii);
    plot(t*WK.f,x(2*ii-1,:),'g');
    hold on;
end

h_v=figure;
for ii=1:2
    subplot(2,1,ii);
    plot(t*WK.f,x_dot(2*ii-1,:),'g');
    hold on;
end

h_v3=figure;
plot(x_dot(1,:),x_dot(3,:),'g');
hold on;

h_theta=figure;
subplot(2,1,1);
plot(t*WK.f,theta_B*180/pi,'g');
hold on;
subplot(2,1,2);
plot(t*WK.f,theta_A*180/pi,'g');
hold on;



%%
load sim_QS_x_with_ab;

figure(h_x3);
plot(x(1,:),x(3,:),'b');
set(gca,'YDir','reverse');
hold on;
xlabel('$x_1$','interpreter','latex');
ylabel('$x_3$','interpreter','latex');

figure(h_x);
for ii=1:2
    subplot(2,1,ii);
    plot(t*WK.f,x(2*ii-1,:),'b');
end
xlabel('$t/T$','interpreter','latex');
subplot(2,1,1);
ylabel('$x_1$','interpreter','latex');
subplot(2,1,2);
ylabel('$x_3$','interpreter','latex');


figure(h_v);
for ii=1:2
    subplot(2,1,ii);
    plot(t*WK.f,x_dot(2*ii-1,:),'b');
end
xlabel('$t/T$','interpreter','latex');
subplot(2,1,1);
ylabel('$\dot x_1$','interpreter','latex');
subplot(2,1,2);
ylabel('$\dot x_3$','interpreter','latex');


figure(h_v3);
plot(x_dot(1,:),x_dot(3,:),'b');
set(gca,'YDir','reverse');
axis equal;

xlabel('$\dot x_1$','interpreter','latex');
ylabel('$\dot x_3$','interpreter','latex');

figure(h_theta);
subplot(2,1,1);
plot(t*WK.f,theta_B*180/pi,'b');
ylabel('$\theta_B$','interpreter','latex');
subplot(2,1,2);
plot(t*WK.f,theta_A*180/pi,'b');
ylabel('$\theta_A$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');


%%
load sim_QS_x_with_opposite_ab;

figure(h_x3);
plot(x(1,:),x(3,:),'m');

figure(h_x);
for ii=1:2
    subplot(2,1,ii);
    plot(t*WK.f,x(2*ii-1,:),'m');
end

figure(h_v);
for ii=1:2
    subplot(2,1,ii);
    plot(t*WK.f,x_dot(2*ii-1,:),'m');
    patch_downstroke(h_v,t*WK.f,Euler_R_dot);
    xlim([0 5]);
end

figure(h_v3);
plot(x_dot(1,:),x_dot(3,:),'m');
set(gca,'YDir','reverse');
axis equal;

figure(h_theta);
subplot(2,1,1);
plot(t*WK.f,theta_B*180/pi,'m');
xlim([0 5]);
patch_downstroke(h_theta,t*WK.f,Euler_R_dot);
subplot(2,1,2);
plot(t*WK.f,theta_A*180/pi,'m');
xlim([0 5]);
patch_downstroke(h_theta,t*WK.f,Euler_R_dot);



%%

bool_print=true;
figname='comp_ab';
if bool_print
    figure(h_x3);print([figname '_x3'],'-depsc');
    figure(h_x);print([figname '_x'],'-depsc');
    figure(h_v);print([figname '_v'],'-depsc');
    figure(h_v3);print([figname '_v3'],'-depsc');
    figure(h_theta);print([figname '_theta'],'-depsc');
end
!mv comp_ab*.eps ../doc/Figs/