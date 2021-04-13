%% File to plot the data

% % print options
% print(h_x_dot, 'hover_control_vel', '-depsc');
% print(h_x_dot, 'hover_control_vel', '-depsc', '-r0');
% print(h_x_dot, 'hover_control_vel', '-depsc', '-r300');
% print(h_x_dot, '-painters', 'hover_control_vel', '-depsc');

addpath('../modules', '../sim_data', '../');
set(0,'DefaultAxesFontName','times');
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultLineLineWidth',1.5);

h_x3=figure;
plot3(x(1,:),x(2,:),x(3,:));
% set(gca,'YDir','reverse');
% set(gca,'fontname','times')
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
zlabel('$x_3$','interpreter','latex');
axis equal;
print(h_x3, 'hover_position_3d', '-depsc');

h_x=figure;
h_x.PaperUnits = 'inches';
h_x.PaperPosition = [0 0 8 8];
for ii=1:3 
    subplot(3,1,ii);
    plot(t*WK.f,x(ii,:));
end
xlabel('$t/T$','interpreter','latex');
subplot(3,1,2);
ylabel('$x$','interpreter','latex');

h_x_dot=figure;
for ii=1:3 
    subplot(3,1,ii);
    plot(t*WK.f,x_dot(ii,:));
    patch_downstroke(h_x_dot,t*WK.f,Euler_R_dot);
end
xlabel('$t/T$','interpreter','latex');
subplot(3,1,2);
ylabel('$\dot x$','interpreter','latex');
print(h_x_dot, 'hover_velocity', '-depsc');

nn = numel(t);
% Euler = zeros(3, nn);
% for j=1:nn
%     Euler(:, j) = wing_kinematics(t(j),WK) * 180/pi;
% end
h_wk = figure;
h_wk.PaperUnits = 'inches';
h_wk.PaperPosition = [0 0 6 6];
subplot(4,1,1);
plot(t*WK.f, Euler_R(1,:) * 180/pi);
patch_downstroke(h_wk,t*WK.f,Euler_R_dot);
ylabel('$\phi$','interpreter','latex');
subplot(4,1,2);
plot(t*WK.f, Euler_R(2,:) * 180/pi);
patch_downstroke(h_wk,t*WK.f,Euler_R_dot);
ylabel('$\theta$','interpreter','latex');
subplot(4,1,3);
plot(t*WK.f, Euler_R(3,:) * 180/pi);
patch_downstroke(h_wk,t*WK.f,Euler_R_dot);
ylabel('$\psi$','interpreter','latex');
subplot(4,1,4);
plot(t*WK.f, theta_A*180/pi);
patch_downstroke(h_wk,t*WK.f,Euler_R_dot);
ylabel('$\theta_A$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');
print(h_wk, 'hover_wing_kinematics', '-depsc');

h_aero = figure;
subplot(3,1,1);
plot(t*WK.f, f_a(1,:));
axis tight;
% patch_downstroke(h_aero,t*WK.f,Euler_R_dot);
subplot(3,1,2);
plot(t*WK.f, f_a(2,:));
axis tight;
% patch_downstroke(h_aero,t*WK.f,Euler_R_dot);
ylabel('$f_a$','interpreter','latex');
subplot(3,1,3);
plot(t*WK.f, f_a(3,:));
axis tight;
% patch_downstroke(h_aero,t*WK.f,Euler_R_dot);
xlabel('$t/T$','interpreter','latex');

h_ab = figure;
subplot(2,1,1);
plot(t*WK.f, theta_B*180/pi);
patch_downstroke(h_ab,t*WK.f,Euler_R_dot);
ylabel('$\theta_B$','interpreter','latex');
subplot(2,1,2);
plot(t*WK.f, W(2,:));
patch_downstroke(h_ab,t*WK.f,Euler_R_dot);
ylabel('$\langle \Omega, e_2 \rangle$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');
print(h_ab, 'hover_abd_body', '-depsc');

[pow, E, E_dot, eff] = compute_power(INSECT.m, t, x, x_dot, tau, Q_R, Q_L, Q_A, W_R, W_L, W_A, W, f_a, f_tau);

h_pow = figure;
sgtitle('Power');
subplot(3,1,1);
plot(t*WK.f,pow(1,:));
patch_downstroke(h_pow,t*WK.f,Euler_R_dot);
ylabel('$P_R$','interpreter','latex');
hold on;
subplot(3,1,2);
plot(t*WK.f,pow(2,:));
patch_downstroke(h_pow,t*WK.f,Euler_R_dot);
ylabel('$P_L$','interpreter','latex');
hold on;
subplot(3,1,3);
plot(t*WK.f,pow(3,:));
patch_downstroke(h_pow,t*WK.f,Euler_R_dot);
ylabel('$P_A$','interpreter','latex');
hold on;
% print(h_pow, 'power_abdomen_osc', '-depsc');

h_E = figure;
subplot(3,1,1);
plot(t*WK.f,E);
patch_downstroke(h_E,t*WK.f,Euler_R_dot);
ylabel('$E$','interpreter','latex');
hold on;
subplot(3,1,2);
plot(t*WK.f,E_dot);
patch_downstroke(h_E,t*WK.f,Euler_R_dot);
ylabel('$\dot{E}$','interpreter','latex');
hold on;
subplot(3,1,3);
plot(t*WK.f,eff);
patch_downstroke(h_E,t*WK.f,Euler_R_dot);
ylabel('$eff$','interpreter','latex');
hold on;
% print(h_E, 'energy_abdomen_osc', '-depsc');
