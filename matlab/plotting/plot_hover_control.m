%% Figures for controlled hover
addpath('../modules', '../sim_data', '../');
set(0,'DefaultAxesFontName','times');
set(0,'DefaultAxesFontSize',18);
% load('sim_QS_x_hover_control.mat');

h_x=figure;
for ii=1:3 
    subplot(3,1,ii);
    plot(t*WK.f,x(ii,:));
    hold on;
    plot(t*WK.f,des.x_fit_t(ii, :), 'k');
    patch_downstroke(h_x,t*WK.f,Euler_R_dot);
end
xlabel('$t/T$','interpreter','latex');
subplot(3,1,2);
ylabel('$x$','interpreter','latex');
print(h_x, 'hover_control_pos', '-depsc', '-r0');

h_x3=figure;
plot3(x(1,:),x(2,:),x(3,:));
hold on;
plot3(des.x_fit_t(1,:),des.x_fit_t(2,:),des.x_fit_t(3,:), 'k');
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
zlabel('$x_3$','interpreter','latex');
axis equal;
xlim([-0.01, 0.01]);
zlim([-0.005, 0.005]);
print(h_x3, 'hover_control_pos_3d', '-depsc', '-r0');

h_x_dot=figure;
for ii=1:3 
    subplot(3,1,ii);
    plot(t*WK.f,x_dot(ii,:));
    hold on;
    plot(t*WK.f,des.x_dot_fit_t(ii, :), 'k');
    patch_downstroke(h_x_dot,t*WK.f,Euler_R_dot);
end
xlabel('$t/T$','interpreter','latex');
subplot(3,1,2);
ylabel('$\dot x$','interpreter','latex');
print(h_x_dot, 'hover_control_vel', '-depsc', '-r0');

h_err = figure;
subplot(3,1,1);
plot(t*WK.f, pos_err(1,:));
patch_downstroke(h_err,t*WK.f,Euler_R_dot);
subplot(3,1,2);
plot(t*WK.f, pos_err(2,:));
patch_downstroke(h_err,t*WK.f,Euler_R_dot);
ylabel('$\Delta x$','interpreter','latex');
subplot(3,1,3);
plot(t*WK.f, pos_err(3,:));
patch_downstroke(h_err,t*WK.f,Euler_R_dot);

des_cont = load('sim_QS_x_hover.mat', 'Euler_R', 'theta_A', 't');
des.Euler_R_fit_t = zeros(3, N);
des.theta_A_fit_t = fit(des_cont.t, des_cont.theta_A', 'fourier8');
des.theta_A_fit_t = des.theta_A_fit_t(t);
for i=1:3
    f = fit(des_cont.t, des_cont.Euler_R(i, :)', 'fourier8');
    des.Euler_R_fit_t(i,:) = f(t);
end
h_wk = figure;
subplot(4,1,1);
plot(t*WK.f, Euler_R(1,:) * 180/pi);
hold on;
plot(t*WK.f, des.Euler_R_fit_t(1,:) * 180/pi, 'k');
patch_downstroke(h_wk,t*WK.f,Euler_R_dot);
ylabel('$\phi$','interpreter','latex');
subplot(4,1,2);
plot(t*WK.f, Euler_R(2,:) * 180/pi);
hold on;
plot(t*WK.f, des.Euler_R_fit_t(2,:) * 180/pi, 'k');
patch_downstroke(h_wk,t*WK.f,Euler_R_dot);
ylabel('$\theta$','interpreter','latex');
subplot(4,1,3);
plot(t*WK.f, Euler_R(3,:) * 180/pi);
hold on;
plot(t*WK.f, des.Euler_R_fit_t(3,:) * 180/pi, 'k');
patch_downstroke(h_wk,t*WK.f,Euler_R_dot);
ylabel('$\psi$','interpreter','latex');
subplot(4,1,4);
plot(t*WK.f, theta_A * 180/pi);
hold on;
plot(t*WK.f, des.theta_A_fit_t * 180/pi, 'k');
patch_downstroke(h_wk,t*WK.f,Euler_R_dot);
ylabel('$\theta_A$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');
% print(h_wk, 'hover_control_wk', '-depsc', '-r0');

h_control = figure;
h_control.PaperUnits = 'inches';
h_control.PaperPosition = [0 0 6 5];
subplot(4,1,1);
plot(t*WK.f, dang(1,:));
patch_downstroke(h_control,t*WK.f,Euler_R_dot);
ylabel('$\Delta \phi_{m_s}$','interpreter','latex');
subplot(4,1,2);
plot(t*WK.f, dang(4,:));
patch_downstroke(h_control,t*WK.f,Euler_R_dot);
ylabel('$\Delta \phi_{m_k}$','interpreter','latex');
subplot(4,1,3);
plot(t*WK.f, dang(2,:));
patch_downstroke(h_control,t*WK.f,Euler_R_dot);
ylabel('$\Delta \theta_0$','interpreter','latex');
subplot(4,1,4);
plot(t*WK.f, dang(3,:));
patch_downstroke(h_control,t*WK.f,Euler_R_dot);
ylabel('$\Delta \theta_{A_m}$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');
print(h_control, 'hover_control_input', '-depsc', '-r0');
