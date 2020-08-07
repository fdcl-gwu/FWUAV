%% Figures for controlled hover
addpath('../modules', '../sim_data', '../');
set(0,'DefaultAxesFontName','times');
set(0,'DefaultAxesFontSize',18);
% load('sim_QS_x_hover_control.mat');

des_cont = load('sim_QS_xR_hover.mat', 'Euler_R', 'theta_A', 'f_a',...
    'Q_R', 'Q_L', 'M_R', 'M_L');
des_cont.fM_a = des_cont.f_a(1:6, :);
fM_a = f_a(1:6, :);
for i=1:N
    fM_a(4:6, i) = fM_a(4:6, i) + Q_R(:,:,i)*M_R(:,i) + Q_L(:,:,i)*M_L(:,i);
    des_cont.fM_a(4:6, i) = des_cont.fM_a(4:6, i) + des_cont.Q_R(:,:,i)*des_cont.M_R(:,i) +...
        des_cont.Q_L(:,:,i)*des_cont.M_L(:,i);
    [e_phi(i), e_psi(i), e_theta(i)] = dcm2angle(R(:,:,i)'*des.R_fit_t(:,:,i), 'xzy');
end

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
% print(h_x, 'hover_control_pos', '-depsc', '-r0');

% h_x3=figure;
% plot3(x(1,:),x(2,:),x(3,:));
% hold on;
% plot3(des.x_fit_t(1,:),des.x_fit_t(2,:),des.x_fit_t(3,:), 'k');
% xlabel('$x_1$','interpreter','latex');
% ylabel('$x_2$','interpreter','latex');
% zlabel('$x_3$','interpreter','latex');
% axis equal;
% xlim([-0.01, 0.01]);
% zlim([-0.005, 0.005]);
% set(gca,'Zdir','reverse','Xdir','reverse');
% print(h_x3, 'hover_control_pos_3d', '-depsc', '-r0');

h_rot=figure;
subplot(3,1,1);
plot(t*WK.f,e_phi);
patch_downstroke(h_rot,t*WK.f,Euler_R_dot);
ylabel('$d\phi$','interpreter','latex');
subplot(3,1,2);
plot(t*WK.f,e_psi);
patch_downstroke(h_rot,t*WK.f,Euler_R_dot);
ylabel('$d\psi$','interpreter','latex');
subplot(3,1,3);
plot(t*WK.f,e_theta);
patch_downstroke(h_rot,t*WK.f,Euler_R_dot);
ylabel('$d\theta$','interpreter','latex');
xlabel('$t/T$','interpreter','latex');
% print(h_rot, 'hover_control_rot', '-depsc', '-r0');

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
% print(h_x_dot, 'hover_control_vel', '-depsc', '-r0');

h_W=figure;
for ii=1:3 
    subplot(3,1,ii);
    plot(t*WK.f,W(ii,:));
    hold on;
    plot(t*WK.f,des.W_fit_t(ii, :), 'k');
    patch_downstroke(h_W,t*WK.f,Euler_R_dot);
end
xlabel('$t/T$','interpreter','latex');
subplot(3,1,2);
ylabel('$\Omega$','interpreter','latex');
% print(h_W, 'hover_control_ang_vel', '-depsc', '-r0');

% W_dot = diff(W, 1, 2);
% W_dot = [W_dot, W_dot(:,end)];
% des.W_dot = diff(des.W, 1, 2);
% des.W_dot = [des.W_dot, des.W_dot(:,end)];
% h_W=figure;
% for ii=1:3 
%     subplot(3,1,ii);
%     plot(t*WK.f,W_dot(ii,:));
%     hold on;
%     plot(t*WK.f,des.W_dot(ii, :), 'k');
%     patch_downstroke(h_W,t*WK.f,Euler_R_dot);
% end
% xlabel('$t/T$','interpreter','latex');
% subplot(3,1,2);
% ylabel('$\dot{\Omega}$','interpreter','latex');
% print(h_W, 'hover_control_ang_acc', '-depsc', '-r0');

% des_cont = load('sim_QS_x_hover.mat', 'Euler_R', 'theta_A', 't');
% des.Euler_R_fit_t = zeros(3, N);
% des.theta_A_fit_t = fit(des_cont.t, des_cont.theta_A', 'fourier8');
% des.theta_A_fit_t = des.theta_A_fit_t(t);
% for i=1:3
%     f = fit(des_cont.t, des_cont.Euler_R(i, :)', 'fourier8');
%     des.Euler_R_fit_t(i,:) = f(t);
% end
% h_wk = figure;
% subplot(4,1,1);
% plot(t*WK.f, Euler_R(1,:) * 180/pi);
% hold on;
% plot(t*WK.f, des.Euler_R_fit_t(1,:) * 180/pi, 'k');
% patch_downstroke(h_wk,t*WK.f,Euler_R_dot);
% ylabel('$\phi$','interpreter','latex');
% subplot(4,1,2);
% plot(t*WK.f, Euler_R(2,:) * 180/pi);
% hold on;
% plot(t*WK.f, des.Euler_R_fit_t(2,:) * 180/pi, 'k');
% patch_downstroke(h_wk,t*WK.f,Euler_R_dot);
% ylabel('$\theta$','interpreter','latex');
% subplot(4,1,3);
% plot(t*WK.f, Euler_R(3,:) * 180/pi);
% hold on;
% plot(t*WK.f, des.Euler_R_fit_t(3,:) * 180/pi, 'k');
% patch_downstroke(h_wk,t*WK.f,Euler_R_dot);
% ylabel('$\psi$','interpreter','latex');
% subplot(4,1,4);
% plot(t*WK.f, theta_A * 180/pi);
% hold on;
% plot(t*WK.f, des.theta_A_fit_t * 180/pi, 'k');
% patch_downstroke(h_wk,t*WK.f,Euler_R_dot);
% ylabel('$\theta_A$','interpreter','latex');
% xlabel('$t/T$','interpreter','latex');
% % print(h_wk, 'hover_control_wk', '-depsc', '-r0');

h_err = figure;
h_err.PaperUnits = 'inches';
h_err.PaperPosition = [0 0 6 8];
for ii=1:6
    subplot(6,1,ii);
    plot(t*WK.f, err_xR(ii,:));
    patch_downstroke(h_err,t*WK.f,Euler_R_dot);
end
xlabel('$t/T$','interpreter','latex');
subplot(6,1,3);
ylabel('Error');

h_control = figure;
h_control.PaperUnits = 'inches';
h_control.PaperPosition = [0 0 6 8];
ylabels = {'$\Delta\phi_{m_s}$', '$\Delta\theta_{0_s}$', '$\Delta\phi_{m_k}$',...
        '$\Delta\phi_{0_s}$', '$\Delta\theta_{0_k}$', '$\Delta\psi_{0_k}$'};
for ii=1:6
    subplot(6,1,ii);
    plot(t*WK.f, dang(ii,:));
    patch_downstroke(h_control,t*WK.f,Euler_R_dot);
    ylabel(ylabels{ii},'interpreter','latex');
end
xlabel('$t/T$','interpreter','latex');
% print(h_control, 'hover_control_input', '-depsc', '-r0');

h_aero = figure;
h_aero.PaperUnits = 'inches';
h_aero.PaperPosition = [0 0 6 8];
for ii=1:6
    subplot(6,1,ii);
    plot(t*WK.f, fM_a(ii,:));
    hold on;
    plot(t*WK.f, des_cont.fM_a(ii, :), 'k');
    patch_downstroke(h_aero,t*WK.f,Euler_R_dot);
end
xlabel('$t/T$','interpreter','latex');
subplot(6,1,3);
ylabel('F, M');
