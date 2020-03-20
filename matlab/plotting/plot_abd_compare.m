%% Comparison to observe the effect of abdomen oscillation
addpath('../modules', '../sim_data', '../');
set(0,'DefaultAxesFontName','times');
set(0,'DefaultAxesFontSize',16);

%% Monte carlo study
load('sim_QS_x_hover_control_mc_pos_long_asy.mat');
err_bound = 1e-4; % 5e-4
% err_bound = squeeze(vecnorm(x_pert, 2, 1))';
idx_no_ab = (abs(err_pos(:, 1)) < err_bound);
x_no_ab = x_pert(:, idx_no_ab);
bound_ix_no_ab = boundary(x_no_ab(1, :)', x_no_ab(3, :)');
idx_with_ab = (abs(err_pos(:, 2)) < err_bound);
x_with_ab = x_pert(:, idx_with_ab);
bound_ix_with_ab = boundary(x_with_ab(1, :)', x_with_ab(3, :)');
x_with_ab_only = x_pert(:, (abs(err_pos(:, 2)) < err_bound) & (abs(err_pos(:, 1)) > err_bound));

% fac = 2;
% plot_ix = abs(x_no_ab(1, :)) < fac*cov_no_ab(1) & abs(x_no_ab(3, :)) < fac*cov_no_ab(3);
h_err = figure;

% % y perturbation
% scatter(x_no_ab(1, :), x_no_ab(2, :), 20, 'r', 'x');
% hold on;
% scatter(x_with_ab(1, :), x_with_ab(2, :), 10, 'b', 'filled');

% % x-z perturbation
scatter(x_no_ab(1, :), x_no_ab(3, :), 20, 'r', 'x');
hold on;
scatter(x_with_ab_only(1, :), x_with_ab_only(3, :), 10, 'b', 'filled');
plot(x_no_ab(1, bound_ix_no_ab), x_no_ab(3, bound_ix_no_ab), 'r', 'LineWidth', 2);
plot(x_with_ab(1, bound_ix_with_ab), x_with_ab(3, bound_ix_with_ab), 'b', 'LineWidth', 2);
legend({'without abdomen effect, $ w = 0 $', 'with abdomen effect, $ w = 0.1 $'},...
    'interpreter','latex','FontSize',12);
xlabel('$x$','interpreter','latex');
ylabel('$z$','interpreter','latex');
axis tight;
print(h_err, 'hover_mc_roa', '-depsc');

N_no_ab = N_conv(idx_no_ab, 1);
N_with_ab = N_conv(idx_no_ab, 2);
h_Nconv = figure;
scatter(vecnorm(x_no_ab, 2, 1), (N_with_ab - N_no_ab), 10, 'filled');
xlabel('Norm of x-z perturbation');
ylabel('Difference in number of cycles to convergence', 'FontSize',12);
print(h_Nconv, 'hover_mc_perf', '-depsc');
%
% N_mesh = 100;
% [x_mesh, z_mesh] = meshgrid(linspace(min(x_no_ab(1,:)),max(x_no_ab(1,:)),N_mesh),...
%     linspace(min(x_no_ab(3,:)),max(x_no_ab(3,:)),N_mesh));
% N_diff_mesh = griddata(x_no_ab(1,:),x_no_ab(3,:), (N_with_ab-N_no_ab), x_mesh,z_mesh, 'natural');
% h_Nconv = figure;
% surf(x_mesh, z_mesh, N_diff_mesh, 'FaceColor', 'blue', 'FaceAlpha', 0.5, 'EdgeColor', 'none');

%% Power and energy
% Without abdomen oscillation
% load('sim_QS_x_hover_ab_no.mat')
load('sim_QS_x_hover_no_ab.mat')
tau_no = tau(4:12,:);
[pow_no, E_no, E_dot_no, eff_no] = compute_power(INSECT.m, t, x, x_dot, tau, Q_R, Q_L, Q_A, W_R, W_L, W_A);

% With abdomen oscillation
% load('sim_QS_x_hover.mat')
load('sim_QS_x_hover.mat')
tau_ab = tau(4:12,:);
[pow_ab, E_ab, E_dot_ab, eff_ab] = compute_power(INSECT.m, t, x, x_dot, tau, Q_R, Q_L, Q_A, W_R, W_L, W_A);

N_period = T*WK.f;
mean_pow_no = sum(mean(abs(pow_no), 2) / N_period);
mean_pow_ab = sum(mean(abs(pow_ab), 2) / N_period);
change_pow = (mean_pow_ab - mean_pow_no) ./ mean_pow_no;
h_pow = figure;
subplot(2,1,1);
plot(t*WK.f,pow_no(1,:),'r');
hold on;
plot(t*WK.f,pow_ab(1,:),'b');
patch_downstroke(h_pow,t*WK.f,Euler_R_dot);
ylabel('$P_R$','interpreter','latex');
subplot(2,1,2);
plot(t*WK.f,pow_no(3,:),'r');
hold on;
plot(t*WK.f,pow_ab(3,:),'b');
patch_downstroke(h_pow,t*WK.f,Euler_R_dot);
ylabel('$P_A$','interpreter','latex');
sgtitle('Reduction in total mean power is ' + string(round(-change_pow*100, 1)) + ' %');
% print(h_pow, 'hover_power_ab', '-depsc');

mean_tau_no = mean(vecnorm(tau_no(1:3,:), 2, 1), 2);
mean_tau_ab = mean(vecnorm(tau_ab(1:3,:), 2, 1), 2);
change_tau = (mean_tau_ab - mean_tau_no) ./ mean_tau_no;
h_tau = figure;
subplot(2,1,1);
plot(t*WK.f,vecnorm(tau_no(1:3,:), 2, 1),'r');
hold on;
plot(t*WK.f,vecnorm(tau_ab(1:3,:), 2, 1),'b');
patch_downstroke(h_tau,t*WK.f,Euler_R_dot);
ylabel('$\tau_R$','interpreter','latex');
subplot(2,1,2);
plot(t*WK.f,vecnorm(tau_no(7:9,:), 2, 1),'r');
hold on;
plot(t*WK.f,vecnorm(tau_ab(7:9,:), 2, 1),'b');
patch_downstroke(h_tau,t*WK.f,Euler_R_dot);
ylabel('$\tau_A$','interpreter','latex');
% print(h_tau, 'hover_torque_ab', '-depsc');

mean_E_no = mean(abs(E_no), 2) / N_period;
mean_E_ab = mean(abs(E_ab), 2) / N_period;
change_E = (mean_E_ab - mean_E_no) ./ mean_E_no;
h_E = figure;
subplot(2,1,1);
plot(t*WK.f,E_no,'r');
patch_downstroke(h_E,t*WK.f,Euler_R_dot);
hold on;
plot(t*WK.f,E_ab,'b');
ylabel('$E$','interpreter','latex');
subplot(2,1,2);
plot(t*WK.f,E_dot_no,'r');
patch_downstroke(h_E,t*WK.f,Euler_R_dot);
hold on;
plot(t*WK.f,E_dot_ab,'b');
ylabel('$\dot{E}$','interpreter','latex');
sgtitle('Reduction in total energy is ' + string(round(-change_E*100, 1)) + ' %');
% print(h_E, 'hover_energy_ab', '-depsc');
