%% Comparison to observe the effect of abdomen oscillation
addpath('../modules', '../sim_data', '../');

%% Monte carlo study
load('sim_QS_x_hover_control_monte_carlo_pos_longitudinal.mat');
err_bound = 5e-4;
% err_bound = squeeze(vecnorm(x_pert, 2, 1))';
x_no_ab = x_pert(:, abs(err_pos(:, 1)) > err_bound);
x_with_ab = x_pert(:, abs(err_pos(:, 2)) > err_bound);
% cov_no_ab = diag(cov(x_no_ab'));
% x_cov = linspace(-cov_no_ab(1), cov_no_ab(1), 50);
% z_cov = cov_no_ab(3) * sqrt(1 - (x_cov/cov_no_ab(1)).^2);
% cov_with_ab = diag(cov(x_with_ab'));
% x_cov_ab = linspace(-cov_with_ab(1), cov_with_ab(1), 50);
% z_cov_ab = cov_with_ab(3) * sqrt(1 - (x_cov_ab/cov_with_ab(1)).^2);

% fac = 2;
h_err = figure;
% plot_ix = abs(x_no_ab(1, :)) < fac*cov_no_ab(1) & abs(x_no_ab(3, :)) < fac*cov_no_ab(3);
scatter(x_no_ab(1, :), x_no_ab(3, :), 20, 'r', 'x');
hold on;
% plot_ix = abs(x_with_ab(1, :)) < fac*cov_with_ab(1) & abs(x_with_ab(3, :)) < fac*cov_with_ab(3);
scatter(x_with_ab(1, :), x_with_ab(3, :), 10, 'b', 'filled');
% hold on;
% plot(x_cov, z_cov, 'r', 'LineWidth', 2);
% hold on;
% plot(x_cov, -z_cov, 'r', 'LineWidth', 2);
% hold on;
% plot(x_cov_ab, z_cov_ab, 'b', 'LineWidth', 2);
% hold on;
% plot(x_cov_ab, -z_cov_ab, 'b', 'LineWidth', 2);
legend({'without abdomen effect, $ w = 0 $', 'with abdomen effect, $ w = 0.1 $'},'interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$z$','interpreter','latex');
axis tight;
% print(h_err, 'hover_monte_carlo', '-depsc');

%% Power and energy
% Without abdomen oscillation
% load('sim_QS_x_hover_ab_no.mat')
load('sim_QS_x_hover_control_no_ab.mat')
[pow_no, E_no, E_dot_no, eff_no] = compute_power(INSECT.m, t, x, x_dot, tau, Q_R, Q_L, Q_A, W_R, W_L, W_A);

% With abdomen oscillation
% load('sim_QS_x_hover.mat')
load('sim_QS_x_hover_control_with_ab.mat')
[pow_ab, E_ab, E_dot_ab, eff_ab] = compute_power(INSECT.m, t, x, x_dot, tau, Q_R, Q_L, Q_A, W_R, W_L, W_A);

h_pow = figure;
subplot(3,1,1);
plot(t*WK.f,pow_no(1,:),'r');
hold on;
plot(t*WK.f,pow_ab(1,:),'b');
patch_downstroke(h_pow,t*WK.f,Euler_R_dot);
ylabel('$P_R$','interpreter','latex');
subplot(3,1,2);
plot(t*WK.f,pow_no(2,:),'r');
hold on;
plot(t*WK.f,pow_ab(2,:),'b');
patch_downstroke(h_pow,t*WK.f,Euler_R_dot);
ylabel('$P_L$','interpreter','latex');
subplot(3,1,3);
plot(t*WK.f,pow_no(3,:),'r');
hold on;
plot(t*WK.f,pow_ab(3,:),'b');
patch_downstroke(h_pow,t*WK.f,Euler_R_dot);
ylabel('$P_A$','interpreter','latex');
% print(h_pow, 'hover_power_ab', '-depsc');

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
% print(h_E, 'hover_energy_ab', '-depsc');
