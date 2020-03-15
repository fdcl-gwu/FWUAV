%% Comparison to observe the effect of abdomen oscillation
addpath('../modules', '../sim_data', '../');

%% Monte carlo study
load('sim_QS_x_hover_control_monte_carlo_pos_longitudinal.mat');
err_bound = 5e-4;
% err_bound = squeeze(vecnorm(x_pert, 2, 1))';
x_no_ab = x_pert(:, (abs(err_pos(:, 1)) < err_bound));
bound_ix_no_ab = boundary(x_no_ab(1, :)', x_no_ab(3, :)');
x_with_ab = x_pert(:, (abs(err_pos(:, 2)) < err_bound));
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
legend({'without abdomen effect, $ w = 0 $', 'with abdomen effect, $ w = 0.1 $'},'interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$z$','interpreter','latex');
axis tight;
print(h_err, 'hover_monte_carlo', '-depsc');

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
