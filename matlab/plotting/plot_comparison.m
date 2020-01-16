%% Comparison of power and energy to observe the effect of abdomen oscillation
addpath('../modules', '../sim_data', '../');

% Without abdomen oscillation
load('sim_QS_x_hover_ab_no.mat')
[pow_no, E_no, E_dot_no, eff_no] = compute_power(INSECT.m, t, x, x_dot, tau, Q_R, Q_L, Q_A, W_R, W_L, W_A);

% With abdomen oscillation
load('sim_QS_x_hover.mat')
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
