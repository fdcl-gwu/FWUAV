%% File to plot the data
addpath('../modules', '../sim_data', '../');
set(0,'DefaultAxesFontName','times');
set(0,'DefaultAxesFontSize',16);
load('parametric_study.mat');

h_f_a = figure;
h_f_a.PaperUnits = 'inches';
h_f_a.PaperPosition = [0 0 13 9];
nr = 3;
nc = 4;
%
ic = 1;

ir = 1;
subplot(nr,nc,ic);
plot(eps, squeeze(f_a_m(ic,ir,:,1)));
ylabel('mean $f_a(1)$','interpreter','latex');
hold on;
yyaxis right;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
% legend('Positive component', 'Negative component');
title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e\n', params.df_a_1_by_dphi_m), 'interpreter','latex');

ir = 2;
subplot(nr,nc,ic+nc);
plot(eps, squeeze(f_a_m(ic,ir,:,1)));
ylabel('mean $f_a(2)$','interpreter','latex');
hold on;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');

ir = 3;
subplot(nr,nc,ic+2*nc);
% plot(eps, squeeze(f_a_m(ic,id,:,1)));
% hold on;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
xlabel(sprintf('$\\epsilon\\ \\vert\\ \\Delta\\phi_{m, R} = \\epsilon$,\n $\\Delta\\phi_{m, L} = \\epsilon$'),...
    'interpreter','latex');
ylabel('mean $f_a(3)$','interpreter','latex');
title(sprintf('$m_n$ = %0.3e\n', params.df_a_3_by_dphi_m),'interpreter','latex');
%
ic = 2;

ir = 1;
subplot(nr,nc,ic);
plot(eps, squeeze(f_a_m(ic,ir,:,1)));
hold on;
yyaxis right;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e\n', params.df_a_1_by_dtheta_0),'interpreter','latex');

ir = 2;
subplot(nr,nc,ic+nc);
plot(eps, squeeze(f_a_m(ic,ir,:,1)));
hold on;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');

ir = 3;
subplot(nr,nc,ic+2*nc);
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
xlabel(sprintf('$\\epsilon\\ \\vert\\ \\Delta\\theta_{0, R} = \\epsilon$,\n $\\Delta\\theta_{0, L} = \\epsilon$'),...
    'interpreter','latex');
title(sprintf('$m_n$ = %0.3e\n', params.df_a_3_by_dtheta_0),'interpreter','latex');
%
ic = 3;

ir = 1;
subplot(nr,nc,ic);
plot(eps, squeeze(f_a_m(ic,ir,:,1)));
hold on;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');

ir = 2;
subplot(nr,nc,ic+nc);
plot(eps, squeeze(f_a_m(ic,ir,:,1)));
hold on;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e\n', params.df_a_2_by_dphi_m),'interpreter','latex');

ir = 3;
subplot(nr,nc,ic+2*nc);
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
xlabel(sprintf('$\\epsilon\\ \\vert\\ \\Delta\\phi_{m, R} = \\epsilon$,\n $\\Delta\\phi_{m, L} = -\\epsilon$'),...
    'interpreter','latex');
%
ic = 4;

ir = 1;
subplot(nr,nc,ic);
plot(eps, squeeze(f_a_m(ic,ir,:,1)));
hold on;
yyaxis right;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e\n', params.df_a_1_by_dtheta_A_m),'interpreter','latex');

ir = 2;
subplot(nr,nc,ic+nc);
plot(eps, squeeze(f_a_m(ic,ir,:,1)));
hold on;
yyaxis right;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');

ir = 3;
subplot(nr,nc,ic+2*nc);
plot(eps, squeeze(f_a_m(ic,ir,:,1)));
hold on;
yyaxis right;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
xlabel('$\epsilon\ \vert\ \Delta\theta_{A_m} = \epsilon$','interpreter','latex');
title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e\n', params.df_a_3_by_dtheta_A_m),'interpreter','latex');
%
print(h_f_a, 'hover_param_study', '-depsc', '-r0');
