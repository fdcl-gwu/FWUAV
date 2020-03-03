%% File to plot the data
addpath('../modules', '../sim_data', '../');
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
title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e', get_slope(eps', squeeze(f_a_m(ic,ir,:,1))), ...
    get_slope(eps', squeeze(f_a_m(ic,ir,:,2)))),'interpreter','latex');

ir = 2;
subplot(nr,nc,ic+nc);
plot(eps, squeeze(f_a_m(ic,ir,:,1)));
ylabel('mean $f_a(2)$','interpreter','latex');
hold on;
yyaxis right;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e', get_slope(eps', squeeze(f_a_m(ic,ir,:,1))), ...
    get_slope(eps', squeeze(f_a_m(ic,ir,:,2)))),'interpreter','latex');

ir = 3;
subplot(nr,nc,ic+2*nc);
% plot(eps, squeeze(f_a_m(ic,id,:,1)));
% hold on;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
xlabel('$\epsilon\ \vert\ \Delta\phi_{m, R} = \epsilon, \Delta\phi_{m, L} = \epsilon$','interpreter','latex');
ylabel('mean $f_a(3)$','interpreter','latex');
title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e', get_slope(eps', squeeze(f_a_m(ic,ir,:,1))), ...
    get_slope(eps', squeeze(f_a_m(ic,ir,:,2)))),'interpreter','latex');
%
ic = 2;

ir = 1;
subplot(nr,nc,ic);
plot(eps, squeeze(f_a_m(ic,ir,:,1)));
hold on;
yyaxis right;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e', get_slope(eps', squeeze(f_a_m(ic,ir,:,1))), ...
    get_slope(eps', squeeze(f_a_m(ic,ir,:,2)))),'interpreter','latex');

ir = 2;
subplot(nr,nc,ic+nc);
plot(eps, squeeze(f_a_m(ic,ir,:,1)));
hold on;
yyaxis right;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e', get_slope(eps', squeeze(f_a_m(ic,ir,:,1))), ...
    get_slope(eps', squeeze(f_a_m(ic,ir,:,2)))),'interpreter','latex');

ir = 3;
subplot(nr,nc,ic+2*nc);
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
xlabel('$\epsilon\ \vert\ \Delta\theta_{0, R} = \epsilon, \Delta\theta_{m0, L} = \epsilon$','interpreter','latex');
title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e', get_slope(eps', squeeze(f_a_m(ic,ir,:,1))), ...
    get_slope(eps', squeeze(f_a_m(ic,ir,:,2)))),'interpreter','latex');
%
ic = 3;

ir = 1;
subplot(nr,nc,ic);
plot(eps, squeeze(f_a_m(ic,ir,:,1)));
hold on;
yyaxis right;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e', get_slope(eps', squeeze(f_a_m(ic,ir,:,1))), ...
    get_slope(eps', squeeze(f_a_m(ic,ir,:,2)))),'interpreter','latex');

ir = 2;
subplot(nr,nc,ic+nc);
plot(eps, squeeze(f_a_m(ic,ir,:,1)));
hold on;
yyaxis right;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e', get_slope(eps(eps'>0)', squeeze(f_a_m(ic,ir,eps'>0,1))), ...
    get_slope(eps(eps'<0)', squeeze(f_a_m(ic,ir,eps'<0,1)))),'interpreter','latex');

ir = 3;
subplot(nr,nc,ic+2*nc);
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
xlabel('$\epsilon\ \vert\ \Delta\phi_{m, R} = \epsilon, \Delta\phi_{m, L} = -\epsilon$','interpreter','latex');
title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e', get_slope(eps', squeeze(f_a_m(ic,ir,:,1))), ...
    get_slope(eps', squeeze(f_a_m(ic,ir,:,2)))),'interpreter','latex');
%
ic = 4;

ir = 1;
subplot(nr,nc,ic);
plot(eps, squeeze(f_a_m(ic,ir,:,1)));
hold on;
yyaxis right;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e', get_slope(eps', squeeze(f_a_m(ic,ir,:,1))), ...
    get_slope(eps', squeeze(f_a_m(ic,ir,:,2)))),'interpreter','latex');

ir = 2;
subplot(nr,nc,ic+nc);
plot(eps, squeeze(f_a_m(ic,ir,:,1)));
hold on;
yyaxis right;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e', get_slope(eps', squeeze(f_a_m(ic,ir,:,1))), ...
    get_slope(eps', squeeze(f_a_m(ic,ir,:,2)))),'interpreter','latex');

ir = 3;
subplot(nr,nc,ic+2*nc);
plot(eps, squeeze(f_a_m(ic,ir,:,1)));
hold on;
yyaxis right;
plot(eps, squeeze(f_a_m(ic,ir,:,2)), 'r');
xlabel('$\epsilon\ \vert\ \Delta\theta_{A_m} = \epsilon$','interpreter','latex');
title(sprintf('$m_p$ = %0.3e \n $m_n$ = %0.3e', get_slope(eps', squeeze(f_a_m(ic,ir,:,1))), ...
    get_slope(eps', squeeze(f_a_m(ic,ir,:,2)))),'interpreter','latex');
%
print(h_f_a, 'hover_param_study', '-depsc', '-r0');

function m = get_slope(x, y)
    X = ones(length(x), 2);
    if size(X, 1) == size(y, 1)
        X(:, 2) = x;
        m = X \ y;
        m = m(2);
    end
end
