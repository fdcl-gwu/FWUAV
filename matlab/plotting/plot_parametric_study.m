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
subplot(nr,nc,ic);
plot(eps, squeeze(f_a_m(ic,1,:,1)));
hold on;
plot(eps, squeeze(f_a_m(ic,1,:,2)), 'r');
legend('Positive component', 'Negative component');
ylabel('mean $f_a(1)$','interpreter','latex');

subplot(nr,nc,ic+nc);
plot(eps, squeeze(f_a_m(ic,2,:,1)));
hold on;
plot(eps, squeeze(f_a_m(ic,2,:,2)), 'r');
ylabel('mean $f_a(2)$','interpreter','latex');

subplot(nr,nc,ic+2*nc);
plot(eps, squeeze(f_a_m(ic,3,:,2)), 'r');
xlabel('$\epsilon\ \vert\ \Delta\phi_{m, R} = \epsilon, \Delta\phi_{m, L} = \epsilon$','interpreter','latex');
ylabel('mean $f_a(3)$','interpreter','latex');
%
ic = 2;
subplot(nr,nc,ic);
plot(eps, squeeze(f_a_m(ic,1,:,1)));
hold on;
plot(eps, squeeze(f_a_m(ic,1,:,2)), 'r');
legend('Positive component', 'Negative component');

subplot(nr,nc,ic+nc);
plot(eps, squeeze(f_a_m(ic,2,:,1)));
hold on;
plot(eps, squeeze(f_a_m(ic,2,:,2)), 'r');

subplot(nr,nc,ic+2*nc);
plot(eps, squeeze(f_a_m(ic,3,:,2)), 'r');
xlabel('$\epsilon\ \vert\ \Delta\theta_{m, R} = \epsilon, \Delta\theta_{m, L} = \epsilon$','interpreter','latex');
%
ic = 3;
subplot(nr,nc,ic);
plot(eps, squeeze(f_a_m(ic,1,:,1)));
hold on;
plot(eps, squeeze(f_a_m(ic,1,:,2)), 'r');
legend('Positive component', 'Negative component');

subplot(nr,nc,ic+nc);
plot(eps, squeeze(f_a_m(ic,2,:,1)));
hold on;
plot(eps, squeeze(f_a_m(ic,2,:,2)), 'r');

subplot(nr,nc,ic+2*nc);
plot(eps, squeeze(f_a_m(ic,3,:,2)), 'r');
xlabel('$\epsilon\ \vert\ \Delta\phi_{m, R} = \epsilon, \Delta\phi_{m, L} = -\epsilon$','interpreter','latex');
%
ic = 4;
subplot(nr,nc,ic);
plot(eps, squeeze(f_a_m(ic,1,:,1)));
hold on;
plot(eps, squeeze(f_a_m(ic,1,:,2)), 'r');
legend('Positive component', 'Negative component');

subplot(nr,nc,ic+nc);
plot(eps, squeeze(f_a_m(ic,2,:,1)));
hold on;
plot(eps, squeeze(f_a_m(ic,2,:,2)), 'r');

subplot(nr,nc,ic+2*nc);
plot(eps, squeeze(f_a_m(ic,3,:,1)));
hold on;
plot(eps, squeeze(f_a_m(ic,3,:,2)), 'r');
xlabel('$\epsilon\ \vert\ \Delta\theta_{A_m} = \epsilon$','interpreter','latex');
%
print(h_f_a, 'hover_param_study', '-depsc', '-r0');
