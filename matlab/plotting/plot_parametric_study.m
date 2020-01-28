%% File to plot the data
addpath('../modules', '../sim_data', '../');
load('parametric_study.mat');

h_f_a = figure;
h_f_a.PaperUnits = 'inches';
h_f_a.PaperPosition = [0 0 12 9];
%
subplot(3,3,1);
plot(eps, squeeze(f_a_m(1,1,:,1)));
hold on;
plot(eps, squeeze(f_a_m(1,1,:,2)), 'r');
legend('Positive component', 'Negative component');
% scatter(eps, squeeze(f_a_m(1,1,:,1)), 10, 'k', 'filled');
ylabel('mean $f_a(1)$','interpreter','latex');

subplot(3,3,4);
plot(eps, squeeze(f_a_m(1,2,:,1)));
hold on;
plot(eps, squeeze(f_a_m(1,2,:,2)), 'r');
ylabel('mean $f_a(2)$','interpreter','latex');

subplot(3,3,7);
% plot(eps, squeeze(f_a_m(1,3,:,1)));
% hold on;
plot(eps, squeeze(f_a_m(1,3,:,2)), 'r');
xlabel('$\epsilon\ \vert\ \Delta\phi_{m, R} = \epsilon, \Delta\phi_{m, L} = \epsilon$','interpreter','latex');
ylabel('mean $f_a(3)$','interpreter','latex');
%
subplot(3,3,2);
plot(eps, squeeze(f_a_m(2,1,:,1)));
hold on;
plot(eps, squeeze(f_a_m(2,1,:,2)), 'r');
legend('Positive component', 'Negative component');

subplot(3,3,5);
plot(eps, squeeze(f_a_m(2,2,:,1)));
hold on;
plot(eps, squeeze(f_a_m(2,2,:,2)), 'r');

subplot(3,3,8);
plot(eps, squeeze(f_a_m(2,3,:,2)), 'r');
xlabel('$\epsilon\ \vert\ \Delta\theta_{m, R} = \epsilon, \Delta\theta_{m, L} = \epsilon$','interpreter','latex');
%
subplot(3,3,3);
plot(eps, squeeze(f_a_m(3,1,:,1)));
hold on;
plot(eps, squeeze(f_a_m(3,1,:,2)), 'r');
legend('Positive component', 'Negative component');

subplot(3,3,6);
plot(eps, squeeze(f_a_m(3,2,:,1)));
hold on;
plot(eps, squeeze(f_a_m(3,2,:,2)), 'r');

subplot(3,3,9);
plot(eps, squeeze(f_a_m(3,3,:,2)), 'r');
% xlabel('$\epsilon\ \vert\ \Delta\psi_{m, R} = \epsilon, \Delta\psi_{m, L} = -\epsilon$','interpreter','latex');
xlabel('$\epsilon\ \vert\ \Delta\phi_{m, R} = \epsilon, \Delta\phi_{m, L} = -\epsilon$','interpreter','latex');
%
% print(h_f_a, 'hover_param_study', '-depsc', '-r0');
