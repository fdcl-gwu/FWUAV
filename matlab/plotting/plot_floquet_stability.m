%% Plot sample perturbations and charecteristic solutions
addpath('../modules', '../sim_data', '../');
set(0,'DefaultAxesFontName','times');
set(0,'DefaultAxesFontSize',16);

% time=t*WK.f;
% figure;
% xlabel('$t/T$','interpreter','latex');
% subplot(2, 1, 1);
% plot(time(1:stop_idx), delta_g_mag(1:stop_idx));
% ylabel('$\delta x$','interpreter','latex');
% hold on;
% subplot(2, 1, 2);
% plot(time(1:stop_idx), delta_xi_mag(1:stop_idx));
% ylabel('$\delta \dot{x}$','interpreter','latex');
% print('sim_QS_x_hover_stability', '-depsc');

%%
h_floq = figure;
h_floq.PaperUnits = 'inches';
h_floq.PaperPosition = [0 0 11 11];
for i=1:3
    subplot(3, 2, 2*i-1);
    plot(t*WK.f, squeeze(char_soln_mat(:, i+3, :)));
    lgd = legend(string(1:6));
%     title(lgd, 'Component index');
    l = "$\mathbf{x}_" + string(i+3) + "(t)$" + " with " + "$\mu_" + string(i+3) + " = " + string(mus(i+3)) + "$";
    ylabel(l, 'interpreter', 'latex');
    xlabel('$t/T$', 'interpreter', 'latex');

    j = 2*i;
    subplot(3, 2, j);
    plot(t*WK.f, squeeze(per_val_mat(:, i+3, :)));
    lgd = legend(string(1:6));
    l = "$\mathbf{p}_" + string(i+3) + "(t)$" + " with " + "$\mu_" + string(i+3) + " = " + string(mus(i+3)) + "$";
    ylabel(l, 'interpreter', 'latex');
    xlabel('$t/T$', 'interpreter', 'latex');
end
print(h_floq, 'hover_char_soln', '-depsc', '-r0');
