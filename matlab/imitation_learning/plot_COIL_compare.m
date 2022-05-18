addpath('../modules', '../sim_data', '../');
set(0,'DefaultAxesFontName','times');
set(0,'DefaultAxesFontSize',16);

cloning = load('control_performance_cloning');
dagger = load('control_performance_dagger');
coil = load('control_performance_coil');

idx_plot = 1:201;
[~, cloning.t_decay] = min(sum(cloning.cost(cloning.cost(:, 1) < 1, idx_plot), 1)); % no ultimate bound
dagger.t_decay = sum(~all(dagger.cost(dagger.cost(:, 1) < 1, idx_plot) < dagger.ult_bound, 1)) + 1;
coil.t_decay = sum(~all(coil.cost(coil.cost(:, 1) < 1, idx_plot) < coil.ult_bound, 1)) + 1;

x0 = [1, 1e-2]; % initial exp point, decay rate
datas = {cloning, dagger, coil};
for i =1:3
    start_idx = 1:datas{i}.t_decay;
    cost = datas{i}.cost(datas{i}.cost(:, 1) < 1, idx_plot);
    bound_area = @(x) sum((log(x(1)*exp(-x(2) * (start_idx-1))) - log(max(cost(:, start_idx), [], 1))) ./ start_idx);
    nonlcon = @(x) deal(-(x(1)*exp(-x(2) * (start_idx-1)) - max(cost(:, start_idx), [], 1)), 0);
    x = fmincon(bound_area, x0, [], [], [], [], [1, 0], [], nonlcon, optimoptions('fmincon','ObjectiveLimit',0));
    datas{i}.decay_e0 = x(1); datas{i}.decay_rate = x(2);
end
[cloning, dagger, coil] = datas{:};

h_comp = figure;
h_comp.PaperUnits = 'inches';
h_comp.PaperPosition = [0 0 7 6];
hold on;
ax = gca;
loglog(ax, idx_plot, cloning.cost(cloning.cost(:, 1) < 1, idx_plot)', 'LineStyle', '-', 'LineWidth', 0.5, 'Color', [0.9290, 0.6940, 0.1250, 0.1]);
loglog(ax, idx_plot, dagger.cost(dagger.cost(:, 1) < 1, idx_plot)', 'LineStyle', '-', 'LineWidth', 0.5, 'Color', [0.8500, 0.3250, 0.0980, 0.1]);
loglog(ax, idx_plot, coil.cost(coil.cost(:, 1) < 1, idx_plot)', 'LineStyle', '-', 'LineWidth', 0.5, 'Color', [0, 0.4470, 0.7410, 0.1]);
b_ult = loglog(ax, idx_plot, dagger.ult_bound * ones(size(idx_plot)), '--', 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980]);
c_ult = loglog(ax, idx_plot, coil.ult_bound * ones(size(idx_plot)), '--', 'LineWidth', 1, 'Color', [0, 0.4470, 0.7410]);
loglog(ax, 1:dagger.t_decay, dagger.decay_e0*exp(-dagger.decay_rate * ((1:dagger.t_decay)-1)), '--', 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980]);
loglog(ax, 1:coil.t_decay, coil.decay_e0*exp(-coil.decay_rate * ((1:coil.t_decay)-1)), '--', 'LineWidth', 1, 'Color', [0, 0.4470, 0.7410]);
ax.XScale = 'log';
ax.YScale = 'log';
xlabel('$t/T$ (log scale)','interpreter','latex');
ylabel('Weighted state error (log scale)');
legend([b_ult, c_ult], {'DAGGER bound', 'COIL bound'});
% print(h_comp, 'state_error_comp', '-depsc');
%%% Too heavy to export complete vector image
% exportgraphics(h_comp,'state_error_comp.pdf','ContentType','vector');

% exportgraphics(h_comp,'state_error_comp.pdf','Resolution', 300);
