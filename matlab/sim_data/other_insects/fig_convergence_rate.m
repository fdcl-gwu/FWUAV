%% Evaluates the convergence rate for various insects
load('sim_QS_x_hover_conv_rate.mat');
set(0,'DefaultAxesFontName','times');
set(0,'DefaultAxesFontSize',14);

conv_rates = [conv_rate_mona, conv_rate_hawk, conv_rate_bumb, conv_rate_fruit];
N = size(conv_rates, 2);

labels = [{'MONARCH, f=10.2'}, {'HAWKMOTH, f=26.3'}, {'BUMBLEBEE, f=116'}, {'FRUITFLY, f=254'}];

f = figure;
ax = gca;
f.PaperUnits = 'inches';
f.PaperPosition = [0 0 8 6];

for c_ix=4:6
%     boxplot(ax, conv_rate, labels);
    ax.XTick = 1:N;
    scatter(ax, ax.XTick, conv_rates(c_ix, :), 75, 'filled');
    ax.XTickLabel = labels;
    grid on;
    hold on;
end

xlim auto;
ylim auto;
xlabel('Insect and its flapping frequency (in Hz)');
ylabel('Characteristic exponents of perturbed dynamics');
print('hover_conv_insects', '-depsc', '-r0');
