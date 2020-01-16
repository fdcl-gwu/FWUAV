%% Evaluates the convergence rate for various insects
load('sim_QS_x_hover_conv_rate.mat');

conv_rates = {conv_rate_osc, conv_rate_hawk, conv_rate_bumb, conv_rate_fruit};
N = size(conv_rates, 2);
N_rand = size(conv_rate_osc, 1);
conv_rate = zeros(N_rand, N);

labels = [{'MONARCH, f=10.2'}, {'HAWKMOTH, f=26.3'}, {'BUMBLEBEE, f=116'}, {'FRUITFLY, f=254'}];

f = figure;
ax = gca;
f.PaperUnits = 'inches';
f.PaperPosition = [0 0 7 6];

for c_ix=4:6
    for i=1:N
    %     sigma = 0.25 * (INSECT{i}.limits(1) - INSECT{i}.limits(2)); % 95 percent confidence
    %     damp_coeff(:, i) = normrnd(INSECT{i}.mu, sigma, [N_rand, 1]);
        conv_rate(:, i) = conv_rates{i}(:, c_ix);
    end
%     boxplot(ax, conv_rate, labels);
    ax.XTick = 1:N;
    scatter(ax, ax.XTick, conv_rate(end, :), 75, 'filled');
    ax.XTickLabel = labels;
    grid on;
    hold on;
end

xlabel('Insect and its flapping frequency (in Hz)');
ylabel('Characteristic exponents of perturbed dynamics');
xlim auto;
ylim auto;
% print('hover_conv_insects', '-depsc', '-r0');