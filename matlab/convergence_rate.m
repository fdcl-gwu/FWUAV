%% Convergence rate
load('sim_QS_x_hover_conv_rate.mat');

conv_rates = {conv_rate_fruit, conv_rate_bumb, conv_rate_hawk, conv_rate_osc};
N = size(conv_rates, 2);
N_rand = size(conv_rate_osc, 1);
conv_rate = zeros(N_rand, N);

labels = [{'FRUITFLY, f=254'}; {'BUMBLEBEE, f=116'}; {'HAWKMOTH, f=26.3'}; {'MONARCH, f=10.2'};];
%y = [{'MONARCH'}, {'m=3.1e-4'}; {'FRUITFLY'}, {'m=7.2e-7'}; {'BUMBLEBEE'}, {'m=1.7e-4'}; {'HAWKMOTH'}, {'m=1.6e-3'};];

f = figure;
ax = gca;
f.PaperSize = [9 6.5];

for c_ix=4:6
    for i=1:N
    %     sigma = 0.25 * (INSECT{i}.limits(1) - INSECT{i}.limits(2)); % 95 percent confidence
    %     damp_coeff(:, i) = normrnd(INSECT{i}.mu, sigma, [N_rand, 1]);
        conv_rate(:, i) = conv_rates{i}(:, c_ix);
    end
    boxplot(ax, conv_rate, labels);
    hold on;
end

xlabel('Insect and its flapping frequency (in Hz)');
ylabel('Convergence rate of perturbed dynamics');
xlim auto;
ylim auto;
% title('');
print('convergence_rate', '-dpdf', '-fillpage');