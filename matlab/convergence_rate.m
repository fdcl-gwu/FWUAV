%% Convergence rate
FRUITFLY.mu = 0.02443;
FRUITFLY.limits = [0.043, 0.002];
BUMBLEBEE.mu = 0.02665;
BUMBLEBEE.limits = [0.041, 0.007];
HAWKMOTH.mu = 0.1632;
HAWKMOTH.limits = [0.31, 0.11];
MONARCH.mu = 0.8558;
MONARCH.limits = [0.92, 0.667];

N = 4;
N_rand = 50;
INSECT = {FRUITFLY, BUMBLEBEE, HAWKMOTH, MONARCH};
damp_coeff = zeros(N_rand, N);

for i=1:N
    sigma = 0.25 * (INSECT{i}.limits(1) - INSECT{i}.limits(2)); % 95 percent confidence
    damp_coeff(:, i) = normrnd(INSECT{i}.mu, sigma, [N_rand, 1]);
end

labels = [{'FRUITFLY, f=254'}; {'BUMBLEBEE, f=116'}; {'HAWKMOTH, f=26.3'}; {'MONARCH, f=10.2'};];
%y = [{'MONARCH'}, {'m=3.1e-4'}; {'FRUITFLY'}, {'m=7.2e-7'}; {'BUMBLEBEE'}, {'m=1.7e-4'}; {'HAWKMOTH'}, {'m=1.6e-3'};];

f = figure;
ax = gca;
f.PaperSize = [9 6.5];
boxplot(ax, damp_coeff, labels);
xlabel('Insect and its flapping frequency (in Hz)');
ylabel('Convergence rate of perturbed dynamics');
% title('');
print('convergence_rate', '-dpdf', '-fillpage');