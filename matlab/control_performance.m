evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
filename = 'control_performance';

load('sim_QS_xR_hover_control_opt_200_WW', 't', 'X_ref0', ...
    'WK_R', 'WK_L', 'INSECT', 'N_single', 'N_iters', 'N_per_iter', ...
    'get_args', 'param_type', 'm', 'N_dang', 'Weights', 'N_sims', ...
    'des', 'opt_complete');

load('iterative_learning_result', 'control_net');
WK = WK_R;

%%
rng default;
N_periods = 200; % 1, 3
N_iters = 10; % 2, 4, 10 per period
m = N_iters; % Control Horizon multiplier of iter
N_single = 100; % per time period
N_per_iter = N_single / N_iters;
N = N_single*N_periods + 1;
T = N_periods /WK.f;
t = linspace(0,T,N);
dt = t(2)-t(1);
X_ref0 = des.X0;

%%
rng(1);
N_sims = 1500;
scale = logspace(1, -2, 13); % 0.5/-2/11, 029, 037, 045
% scale = logspace(0, -0.001, 11); % 029, 037, 045
% scale = linspace(1, 0.1, 10);
N_scale = length(scale);
dX = zeros(N_sims*N_scale, 12);
dX_scale = zeros(N_sims*N_scale, 1);

for i = 1:N_sims
    dx = 2*rand(1,3)-1; dx = rand(1) * dx / norm(dx);
    dtheta = 2*rand(1,3)-1; dtheta = rand(1) * dtheta / norm(dtheta);
    dx_dot = 2*rand(1,3)-1; dx_dot = rand(1) * dx_dot / norm(dx_dot);
    domega = 2*rand(1,3)-1; domega = rand(1) * domega / norm(domega);
    for j = 1:N_scale
        dX(i+(j-1)*N_sims,:) = scale(j) * Weights.PerturbVariables .* [dx, dtheta, dx_dot, domega];
        dX_scale(i+(j-1)*N_sims) = scale(j) * norm([dx, dtheta, dx_dot, domega]);
    end
end

%%
tic;
cost = zeros(N_sims*N_scale, N_periods+1);
% X_T = zeros(N_sims*N_scale, N_periods+1, 18);

parfor i=1:N_sims*N_scale

    X = zeros(1+N_single*N_periods, 18);
    dang0 = zeros(N_dang, 1);
    cost_arr = zeros(1+N_periods, 1);

    dX0 = dX(i, :)';
    X0 = [X_ref0(1:3)+ dX0(1:3); reshape(reshape(X_ref0(4:12),3,3)*expmhat(dX0(4:6)), 9, 1); ...
        X_ref0(13:18) + dX0(7:12)];
    cost_arr(1) = sqrt(sum((Weights.OutputVariables' .* (X0 - X_ref0)).^2));

    for period=1:N_periods
        idx = (1+(period-1)*N_single):(1+period*N_single);
        param = control_net((1 ./ Weights.PerturbVariables)' .* get_error(X0, X_ref0));

        X(idx, :) = crgr_xR_control_mex(INSECT, WK_R, WK_L, t(idx), X0, dang0, param, N_dang, m, N_per_iter);
        X0 = X(idx(end),:)';
        dang0 = dang0 + sum(reshape(param, 6, 10), 2) / WK_R.f / m;
        cost_arr(period+1) = sqrt(sum((Weights.OutputVariables .* (X(idx(end), :) - X_ref0')).^2));
    end
    cost(i, :) = cost_arr;
	% X_T(i, :, :) = X(1:N_single:end, :);

end

time_taken = toc;

%% Change of weights
% WW = Weights.OutputVariables;
% % WW = ones(1, 18);
% WW = (1 ./ Weights.PerturbVariables); WW = [WW(1:3) WW(4:6) WW(4:6) WW(4:6) WW(7:12)];
% for i=1:N_scale*N_sims
% 	    cost(i, :) = sqrt(sum((WW .* (squeeze(X_T(i, :, :)) - X_ref0')).^2, 2));
% end

%% Performance measures
idx_cost_good = cost(:, 1) < 1;
log_cost = log(cost);
% d_log = (log_cost(:, 1) - log_cost(:, 2:end)) ./ (1:N_periods);
d_log = log_cost(:, 1:end-1) - log_cost(:, 2:end);

ult_bound = max(cost(idx_cost_good, 10:end), [], 'all');
decay_rate = -Inf;
while decay_rate < 0
    ult_bound = ult_bound * 1.01;
    idx_ult_bound = (cost(:, 1:end-1) > ult_bound) & (cost(:, 2:end) > ult_bound) & idx_cost_good;
    decay_rate = min(d_log(idx_ult_bound));
end

% x_fit = start_idx';
% conv_rate = zeros(N_sims*N_scale, 1);
% for i=1:N_sims*N_scale
%     if cost(i, 1) < 1
%         y_fit = cost(i, start_idx)';
%         f = fit(x_fit-1, y_fit, 'exp1');
%         conv_rate(i) = f.b;
%     end
% end
% figure;
% plot(conv_rate(cost(:, 1) > 0.01 & cost(:, 1) < 1));
% decay_rate = max(conv_rate(cost(:, 1) < 1 & cost(:, 1) > ult_bound));

%%
% all(diff(max(cost(cost(:, 1) < 1, 10:end), [], 1)) < 0)
% a = diff(mean(cost(cost(:, 1) < 1, 10:end), 1));
% sum(a(a>0)), sum(a(a<0))
% 
% start_idx = 1:10;
% decay_rate = 1; decay_factor = 10;
% for k=1:10
%     conv_rates = [decay_rate/decay_factor, decay_rate, decay_rate * decay_factor];
%     c_r_bool = ones(size(conv_rates), 'logical');
%     for j = 1:length(conv_rates)
%         c_r = conv_rates(j);
%         for i=1:N_sims*N_scale
%             cost_idx = (cost(i, start_idx) >=  ult_bound) & (cost(i, 1) * exp(-c_r * (start_idx - 1)) >= ult_bound);
%             if (cost(i, 1) < 1) && (cost(i, 1) > ult_bound) ...
%                     && any(cost(i, start_idx(cost_idx)) > cost(i, 1) * exp(-c_r * (start_idx(cost_idx) - 1)))
%                 c_r_bool(j) = false;
%                 break;
%             end
%         end
%         if ~c_r_bool(j)
%             try
%                 c_r_bool((j+1):end) = false;
%             end
%             break;
%         end
%     end
%     decay_rate = conv_rates(c_r_bool);
%     decay_rate = decay_rate(end);
%     decay_factor = sqrt(decay_factor);
% end

allvars = whos;
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
save(filename, allvars(tosave).name)
evalin('base',['load ' filename]);

%%
function err = get_error(X, Xd)
err = zeros(12, 1);
R = reshape(X(4:12), 3, 3); Rd = reshape(Xd(4:12), 3, 3);
err(1:3) = X(1:3) - Xd(1:3);
err(4:6) = 0.5*vee(Rd'*R - R'*Rd);
% [err(4), err(6), err(5)] = dcm2angle(R'*Rd, 'xzy');
err(7:9) = X(13:15) - Xd(13:15);
err(10:12) = X(16:18) - (R'*Rd*Xd(16:18));
end

%% Initial studies
% L_scale = zeros(N_sims*(N_scale -1), N_periods+1);
% for j=1:N_scale-1
%     L_scale((1+(j-1)*N_sims):(j*N_sims), :) = ...
%         (cost((1+(j-1)*N_sims):(j*N_sims), :) - cost((1+(j)*N_sims):((j+1)*N_sims), :)) ...
%         / (scale(j) - scale(j+1));
% end

% cutoff_cost = 0.01;
% idx_scale = ones(N_sims*(N_scale -1), N_periods+1, 'logical');
% idx_scale((cost(1:N_sims*(N_scale-1), :) < cutoff_cost)) = false;
% % Initial scales not reliable since OutputVariables & PerturbVariables don't match
% idx_scale((cost(1:N_sims*(N_scale-1), 1) < cost((1+N_sims):N_sims*N_scale, 1)), :) = false;
% disp('Ratio of elements not cutoff because of cost');
% disp(sum(idx_scale, 'all') / numel(idx_scale));
% disp('Number of elements for which cost increases when scale decreases');
% disp(sum(L_scale(idx_scale) < 0));
% disp('Ratio of such elements');
% disp(sum(L_scale(idx_scale) < 0) / numel(L_scale(idx_scale)));
% idx_bad = (L_scale < 0) & (idx_scale == true);
% disp(min(L_scale(idx_bad)));
% [row_idx, col_idx] = find(idx_bad);

% int_L = cumsum(L_scale(:, 2:end), 2) * 1/WK.f;
% any((int_L <= 0) & (idx_scale(:, 2:end)), 'all')

% gap_period = 1; % 1 to N_periods
% L_time = diff(cost(:, 1:gap_period:(N_periods+1))')' * WK.f / gap_period;

% cutoff_cost = 0.01;
% idx_cost = ones(N_sims*N_scale, N_periods/gap_period, 'logical');
% idx_cost((cost(:, 1:gap_period:N_periods) < cutoff_cost)) = false;
% disp('Ratio of elements not cutoff because of cost');
% disp(sum(idx_cost, 'all') / numel(idx_cost));
% disp('Number of elements for which cost increases when time increases');
% disp(sum(L_time(idx_cost) > 0));
% disp('Ratio of such elements');
% disp(sum(L_time(idx_cost) > 0) / numel(L_time(idx_cost)));
% idx_bad = (L_time > 0) & (idx_cost == true);
% find(any(idx_bad, 2) == 1);

% figure;
% plot(L_time');
% xlabel('time period');
% ylabel('$\frac{\Delta cost}{\Delta t}$', 'Interpreter', 'latex');
