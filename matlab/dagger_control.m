evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
filename = 'iterative_learning';

load('sim_QS_xR_hover_control_opt_200_WW', 't', 'X_ref0', ...
    'WK_R', 'WK_L', 'INSECT', 'N_single', 'N_iters', 'N_per_iter', ...
    'get_args', 'param_type', 'm', 'N_dang', 'Weights', 'N_sims', ...
    'X0_pert', 'opt_param', 'des', 'opt_complete', 'cost', 'problem', 'solver', 'X_ref');

X0_pert_flat = squeeze(X0_pert(:, 1, :));
opt_param_flat = reshape(permute(opt_param, [1, 3, 2]), 3*N_sims, 6*N_iters);
err = zeros(3*N_sims, 12);
for i=1:(3*N_sims)
    err(i, :) = get_error(X0_pert_flat(i, :)', des.X0)';
end
inputs = (1 ./ Weights.PerturbVariables)' .* err'; targets = opt_param_flat';
opt_complete = opt_complete & (cost(:, end) < cost(:, 1));
% inputs = inputs(:, opt_complete); targets = targets(:, opt_complete);
idx = 1:1000;
inputs = inputs(:, opt_complete([idx, N_sims+idx, 2*N_sims+idx]));
targets = targets(:, opt_complete([idx, N_sims+idx, 2*N_sims+idx]));

% N_zero = 1000; % 25, 100
N_data = size(inputs, 2);
N_zero = round(N_data/15);
inputs(:, (N_data+1):(N_data+N_zero)) = zeros(12, N_zero);
targets(:, (N_data+1):(N_data+N_zero)) = zeros(60, N_zero);

N_features = size(inputs, 1);
[N_outputs, N_data] = size(targets);
N_dagger_iters = 5;

%% NN Model
control_net = cascadeforwardnet([60]); % 36, 60
control_net = configure(control_net, inputs, targets);
trainFcn = 'trainbr'; % trainbr, trainrp
control_net.trainFcn = trainFcn; % trainbr, trainrp
switch trainFcn
    case 'trainbr'
        control_net.trainParam.epochs = 25;
    case 'trainrp'
        control_net.trainParam.epochs = 20000;
        control_net.trainParam.max_fail = 20;
        control_net.performParam.regularization = 0.1; % if not trainbr
end
control_net.performParam.normalization = 'standard';
control_net.inputs{1}.processFcns = {'mapminmax', 'processpca'};
control_net.layers{1:(end-1)}.transferFcn = 'leakyrelu'; % leakyrelu, tansig, poslin, purelin
control_net.divideParam.trainRatio = 80/100;
control_net.divideParam.valRatio = 10/100;
control_net.divideParam.testRatio = 10/100;

%%
rng default;
WK = WK_R;
N_periods = 5; % 1, 3
N_iters = 10; % 2, 4, 10 per period
m = N_iters; % Control Horizon multiplier of iter
p = 2*m;
N_single = 100; % per time period
N_per_iter = N_single / N_iters;
N = N_single*N_periods + 1;
T = N_periods /WK.f;
t = linspace(0,T,N);
dt = t(2)-t(1);
X_ref0 = des.X0;

%%
rng(1);
N_sims = 10;
scale = logspace(0, -2, 3);
N_scale = length(scale);
dX = zeros(N_sims*N_scale, 12);
dX_scale = zeros(N_sims*N_scale, 1);

% for i = 1:N_sims
%     dx = 2*rand(1,3)-1; dx = rand(1) * dx / norm(dx);
%     dtheta = 2*rand(1,3)-1; dtheta = rand(1) * dtheta / norm(dtheta);
%     dx_dot = 2*rand(1,3)-1; dx_dot = rand(1) * dx_dot / norm(dx_dot);
%     domega = 2*rand(1,3)-1; domega = rand(1) * domega / norm(domega);
%     for j = 1:N_scale
%         dX(i+(j-1)*N_sims,:) = scale(j) * Weights.PerturbVariables .* [dx, dtheta, dx_dot, domega];
%         dX_scale(i+(j-1)*N_sims) = scale(j) * norm([dx, dtheta, dx_dot, domega]);
%     end
% end

%% Algorithm
y_star = targets;
perf = zeros(N_dagger_iters, 1); perf_z = zeros(N_dagger_iters, 1);
cons = zeros(N_dagger_iters, 1); cons_z = zeros(N_dagger_iters, 1);
cost = zeros(N_sims*N_scale, N_periods+1);
X_T = zeros(N_sims*N_scale, N_periods+1, 18);

tic;
for iter=1:N_dagger_iters
    %% TRAINING
    control_net = init(control_net);
	% control_net = configure(control_net, inputs, targets);
    control_net = train(control_net, inputs, targets, 'useParallel', 'yes'); % 'useParallel', 'yes'

	%% Rollout
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
		X_T(i, :, :) = X(1:N_single:end, :);
	end

	%% Expert labels
	X0_pi = reshape(X_T(:, 2:end, :), [], 18);
	N_pi = size(X0_pi, 1);
	is_bounded = zeros(N_pi, 1, 'logical');
	opt_complete_pi = zeros(N_pi, 1, 'logical');
	inputs_pi = zeros(N_features, N_pi);
	targets_pi = zeros(N_outputs, N_pi);
	N_period_pi = 1;
	parfor i=1:N_pi
		X0 = X0_pi(i, :)';
		inputs_pi(:, i) = (1 ./ Weights.PerturbVariables)' .* get_error(X0, X_ref0);
		if norm(inputs_pi(:, i)) <= 1
			is_bounded(i) = true;
			[cost, opt_param, opt_fval, opt_iter, opt_firstorder, opt_complete] = ...
				simulate_opt_control(t, X0, X_ref, WK_R, WK_L, INSECT, ...
				N, N_period_pi, N_single, N_iters, N_per_iter, problem, solver, Weights, ...
				get_args, param_type, p, m, 0);
			opt_complete_pi(i) = opt_complete;
			targets_pi(:, i) = opt_param;
		end
	end

	inputs = [inputs, inputs_pi(:, opt_complete_pi & is_bounded)];
	targets = [targets, targets_pi(:, opt_complete_pi & is_bounded)];

    %% Performance
    y = control_net(inputs);
    cons(iter) = norm(control_net(zeros(N_features,1)));
    perf(iter) = perform(control_net, targets, y);
end

time_taken = toc;

%%
allvars = whos;
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
save(filename, allvars(tosave).name)
evalin('base',['load ' filename]);

%%
function [cost, opt_param, opt_fval, opt_iter, opt_firstorder, opt_complete, dang, ...
    X, x, x_dot, R, Q_R, Q_L, Q_A, theta_A, W, W_R, W_R_dot, W_L, ...
    W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, f_tau, tau, Euler_R, ...
    Euler_R_dot] = simulate_opt_control(t, X0, X_ref, WK_R, WK_L, INSECT, ...
    N, N_periods, N_single, N_iters, N_per_iter, problem, solver, Weights, ...
    get_args, param_type, p, m, control_net)
%%
cost = zeros(1, 1+N_periods*N_iters);
cost(1) = sqrt(sum((Weights.OutputVariables .* (X0' - X_ref(1, :))).^2));
X = zeros(1+N_single*N_periods, 18);
des_X0 = X_ref(1, :)';
N_p = length(problem.x0) / p;
dang = zeros(N_p, 1+N_single*N_periods);
opt_param = zeros(N_p*m, N_periods*N_iters/m); opt_fval = zeros(1, N_periods*N_iters/m);
opt_iter = zeros(1, N_periods*N_iters/m); opt_firstorder = zeros(1, N_periods*N_iters/m);
opt_complete = zeros(1, N_periods*N_iters/m, 'logical');

if isa(control_net, 'network')
    X_temp0 = X0;
    param = sim(control_net, (1 ./ Weights.PerturbVariables)' .* get_error(X_temp0, des_X0));
    problem.x0(1:(N_p*(p-m))) = param;
    [~, ~, ~, X_temp] = simulate_control(t, X_temp0, X_ref, WK_R, WK_L, INSECT, ...
    1, N_single, N_iters, N_per_iter, problem, Weights, ...
    get_args, param_type, p, m, param, 'none', control_net, des_X0);
    X_temp0 = X_temp(end, :)';
    param = sim(control_net, (1 ./ Weights.PerturbVariables)' .* get_error(X_temp0, des_X0));
    problem.x0((1+(N_p*(p-m))):end) = param;
else
    disp('Not using neural net');
end

for period=1:N_periods
    for iter=1:m:N_iters
        idx_opt = (1+(period-1)*N_single+(iter-1)*N_per_iter):(1+(period-1)*N_single+(iter+p-1)*N_per_iter);
        dang0 = dang(:,idx_opt(1));
        varargin = get_args(t(idx_opt(1)), dang(:,idx_opt(1)));
        problem.objective = @(param) cost_fun(param, t(idx_opt), X0, ...
            INSECT, WK_R, WK_L, X_ref(idx_opt,:), Weights, param_type, ...
            get_args, p, N_p, N_per_iter, varargin{:});
%         problem.nonlcon = @(param) param_cons(param, p, N_iters);
		for sol_iter=1:5
            [param, fval, exitflag, output] = solver(problem);
            if ( ((period == 1 && output.iterations >= 15) || (period > 1 && output.iterations >= 5)) ...
					&& (exitflag ~= -2) && (fval < 1e0) && (output.firstorderopt < 1e10))
                fprintf('Exit flag is %d\n', exitflag);
				opt_complete(1+((period-1)*N_iters+(iter-1))/m) = true;
				break;
            else
                if period == 1
                    problem.x0 = 1e-2 * (2*rand(N_p*p, 1)-1) .* problem.ub;
                else
                    problem.x0((1+N_p*(p-m)):end) = 5e-3 * (2*rand(N_p*m, 1)-1) .* problem.ub((1+N_p*(p-m)):end);
                end
            end
        end
        opt_fval(1+((period-1)*N_iters+(iter-1))/m) = fval;
        opt_iter(1+((period-1)*N_iters+(iter-1))/m) = output.iterations;
        opt_firstorder(1+((period-1)*N_iters+(iter-1))/m) = output.firstorderopt;
        
        if p > 1 
            problem.x0 = [param((1+N_p*m):end); zeros(N_p*m,1)];
        end
        param = param(1:(N_p*m));
        opt_param(:,1+((period-1)*N_iters+(iter-1))/m) = param;
        idx_con = (1+(period-1)*N_single+(iter-1)*N_per_iter):(1+(period-1)*N_single+(iter+m-1)*N_per_iter);
        
        idx = idx_con(1:(1+m*N_per_iter));
        X(idx, :) = crgr_xR_control_mex(INSECT, WK_R, WK_L, t(idx), X0, dang0, param, N_p, m, N_per_iter);
        for con=1:m
            param_idx = (1+(con-1)*N_p):(con*N_p);
            param_m = param(param_idx);
            idx = idx_con((1+(con-1)*N_per_iter):(1+con*N_per_iter));

            % [~, X(idx,:)] = ode45(@(t,X) eom_param(INSECT, WK_R, WK_L, t, X, param_m, param_type, varargin{:}), ...
            % t(idx), X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
            for k=idx
                dang(:,k) = param_type(param_m, t(k), varargin{:});
                [WK_R_new, WK_L_new] = get_WK(WK_R, WK_L, dang(:,k));
                [X_dot(:,k), R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k),...
                    theta_A(k), W(:,k), W_R(:,k), W_R_dot(:,k), W_L(:,k),...
                    W_L_dot(:,k), W_A(:,k), W_A_dot(:,k), F_R(:,k), F_L(:,k), M_R(:,k),...
                    M_L(:,k), f_a(:,k), f_g(:,k), f_tau(:,k), tau(:,k), ...
                    Euler_R(:,k), Euler_R_dot(:,k)] = ...
                    eom_QS_xR(INSECT, WK_R_new, WK_L_new, t(k), X(k, :)');
            end
            
            X0 = X(idx(end),:)';
            dang0 = param_type(param_m, t(idx(end)), varargin{:});
            varargin = get_args(t(idx(end)), dang0);
            cost(1+(period-1)*N_iters+(iter-1)+con) = sqrt(...
                sum((Weights.OutputVariables .* (X(idx(end), :) - X_ref(idx(end), :))).^2));
        end
        if isa(control_net, 'network')
            X_temp0 = X0;
            [~, ~, ~, X_temp] = simulate_control(t, X_temp0, X_ref, WK_R, WK_L, INSECT, ...
            1, N_single, N_iters, N_per_iter, problem, Weights, ...
            get_args, param_type, p, m, problem.x0(1:(N_p*(p-m))), 'none', control_net, des_X0);
            X_temp0 = X_temp(end, :)';
            problem.x0((1+(N_p*(p-m))):end) = sim(control_net, (1 ./ Weights.PerturbVariables)' .* get_error(X_temp0, des_X0));
        end
    end
    fprintf('Period number %d completed\n', period);
end

x=X(:,1:3)';
x_dot=X(:,13:15)';

end

function err = get_error(X, Xd)
err = zeros(12, 1);
R = reshape(X(4:12), 3, 3); Rd = reshape(Xd(4:12), 3, 3);
err(1:3) = X(1:3) - Xd(1:3);
err(4:6) = 0.5*vee(Rd'*R - R'*Rd);
% [err(4), err(6), err(5)] = dcm2angle(R'*Rd, 'xzy');
err(7:9) = X(13:15) - Xd(13:15);
err(10:12) = X(16:18) - (R'*Rd*Xd(16:18));
end

function J = cost_fun(param, t, X0, INSECT, WK_R, WK_L, X_ref, Weights, param_type, ...
    get_args, p, N_p, N_per_iter, varargin)
%%
dt = t(2) - t(1);

X = crgr_xR_control_mex(INSECT, WK_R, WK_L, t(1:(1+p*N_per_iter)), X0, zeros(N_p,1), param, N_p, p, N_per_iter);
int_dx = cumsum([zeros(1, 3); (X(:, 1:3) - X_ref(1:(1+p*N_per_iter), 1:3))*dt], 1);
int_dx = int_dx(1:end-1, :);

if (~all(isfinite(X), 'all'))% || any(all(X == 0, 2))
	X_len = sum(all(isfinite(X) == 1, 2));
	if X_len > 1
		J = sqrt(sum((Weights.OutputVariables .* (X(X_len-1, :) - X_ref(X_len-1,:))).^2, 2));
		J = J*(exp(p*N_per_iter+1-X_len) - 1);
	else
		J = 1e100;
	end
    return;
end

preds = (1+N_per_iter):N_per_iter:(1+p*N_per_iter);
J = sum(Weights.PredictionHorizon' .* sqrt(...
        sum((Weights.OutputVariables .* (X(preds, :) - X_ref(preds,:))).^2, 2) + ...
        sum((Weights.OutputVariables(1:3) .* int_dx(preds, :)).^2, 2)));

end

function dang = linear_func(param, t, varargin)
t0 = varargin{1};
dang0 = varargin{2};
dang = dang0 + param*(t - t0);
end

function args = linear_args(t0, dang0)
args{1} = t0;
args{2} = dang0;
end

function [WKR_new, WKL_new] = get_WK(WKR_old, WKL_old, dang)
%%
    WKR_new = WKR_old; WKL_new = WKL_old;
    WKR_new.phi_m = WKR_new.phi_m + dang(1) + dang(3);
    WKL_new.phi_m = WKL_new.phi_m + dang(1) - dang(3);
    WKR_new.theta_0 = WKR_new.theta_0 + dang(2) + dang(5);
    WKL_new.theta_0 = WKL_new.theta_0 + dang(2) - dang(5);
    WKR_new.phi_0 = WKR_new.phi_0 + dang(4);
    WKL_new.phi_0 = WKL_new.phi_0 + dang(4);
    WKR_new.psi_m = WKR_new.psi_m + dang(6);
    WKL_new.psi_m = WKL_new.psi_m - dang(6);
    
    % Include more wing parameters?
    
%     WKR_new.psi_m = WKR_new.psi_m + dang(4) + dang(6);
%     WKL_new.psi_m = WKL_new.psi_m + dang(4) - dang(6);
    
%     WKR_new.theta_A_m = WKR_new.theta_A_m + dang(7);
%     WKL_new.theta_A_m = WKL_new.theta_A_m + dang(7);
end
