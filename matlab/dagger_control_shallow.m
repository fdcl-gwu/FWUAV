evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
filename = 'dagger_control_shallow_sectional';

load('sim_QS_xR_hover_control_opt_bigxR_T2_N30000', 't', 'X_ref0', ...
    'WK_R', 'WK_L', 'INSECT', 'N_single', 'N_iters', 'N_per_iter', ...
    'get_args', 'param_type', 'm', 'N_dang', 'Weights', 'N_sims', ...
    'X0_pert', 'opt_param', 'des', 'opt_complete', 'cost', 'problem', 'solver', 'X_ref');

% Weights.OutputVariables(1:3) = 40; % Playing with the tuning of weights
% in optimal control for better position tracking

% This NN controller can be used instead of the expert optimal control,
% speeds up the whole process
expert = load('iterative_learning_coil_N5000.mat', 'control_net');
control_net_expert = expert.control_net;

expert = load('dagger_control_shallow_sectional_predicttrain_worldpos_regW_varylongroll_newbuff_25ep_buffer_estdata.mat', 'control_net');
control_net_1env = expert.control_net;

% load('python_training_data.mat');%, 'input_coil', 'target_coil', 'input_zeros');
idx = 1:5000;
% idx = 1:2500;
load('python_traj_data_worldpos.mat');

input_coil = reshape(input_coil, [], 12);
input_est = reshape(input_est, [], 12);
target_coil = reshape(target_coil, [], 60);
% idx = 1:1500;
inputs = input_coil(idx, :)'; targets = double(target_coil(idx, :)');
inputs = input_est(idx, :)';

% N_zero = 1000; % 25, 100
N_data = size(inputs, 2);
N_zero = round(N_data/15);
N_zero = 300;
% N_zero = 900;
N_zero = 0;

inputs(:, (N_data+1):(N_data+N_zero)) = zeros(12, N_zero);
targets(:, (N_data+1):(N_data+N_zero)) = zeros(60, N_zero);
% N_zero = size(input_zeros', 2);
% inputs = [inputs, input_zeros'];
% targets = [targets, zeros(60, N_zero)];

N_features = size(inputs, 1);
[N_outputs, N_data] = size(targets);

% inputs_orig = inputs; targets_orig = targets;
rand_idx = randperm(N_data);
inputs = inputs(:, rand_idx); targets = targets(:, rand_idx);

large_inp_idx = find(vecnorm(inputs, 2, 1) > 1);

%% Pose server
client = rossvcclient('/get_estimate', 'DataFormat', 'struct');
req = rosmessage(client);
gravity = [0, 0, INSECT.g]';

client_train = rossvcclient('/train_model', 'DataFormat', 'struct');
req_train = rosmessage(client_train);

%% NN Model
control_net = cascadeforwardnet([36]); % 36, 60
% control_net = feedforwardnet([32, 96]);
% control_net.inputConnect(end) = 1;
control_net = configure(control_net, inputs, targets);
trainFcn = 'trainbr'; % trainbr, trainrp
control_net.trainFcn = trainFcn; % trainbr, trainrp
switch trainFcn
    case 'trainbr'
        control_net.trainParam.epochs = 15; % 25
    case 'trainrp'
        control_net.trainParam.epochs = 20000;
        control_net.trainParam.max_fail = 20;
        control_net.performParam.regularization = 0.1; % if not trainbr
    case 'trainlm'
        control_net.trainParam.epochs = 100;
%         control_net.trainParam.max_fail = 20;
        control_net.performParam.regularization = 1e-2;
end
control_net.performParam.normalization = 'standard';
% for i=1:control_net.numInputs
%     control_net.inputs{i}.processFcns = {'mapminmax', 'processpca'};
% end
control_net.layers{1:(end-1)}.transferFcn = 'leakyrelu'; % leakyrelu, tansig, poslin, purelin
control_net.divideParam.trainRatio = 80/100;
control_net.divideParam.valRatio = 10/100;
control_net.divideParam.testRatio = 10/100;

%%
N_dagger_iters = 10;
N_period_pi = 3; % For optimal control expert labeling
% N_period_pi = 10;
% N_period_pi = 1;

rng default;
WK = WK_R;
% N_periods = 10; % 5, 10, 15
N_periods = 100;
N_iters = 10; % 2, 4, 10 per period
m = N_iters; % Control Horizon multiplier of iter
p = 2*m;
N_single = 100; % per time period
N_per_iter = N_single / N_iters;
N = N_single*(max(N_periods, N_period_pi)-1) + 1 + N_per_iter*p;
T = (max(N_periods, N_period_pi) - 1 + p/N_iters)/WK.f;
t = linspace(0,T,N);
dt = t(2)-t(1);
X_ref0 = des.X0;

X_ref = zeros(N, 18);
idx_ref = 1:(1+N_single);
X_ref(idx_ref, :) = crgr_xR_mex(INSECT, WK, WK, t(idx_ref), X_ref0);
idx_ref = (1+N_single):N;
X_ref(idx_ref, :) = X_ref(mod(idx_ref-1, N_single)+1, :);

[~, acc] = crgr_xR_control_mex(INSECT, WK_R, WK_L, t(1:2), X_ref0, zeros(6,1), 0, 0, 0, 0);
acc_ref = reshape(X_ref0(4:12),3,3)' * (acc(1, 1:3)' -  gravity);

%%
rng(1);
% N_sims = 10;
% N_sims = 5;
N_environments = 1;
N_sims_per_env = 5;
N_sims = N_sims_per_env * N_environments;
scale = logspace(0, -2, 3);
% scale = logspace(-1, -3, 3); % N_sims = 30
N_scale = length(scale);

N_val_sims = 150;
val_X = zeros(18, N_val_sims*N_scale);
val_input = zeros(12, N_val_sims*N_scale);
for i = 1:N_val_sims
    dx = 2*rand(1,3)-1; dx = rand(1) * dx / norm(dx);
    dtheta = 2*rand(1,3)-1; dtheta = rand(1) * dtheta / norm(dtheta);
    dx_dot = 2*rand(1,3)-1; dx_dot = rand(1) * dx_dot / norm(dx_dot);
    domega = 2*rand(1,3)-1; domega = rand(1) * domega / norm(domega);
    for j = 1:N_scale
	    dX0 = (scale(j) * Weights.PerturbVariables .* [dx, dtheta, dx_dot, domega])';
        X0 = [X_ref0(1:3)+ dX0(1:3); reshape(reshape(X_ref0(4:12),3,3)*expmhat(dX0(4:6)), 9, 1); ...
			X_ref0(13:18) + dX0(7:12)];
        val_X(:, i+(j-1)*N_val_sims) = X0;
        val_input(:, i+(j-1)*N_val_sims) = (1 ./ Weights.PerturbVariables)' .* get_error(X0, X_ref0);
    end
end
dX_order = randperm(N_val_sims*N_scale);
val_X = val_X(:, dX_order);
val_input = val_input(:, dX_order);
val_est_input = val_input;

%% Algorithm
perf = zeros(N_dagger_iters, 1); perf_z = zeros(N_dagger_iters, 1);
cons = zeros(N_dagger_iters, 1); cons_z = zeros(N_dagger_iters, 1);
diff_score = cell(N_dagger_iters, 1); diff_pose_score = cell(N_dagger_iters, 1);

ttotal = tic;

disp('Starting initial training');
control_net = init(control_net);
% control_net = configure(control_net, inputs, targets);
control_net = train(control_net, inputs, targets, 'useParallel', 'yes'); % 'useParallel', 'yes'
disp('Finished initial training');
% load('dagger_control_shallow_sectional.mat', 'control_net', 'inputs', 'targets');
% control_net = control_net_1env;
cs_init = control_net;
y_val_init = control_net(val_est_input);

for iter=1:N_dagger_iters
    N_periods = 15 + round((iter-1)/N_dagger_iters * 50);
%     N_periods = 36;
%     N_periods = 15 + round((10+iter-1)/N_dagger_iters * 50);

	%% Rollout
    dX = zeros(N_sims_per_env*N_scale, 12);
    dX_scale = zeros(N_sims_per_env*N_scale, 1);
    for i = 1:N_sims_per_env
	    dx = 2*rand(1,3)-1; dx = rand(1) * dx / norm(dx);
	    dtheta = 2*rand(1,3)-1; dtheta = rand(1) * dtheta / norm(dtheta);
	    dx_dot = 2*rand(1,3)-1; dx_dot = rand(1) * dx_dot / norm(dx_dot);
	    domega = 2*rand(1,3)-1; domega = rand(1) * domega / norm(domega);
	    for j = 1:N_scale
		    dX(i+(j-1)*N_sims_per_env,:) = scale(j) * Weights.PerturbVariables .* [dx, dtheta, dx_dot, domega];
		    dX_scale(i+(j-1)*N_sims_per_env) = scale(j) * norm([dx, dtheta, dx_dot, domega]);
	    end
    end
    dX = repmat(dX, N_environments, 1);
%     dX_order = randperm(N_sims*N_scale);
%     dX = dX(dX_order, :);
%     dX_scale = dX_scale(dX_order);

	disp('Starting rollout for iteration ' + string(iter));
	troll = tic;
    cost = zeros(N_sims*N_scale, N_periods+1);
    X_T = zeros(N_sims*N_scale, N_periods+1, 18);
    inputs0_pi = zeros(N_sims*N_scale, N_periods, 12);
    acc0_pi = zeros(N_sims*N_scale, N_periods, 3);
    WorldIdx = zeros(N_sims*N_scale, N_periods);
	for i=1:N_sims*N_scale
		X = zeros(1+N_single*N_periods, 18);
		dang0 = zeros(N_dang, 1);
		cost_arr = zeros(1+N_periods, 1);

		dX0 = dX(i, :)';
		X0 = [X_ref0(1:3)+ dX0(1:3); reshape(reshape(X_ref0(4:12),3,3)*expmhat(dX0(4:6)), 9, 1); ...
			X_ref0(13:18) + dX0(7:12)];
        [~, acc] = crgr_xR_control_mex(INSECT, WK_R, WK_L, t(1:2), X0, dang0, 0, 0, 0, 0);
        Acc0 = reshape(X0(4:12),3,3)' * (acc(1, 1:3)' -  gravity);
		cost_arr(1) = sqrt(sum((Weights.OutputVariables' .* (X0 - X_ref0)).^2));

		for period=1:N_periods
			idx = (1+(period-1)*N_single):(1+period*N_single);

			req.X = X0;
            req.Acc = Acc0;
%             req.Period = int64(period);
            req.Period = int64(N_periods);
            resp = call(client, req);
            Xest = reshape(resp.Xest, 18, 1);
%             param = reshape(resp.ControlParam, 60, 1);

			param = control_net((1 ./ Weights.PerturbVariables)' .* get_error(Xest, X_ref0));

			[X(idx, :), acc] = crgr_xR_control_mex(INSECT, WK_R, WK_L, t(idx), X0, dang0, param, N_dang, m, N_per_iter);

			X0 = X(idx(end),:)';
            Acc0 = reshape(X0(4:12),3,3)' * (acc(end, 1:3)' - gravity);
            acc0_pi(i, period, :) = Acc0;
            inputs0_pi(i, period, :) = (1 ./ Weights.PerturbVariables)' .* get_error(X0, X_ref0);

			dang0 = dang0 + sum(reshape(param, 6, 10), 2) / WK_R.f / m;
			cost_arr(period+1) = sqrt(sum((Weights.OutputVariables .* (X(idx(end), :) - X_ref0')).^2));
		end
		cost(i, :) = cost_arr;
		X_T(i, :, :) = X(1:N_single:end, :);
        WorldIdx(i, :) = resp.WorldIdx; % assumption is rollout of all periods happens in the same environment
	end
	time_roll = toc(troll);
	disp('Finished rollout for iteration ' + string(iter));
    % To switch to a new environment in next steps
    req.Period = int64(-1);
    resp = call(client, req);

	%% Expert labels
% 	X0_pi = reshape(X_T(:, 2:end, :), [], 18);
%     N_pi = size(X0_pi, 1);
    X0_pi = reshape(permute(X_T(:, 2:end, :), [3,2,1]), 18, []);
    acc0_pi = reshape(permute(acc0_pi, [3,2,1]), 3, []);
    WorldIdx = reshape(WorldIdx', 1, []);
    N_pi = size(X0_pi, 2);
	is_bounded = zeros(N_pi, 1, 'logical');
    for i=1:N_pi
%         X0 = X0_pi(i, :)';
        X0 = X0_pi(:, i);
        if norm((1 ./ Weights.PerturbVariables)' .* get_error(X0, X_ref0)) <= 1
            is_bounded(i) = true;
        end
    end

%     X0_pi = X0_pi(is_bounded, :);
%     N_pi = size(X0_pi, 1);
    X0_pi = X0_pi(:, is_bounded);
    acc0_pi = acc0_pi(:, is_bounded);
    WorldIdx = WorldIdx(is_bounded);
    N_pi = size(X0_pi, 2);

    %%
	opt_complete_pi = zeros(N_pi, N_period_pi, 'logical');
	inputs_pi = zeros(N_features, N_pi, N_period_pi);
	targets_pi = zeros(N_outputs, N_pi, N_period_pi);
    X_pi = zeros(18, N_pi, N_period_pi);
    acc_pi = zeros(3, N_pi, N_period_pi);

	disp('Starting expert labeling for iteration ' + string(iter));
	titer = tic;
	parfor i=1:N_pi
% 		X0 = X0_pi(i, :)';
        X0 = X0_pi(:, i);

		[~, opt_param, ~, ~, ~, opt_complete, ~, X_i, ~, ~, x_ddot_i] = ...
			simulate_opt_control(t, X0, X_ref, WK_R, WK_L, INSECT, ...
			N, N_period_pi, N_single, N_iters, N_per_iter, problem, solver, Weights, ...
			get_args, param_type, p, m, 0);

        opt_complete_pi(i, :) = opt_complete;
	    targets_pi(:, i, :) = opt_param;
        for j=1:N_period_pi
            Xper = X_i(1+(j-1)*N_single, :)';
			inputs_pi(:, i, j) = (1 ./ Weights.PerturbVariables)' .* get_error(Xper, X_ref0);
            X_pi(:, i, j) = Xper;
            acc_pi(:, i, j) = reshape(Xper(4:12),3,3)' * (x_ddot_i(:, 1+(j-1)*N_single) -  gravity);
        end
    end
    % Approximate expert labeling
% 	parfor i=1:N_pi
% 		X = zeros(1+N_single*N_period_pi, 18);
% 		dang0 = zeros(N_dang, 1);
%         X0 = X0_pi(:, i);
% 
% 		for period=1:N_period_pi
% 			idx = (1+(period-1)*N_single):(1+period*N_single);
% 
% 			input_X0 = (1 ./ Weights.PerturbVariables)' .* get_error(X0, X_ref0);
%             param = control_net_expert(input_X0);
% 			[X(idx, :), acc] = crgr_xR_control_mex(INSECT, WK_R, WK_L, t(idx), X0, dang0, param, N_dang, m, N_per_iter);
% 
%             X_pi(:, i, period) = X0;
%             acc_pi(:, i, period) = reshape(X0(4:12),3,3)' * (acc(1, 1:3)' - gravity);
% 
%             inputs_pi(:, i, period) = input_X0;
% 	        targets_pi(:, i, period) = param;
%             opt_complete_pi(i, period) = true;
% 
%             X0 = X(idx(end),:)';
% 			dang0 = dang0 + sum(reshape(param, 6, 10), 2) / WK_R.f / m;
% 		end
%     end

	time_iter = toc(titer);
	disp('Finished expert labeling for iteration ' + string(iter));
    opt_complete_pi = all(opt_complete_pi, 2);
    inputs_pi = inputs_pi(:, opt_complete_pi, :);
	targets_pi = targets_pi(:, opt_complete_pi, :);
    WorldIdx = WorldIdx(opt_complete_pi);
    WorldIdx = repelem(WorldIdx, N_period_pi);

    X_pi = reshape(permute(X_pi(:, opt_complete_pi, :), [1,3,2]), size(X_pi, 1), []);
    acc_pi = reshape(permute(acc_pi(:, opt_complete_pi, :), [1,3,2]), size(acc_pi, 1), []);
    inputs_est_pi = zeros(12, size(X_pi, 2));
    targets_est_pi = reshape(permute(targets_pi, [1,3,2]), size(targets_pi, 1), []);
    N_pi = size(X_pi, 2);

    %%
    header_data = {'N_periods', N_pi, 'N_sims', 1, 'N_dang', 6, 'N_iters', 10};
    file_to_save = 'poses_dagger.csv';
    writecell(header_data, file_to_save)
    writematrix(char("---------------------"), file_to_save, 'WriteMode', 'append')
    writematrix(X_ref0', file_to_save, 'WriteMode', 'append')
    writematrix(char("---------------------"), file_to_save, 'WriteMode', 'append')
    writematrix(X_pi', file_to_save, 'WriteMode','append')
    %
    file_to_save = 'poses_acc_dagger.csv';
    writecell(header_data, file_to_save)
    writematrix(char("---------------------"), file_to_save, 'WriteMode', 'append')
    writematrix(acc_ref', file_to_save, 'WriteMode', 'append')
    writematrix(char("---------------------"), file_to_save, 'WriteMode', 'append')
    writematrix(acc_pi', file_to_save, 'WriteMode','append')
    %
    file_to_save = 'poses_control_param_dagger.csv';
    writecell(header_data, file_to_save)
    writematrix(char("---------------------"), file_to_save, 'WriteMode', 'append')
    writematrix(zeros(1, 60), file_to_save, 'WriteMode', 'append')
    writematrix(char("---------------------"), file_to_save, 'WriteMode', 'append')
    writematrix(targets_est_pi', file_to_save, 'WriteMode','append')
    %
    save('poses_world.mat', 'WorldIdx');

    tpose = tic;
    req_train.Data = true;
    resp_train = call(client_train, req_train);
    assert(resp_train.Success);
    time_pose = toc(tpose);
    est_data = load('dagger_python.mat', 'inputs_est', 'X_est');

    %%
    for i=1:N_pi
% 	    req.X = X_pi(:, i);
%         req.Acc = acc_pi(:, i);
%         req.Period = int64(0);
%         req.WorldIdx = int64(WorldIdx(i));
%         resp = call(client, req);
%         Xest = reshape(resp.Xest, 18, 1);
        Xest = est_data.X_est(i, :)'; % From the trained model directly
        inputs_est_pi(:, i) = (1 ./ Weights.PerturbVariables)' .* get_error(Xest, X_ref0);
    end
%     req.Period = int64(-1);
%     resp = call(client, req);
  
    inputs = [inputs, inputs_est_pi];
%     inputs = [inputs, reshape(inputs_pi, size(inputs_pi, 1), [])];
    targets = [targets, targets_est_pi];

    %%
%     N_old = 500; N_end = 4000;
    N_old = 3000; N_end = size(inputs_est_pi, 2);
%     N_old = 2000; N_end = 2000;
    N_input_data = size(inputs, 2);
    rand_idx = randperm(max(N_input_data - N_end, N_old));
    N_buffer_start = max(N_input_data-N_end+1, length(rand_idx)+1);
    inputs_new = [inputs(:, rand_idx(1:N_old)), inputs(:, N_buffer_start:N_input_data)];
    targets_new = [targets(:, rand_idx(1:N_old)), targets(:, N_buffer_start:N_input_data)];

    dang0 = zeros(N_dang, 1);
    val_est_prev = val_est_input;
    val_est_input = zeros(size(val_input));
    for i=1:size(val_input, 2)
	    req.X = val_X(:, i);
        [~, acc] = crgr_xR_control_mex(INSECT, WK_R, WK_L, t(1:2), req.X, dang0, 0, 0, 0, 0);
        req.Acc  = reshape(req.X(4:12),3,3)' * (acc(1, 1:3)' -  gravity);
        req.Period = int64(N_periods);
%         req.WorldIdx = int64(0);
        resp = call(client, req);
        Xest = reshape(resp.Xest, 18, 1);
        val_est_input(:, i) = (1 ./ Weights.PerturbVariables)' .* get_error(Xest, X_ref0);
    end
    req.Period = int64(-1);
    resp = call(client, req);

    diff_pose_score{iter} = [
        mean(vecnorm(val_est_input - val_est_prev, 2, 1)),...
        mean(vecnorm(val_est_input - val_est_prev, 2, 1) ./ vecnorm(val_est_prev, 2, 1)),...
        mse(val_est_input, val_est_prev),...
        1 - norm(val_est_input * val_est_prev', 'fro')^2 / (norm(val_est_input * val_est_input', 'fro') * norm(val_est_prev * val_est_prev', 'fro')),...
        ];
    
    control_net_prev = control_net;
    y_prev = control_net_prev(val_est_input);

    %% TRAINING
	disp('Starting training for iteration ' + string(iter));
	ttrain = tic;
%     control_net.trainParam.epochs = 5;
    control_net = init(control_net);
%     control_net = train(control_net, inputs_new, targets_new, 'useParallel', 'yes'); %
%     reg_weight = 10;
%     reg_weight = 10^0.5;
%     EW = [ones(1, size(inputs_new, 2)), reg_weight*ones(1, size(val_est_input, 2))];
%     control_net = train(control_net, [inputs_new, val_est_input], [targets_new, y_prev],...
%         [], [], EW, 'useParallel', 'yes'); %
    %
% %     reg_weight = 0.1;
    reg_weight = 0.25;

    EW = [reg_weight*ones(1, N_old), ones(1, N_end)];
    control_net = train(control_net, ...
        [inputs(:, rand_idx(1:N_old)), inputs(:, N_buffer_start:N_input_data)], ...
        [control_net_prev(inputs(:, rand_idx(1:N_old))), targets(:, N_buffer_start:N_input_data)], ...
        [], [], EW, 'useParallel', 'yes'); %
    
%     N_old = 2500; N_end = 2500;
%     N_est = size(inputs_est_pi, 2);
%     rand_est = randperm(N_est);
%     rand_est = rand_est(1:min(N_est, N_end));
%     N_input_data = size(inputs, 2);
%     rand_idx = randperm(max(N_input_data - N_est, N_old));
% 
%     EW = [reg_weight*ones(1, N_old), ones(1, length(rand_est))];
%     control_net = train(control_net, ...
%         [inputs(:, rand_idx(1:N_old)), inputs_est_pi(:, rand_est)], ...
%         [control_net_prev(inputs(:, rand_idx(1:N_old))), targets_est_pi(:, rand_est)], ...
%         [], [], EW, 'useParallel', 'yes'); %

	time_train = toc(ttrain);
	disp('Finished training for iteration ' + string(iter));
    save('dagger_control_shallow_sectional.mat', 'control_net', '-append');

    %% Performance
    y = control_net(inputs);
    cons(iter) = norm(control_net(zeros(N_features,1)));
    perf(iter) = mse(y, targets);

    y_val = control_net(val_est_input);
    diff_score{iter} = [
        mean(vecnorm(y_val - y_prev, 2, 1)),...
        mean(vecnorm(y_val - y_prev, 2, 1) ./ vecnorm(y_prev, 2, 1)),...
        mse(y_val, y_prev),...
        1 - norm(y_val * y_prev', 'fro')^2 / (norm(y_val * y_val', 'fro') * norm(y_prev * y_prev', 'fro')),...
        mean(vecnorm(y_val - y_val_init, 2, 1)),...
        mean(vecnorm(y_val - y_val_init, 2, 1) ./ vecnorm(y_val_init, 2, 1)),...
        mse(y_val, y_val_init),...
        1 - norm(y_val * y_val_init', 'fro')^2 / (norm(y_val * y_val', 'fro') * norm(y_val_init * y_val_init', 'fro'))
        ];
end

time_taken = toc(ttotal);

%%
allvars = whos;
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
save(filename, allvars(tosave).name)
evalin('base',['load ' filename]);

%%
function [cost, opt_param, opt_fval, opt_iter, opt_firstorder, opt_complete, dang, ...
    X, x, x_dot, x_ddot, R, Q_R, Q_L, Q_A, theta_A, W, W_R, W_R_dot, W_L, ...
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
					&& (exitflag ~= -2) && (fval < 2e0) && (output.firstorderopt < 1e10))
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
x_ddot=X_dot(13:15,:);

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
