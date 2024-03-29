evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
filename = 'dart_control_dnn';
dnnTrain = dnnTrainFuns;

load('sim_QS_xR_hover_control_opt_bigxR_T2_N30000', 't', 'X_ref0', ...
    'WK_R', 'WK_L', 'INSECT', 'N_single', 'N_iters', 'N_per_iter', ...
    'get_args', 'param_type', 'm', 'N_dang', 'Weights', 'N_sims', ...
    'X0_pert', 'opt_param', 'des', 'opt_complete', 'cost', 'problem', 'solver', 'X_ref');

% X0_pert_flat = squeeze(X0_pert(:, 1, :));
% opt_param_flat = reshape(permute(opt_param, [1, 3, 2]), 3*N_sims, 6*N_iters);
% err = zeros(3*N_sims, 12);
% for i=1:(3*N_sims)
%     err(i, :) = get_error(X0_pert_flat(i, :)', des.X0)';
% end
% inputs = (1 ./ Weights.PerturbVariables)' .* err'; targets = opt_param_flat';
% opt_complete = opt_complete & (cost(:, end) < cost(:, 1));
% % inputs = inputs(:, opt_complete); targets = targets(:, opt_complete);
% % idx = 1:1000;
% idx = 1:450;
% inputs = inputs(:, opt_complete([idx, N_sims+idx, 2*N_sims+idx]));
% targets = targets(:, opt_complete([idx, N_sims+idx, 2*N_sims+idx]));

is_completely_new = false;
alpha_dart = 1e-3; % 1e-4 works, 1e-2 or 1e-6 don't

% load('python_training_data.mat');%, 'input_coil', 'target_coil', 'input_zeros');
idx = 1:5000;
load('python_multiple_input.mat');
% idx = 1:3000;
inputs = input_coil(idx, :)'; targets = target_coil(idx, :)';
inputs_multiple = {squeeze(input_coil(idx, 1, :))'; squeeze(input_coil(idx, 2, :))'};

% N_zero = 1000; % 25, 100
N_data = size(inputs, 2);
N_zero = round(N_data/15);
% N_zero = 500;
N_zero = 900;

inputs(:, (N_data+1):(N_data+N_zero)) = zeros(24, N_zero);
for i=1:length(inputs_multiple)
    inputs_multiple{i}(:, (N_data+1):(N_data+N_zero)) = zeros(12, N_zero);
end
targets(:, (N_data+1):(N_data+N_zero)) = zeros(60, N_zero);
% N_zero = size(input_zeros', 2);
% inputs = [inputs, input_zeros'];
% targets = [targets, zeros(60, N_zero)];

N_features = size(inputs_multiple{1}, 1);
[N_outputs, N_data] = size(targets);

% inputs_orig = inputs; targets_orig = targets;
rand_idx = randperm(N_data);
inputs = inputs(:, rand_idx); targets = targets(:, rand_idx);
% inputs_orig = inputs_multiple;
inputs_mean = {}; inputs_scale = {};
for i=1:length(inputs_multiple)
    inputs_multiple{i} = inputs_multiple{i}(:, rand_idx);
    [inputs_multiple{i}, inputs_mean{i}, inputs_scale{i}] = normalize(inputs_multiple{i}, 2, 'range', [-1, 1]);
end
[targets, targets_mean, targets_scale] = normalize(targets, 2, 'range', [-1, 1]);
large_inp_idx = find(vecnorm(inputs_multiple{1}, 2, 1) > 1);

%% NN Model
drop_prob = 0.0; % 0.25;
layers = [featureInputLayer(12, 'Normalization', 'none')
    
    fullyConnectedLayer(36) % 'Name','fc1'
    leakyReluLayer
    dropoutLayer(drop_prob)
    
%     fullyConnectedLayer(96)
%     leakyReluLayer
%     dropoutLayer(drop_prob)

%     concatenationLayer(1,2)
%     fullyConnectedLayer(60)

%     regressionLayer
    ];
lGraph = layerGraph(layers);
% lGraph = connectLayers(lGraph, 'input', 'concat/in2');

numEpochs = 150; % 50, 150;
learningRate = 1e-3;
executionEnvironment = "cpu";
validationFrequency = 5;
validationCheck = false;
regularizationCoeff = 1e-3;

miniBatchSize = N_data; % numEpochs = 250 (more needed with a larger batch)

[zero_inp1, zero_inp2, zero_out] = dnnTrain.norm_data(dlarray(zeros(12,1), "CB"), ...
    dlarray(zeros(12,1), "CB"), dlarray(zeros(60,1), "CB"), ...
    inputs_mean, inputs_scale, targets_mean, targets_scale);

%%
N_dagger_iters = 5;

rng default;
WK = WK_R;
N_periods = 5; % 5
N_iters = 10; % 2, 4, 10 per period
m = N_iters; % Control Horizon multiplier of iter
p = 2*m;
N_single = 100; % per time period
N_per_iter = N_single / N_iters;
N = N_single*(N_periods-1) + 1 + N_per_iter*p;
T = (N_periods - 1 + p/N_iters)/WK.f;
t = linspace(0,T,N);
dt = t(2)-t(1);
X_ref0 = des.X0;

X_ref = zeros(N, 18);
idx_ref = 1:(1+N_single);
X_ref(idx_ref, :) = crgr_xR_mex(INSECT, WK, WK, t(idx_ref), X_ref0);
idx_ref = (1+N_single):N;
X_ref(idx_ref, :) = X_ref(mod(idx_ref-1, N_single)+1, :);

%%
rng(1);
if is_completely_new
	N_sims = 55;
	scale = 1;
else
	% N_sims = 10;
	% scale = logspace(0, -2, 3);
	% This works
	N_sims = 50;
	scale = logspace(0, -3, 4);
end
N_scale = length(scale);
dX = zeros(N_sims*N_scale, 12);
dX_scale = zeros(N_sims*N_scale, 1);

%% Algorithm
y_star = targets;
perf = zeros(N_dagger_iters, 1); perf_z = zeros(N_dagger_iters, 1);
cons = zeros(N_dagger_iters, 1); cons_z = zeros(N_dagger_iters, 1);
cost = zeros(N_sims*N_scale, N_periods+1);
X_T = zeros(N_sims*N_scale, N_periods+1, 18);

% INITIAL TRAINING
ttotal = tic;

disp('Starting initial training');
[mbq, mbq_val, mbq_full] = dnnTrain.train_val_ds(inputs_multiple, targets, 1, miniBatchSize, "cpu");
[dlnet, netParams] = dnnTrain.initialize_nn_param(layers);
[dlnet, netParams, loss] = dnnTrain.train_custom_network(dlnet, netParams, mbq, ...
    numEpochs, learningRate, executionEnvironment, ...
    mbq_val, validationFrequency, validationCheck,...
    regularizationCoeff);
disp('Finished initial training');

for iter=1:N_dagger_iters
    %% Noise parameter distribution
	if is_completely_new && (iter==1)
		covariance(iter, :, :) = 0.1 * alpha_dart * eye(N_outputs);
	else
		tcovar = tic;
		cov_iter = zeros(N_outputs);
		parfor j = 1:size(inputs_multiple{1}, 2)
            inp1 = dlarray(inputs_multiple{1}(:, j), "CB"); inp2 = dlarray(inputs_multiple{2}(:, j), "CB");
            err_target = targets_scale .* (extractdata(dnnTrain.predictnet(dlnet, netParams, inp1, inp2)) - targets(:, j));
			cov_iter = cov_iter + err_target * err_target';
		end
		covariance(iter, :, :) = alpha_dart * cov_iter / size(inputs_multiple{1}, 2);
		time_covar = toc(tcovar);
    end

    %% New data from noisy expert
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

    opt_complete_iter = zeros(N_sims*N_scale, N_periods, 1, 'logical');
	inputs_iter = zeros(N_features, N_sims*N_scale, N_periods);
	targets_iter = zeros(N_outputs, N_sims*N_scale, N_periods);
    
%     sc = parallel.pool.Constant(RandStream('Threefry', 'Seed', 0));
	cov_iter = double(squeeze(covariance(iter, :, :)));
	disp('Starting data collection for iteration ' + string(iter));
	titer = tic;
	parfor i=1:N_sims*N_scale
%         stream = sc.Value;        % Extract the stream from the Constant
%         stream.Substream = i;
        % rng(i + N_sims*N_scale*(iter-1));

		dX0 = dX(i, :)';
		X0 = [X_ref0(1:3)+ dX0(1:3); reshape(reshape(X_ref0(4:12),3,3)*expmhat(dX0(4:6)), 9, 1); ...
			X_ref0(13:18) + dX0(7:12)];
        
        [cost_i, opt_param_i, ~, ~, ~, opt_complete_i, ~, X_i] = ...
            simulate_opt_control(t, X0, X_ref, WK_R, WK_L, INSECT, ...
            N, N_periods, N_single, N_iters, N_per_iter, problem, solver, Weights, ...
            get_args, param_type, p, m, 0, cov_iter);

		opt_complete_iter(i, :) = opt_complete_i;
		targets_iter(:, i, :) = opt_param_i;
		for j=1:N_periods
			inputs_iter(:, i, j) = (1 ./ Weights.PerturbVariables)' .* get_error(X_i(1+(j-1)*N_single, :)', X_ref0);
		end
        
		% cost(i, :) = cost_arr;
		% X_T(i, :, :) = X(1:N_single:end, :);
    end
	time_iter = toc(titer);
	disp('Finished data collection for iteration ' + string(iter));
	
    opt_complete_iter = all(opt_complete_iter, 2);
    inputs_iter = inputs_iter(:, opt_complete_iter, :);
	targets_iter = targets_iter(:, opt_complete_iter, :);
    
    inputs_iter1 = reshape(inputs_iter(:, :, 1:(end-1)), size(inputs_iter, 1), []);
    inputs_iter2 = reshape(inputs_iter(:, :, 2:end), size(inputs_iter, 1), []);
    targets_iter = reshape(targets_iter(:, :, 2:end), size(targets_iter, 1), []);

    [inputs_iter1, inputs_iter2, targets_iter] = dnnTrain.norm_data(inputs_iter1, inputs_iter2, targets_iter, ...
        inputs_mean, inputs_scale, targets_mean, targets_scale);
	inputs_multiple{1} = [inputs_multiple{1}, extractdata(inputs_iter1)];
    inputs_multiple{2} = [inputs_multiple{2}, extractdata(inputs_iter2)];
	targets = [targets, extractdata(targets_iter)];

    %% Training
	disp('Starting training for iteration ' + string(iter));
	ttrain = tic;
    [mbq, mbq_val, mbq_full] = dnnTrain.train_val_ds(inputs_multiple, targets, 1, size(targets,2), "cpu");
    [dlnet, netParams] = dnnTrain.initialize_nn_param(layers);
    [dlnet, netParams, loss] = dnnTrain.train_custom_network(dlnet, netParams, mbq, ...
        numEpochs, learningRate, executionEnvironment, ...
        mbq_val, validationFrequency, validationCheck,...
        regularizationCoeff);
	time_train = toc(ttrain);
	disp('Finished training for iteration ' + string(iter));

    y = dnnTrain.predictnet(dlnet, netParams, dlarray(inputs_multiple{1}, "CB"), dlarray(inputs_multiple{2}, "CB"));
    zero_pred = dnnTrain.predictnet(dlnet, netParams, zero_inp1, zero_inp2);
    cons(iter) = mse(zero_pred, zero_out);
    perf(iter) = mse(y, targets);
end

time_taken = toc(ttotal);

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
    get_args, param_type, p, m, control_net, covariance)
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
		for sol_iter=1:3
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
        % DART : adding noise to expert control values
        param_temp = mvnrnd(param, covariance);
		param = reshape(param_temp - mean(param_temp), size(param));
        
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
