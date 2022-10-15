%% Simple NN Model with algorithms : Multilayer Perceptron (2 layer linear activation)
evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
filename = 'iterative_learning_coil';

load('sim_QS_xR_hover_control_opt_200_WW', 't', 'X_ref0', ...
    'WK_R', 'WK_L', 'INSECT', 'N_single', 'N_iters', 'N_per_iter', ...
    'get_args', 'param_type', 'm', 'N_dang', 'Weights', 'N_sims', ...
    'X0_pert', 'opt_param', 'des', 'opt_complete', 'cost', 'X_ref');

X0_pert_flat = squeeze(X0_pert(:, 1, :));
opt_param_flat = reshape(permute(opt_param, [1, 3, 2]), 3*N_sims, 6*N_iters);
err = zeros(3*N_sims, 12);
for i=1:(3*N_sims)
    err(i, :) = get_error(X0_pert_flat(i, :)', des.X0)';
end
inputs = (1 ./ Weights.PerturbVariables)' .* err'; targets = opt_param_flat';
opt_complete = opt_complete & (cost(:, end) < cost(:, 1));
% inputs = inputs(:, opt_complete); targets = targets(:, opt_complete);
% idx = 1:1000;
idx = 1:450;
inputs = inputs(:, opt_complete([idx, N_sims+idx, 2*N_sims+idx]));
targets = targets(:, opt_complete([idx, N_sims+idx, 2*N_sims+idx]));
large_inp_idx = find(vecnorm(inputs, 2, 1) > 1);

% N_zero = 1000; % 25, 100
N_data = size(inputs, 2);
N_zero = round(N_data/15);
N_zero = 300;
inputs(:, (N_data+1):(N_data+N_zero)) = zeros(12, N_zero);
targets(:, (N_data+1):(N_data+N_zero)) = zeros(60, N_zero);

N_features = size(inputs, 1);
[N_outputs, N_data] = size(targets);
N_iters = 5;

%% NN Model
control_net = cascadeforwardnet([36]); % 36, 60
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

%% Setup
% rng default;
% alpha = 0.5;
% alphas = linspace(0.1, 0.5, N_iters);
alphas = 0.75 * ones(1, N_iters);
N_inp = N_features; [N_out, N_neu] = size(control_net.LW{2,1});
solver = @fmincon;
problem.solver = func2str(solver);
problem.options = optimoptions(solver);
problem.options.Algorithm = 'interior-point'; % interior-point, sqp
% problem.options.SubproblemAlgorithm = 'cg';
problem.options.Display = 'iter';
problem.options.UseParallel = true;
problem.options.MaxFunctionEvaluations = 1e6;
problem.options.MaxIterations = 10; % 10, 25
% problem.options.ConstraintTolerance = 1e-10;
% problem.options.StepTolerance = 1e-10;
problem.options.ObjectiveLimit = 0;

%% Algorithm
y_star = targets;
perf = zeros(N_iters, 1); perf_z = zeros(N_iters, 1);
perf_yz = zeros(N_iters, 1);
cons = zeros(N_iters, 1); cons_z = zeros(N_iters, 1);
cs = cell(N_iters, 1); es = cell(N_iters, 1); os = cell(N_iters, 1);

ttotal = tic;
% diary filename
% PRETRAINING
control_net = init(control_net);
% control_net = configure(control_net, inputs, targets);
disp('Starting initial training');
tinit = tic;
control_net = train(control_net, inputs, targets, 'useParallel', 'yes'); %
time_init = toc(tinit);
cs_init = control_net;
disp('Finished initial training');
hints = nn.wb_indices(control_net);
y = control_net(inputs);
% w_star = get_weights(control_net);
w_star = getwb(control_net, hints);
w_y = w_star;

for i=1:N_iters
    %% CONSTRAINTS / MASTER STEP
	disp('Starting adjustment training for iteration ' + string(i));
    control_net = init(control_net);
    control_net = train(control_net, inputs, (1-alphas(i))*y_star + alphas(i)*y, 'useParallel', 'yes'); % 'useParallel', 'yes'
	disp('Finished adjustment training for iteration ' + string(i));
%         w0 = (1-alpha)*w_star + alpha*w_y;

%% Period error
    % % w0 = get_weights(control_net);
	% w0 = getwb(control_net, hints);
    % N_w = length(w0);
    % del = 0.05; % 5e-2
% %             problem.x0 = w0 + del * rand(N_w, 1) .* abs(w0);
    % problem.x0 = w0;
    % problem.lb = w0 - del * abs(w0);
    % problem.ub = w0 + del * abs(w0);
    % % problem.nonlcon = @(w) zero_cons(w, control_net, N_inp, N_neu, N_out, hints);
    % n = 12; epsilon = 1e-6;
    % delta0 = epsilon * Weights.PerturbVariables .* eye(n);
% %     n = 1; epsilon = 1e-6;
% %     delta0 = epsilon * Weights.PerturbVariables';
    % problem.objective = @(w) period_error(w, control_net, N_inp, N_neu, N_out, ...
    %     t, X_ref0', WK_R, WK_L, INSECT, N_single, N_iters, N_per_iter, N_dang, ...
    %     get_args, param_type, m, Weights, delta0, n, X_ref, hints);
	% problem.options.MaxIterations = 10; % 10, 25
	% problem.options.SpecifyObjectiveGradient = false;
    % [w, fval, exitflag, output] = solver(problem);
	% es{i} = exitflag; os{i} = output;
    % % control_net = update_weights(w, control_net, N_inp, N_neu, N_out);
	% control_net = setwb(control_net, w, hints);
    % control_net_s = control_net;

%%
    % w0 = get_weights(control_net);
	w0 = getwb(control_net, hints);
    N_w = length(w0);
    del = 0.2;
    problem.x0 = w0;
    % problem.lb = w0 - del * abs(w0);
    % problem.ub = w0 + del * abs(w0);

	% range = control_net.outputs{end}.range;
	% ratio = 2 ./ (range(:,2) - range(:,1));
	% corr_mat = diag(ratio)^2 / N_out;
	corr_mat = eye(60) / N_out;
	F = zeros(control_net.numWeightElements, control_net.numWeightElements);
	parfor k = 1:length(large_inp_idx)
		j = large_inp_idx(k);
	% parfor j = 1:N_data
		de_dwb = defaultderiv('de_dwb', control_net, inputs(:, j), targets(:, j));
		grad_log_pi = - de_dwb * corr_mat * (control_net(inputs(:, j)) - targets(:, j));
		F = F + grad_log_pi * grad_log_pi';
	end
	F = F / length(large_inp_idx);

	% [eig_vec, eig_val] = eig(F, 'vector');
	% idx_eig = (cumsum(eig_val, 'reverse') / sum(eig_val)) <= 0.9; % 0.85, 0.9, 0.95
	% eig_vals = eig_val(idx_eig);
	% eig_vecs = eig_vec(:, idx_eig);
	% problem.objective = @(w) fisher_norm_eig(w, w0, eig_vals, eig_vecs);

	problem.objective = @(w) fisher_norm(w, w0, F);
    % problem.objective = @(w) weight_fun(w, w0);
    problem.nonlcon = @(w) zero_cons(w, control_net, N_inp, N_neu, N_out, hints);

	% zero_inp = zeros(N_inp, 1);
	% problem.objective = @(w) zero_obj(w, control_net, hints, zero_inp);
	% eps_cons = problem.options.ConstraintTolerance;
    % problem.nonlcon = @(w) fisher_norm_cons(w, w0, eig_vals, eig_vecs, eps_cons);

	problem.options.MaxIterations = 20; % 25, 50
	problem.options.SpecifyObjectiveGradient = true;
	% problem.options.SpecifyConstraintGradient = true;
	disp('Starting optimization for iteration ' + string(i));
	topt = tic;
    [w, fval, exitflag, output] = solver(problem);
	es{i} = exitflag; os{i} = output;
	time_opt = toc(topt)
	disp('Finished optimization for iteration ' + string(i));
    % control_net = update_weights(w, control_net, N_inp, N_neu, N_out);
	control_net = setwb(control_net, w, hints);
    control_net_z = control_net;

%% TRAINING / LEARNING
    z = control_net(inputs);
    cons_z(i) = norm(control_net(zeros(N_features,1)));
    perf_z(i) = perform(control_net, y_star, z);

	disp('Starting unconstrained training for iteration ' + string(i));
    control_net = init(control_net);
    control_net = train(control_net, inputs, z, 'useParallel', 'yes'); %
	disp('Finished unconstrained training for iteration ' + string(i));
	cs{i} = control_net;
    y = control_net(inputs);
    % w_y = get_weights(control_net);
	w_y = getwb(control_net, hints);
    cons(i) = norm(control_net(zeros(N_features,1)));
    perf(i) = perform(control_net, y_star, y);
    perf_yz(i) = perform(control_net, z, y);
end

time_taken = toc(ttotal);
disp('Number of weight params near their bounds are,');
N_bad_param = sum(abs(abs((w-w0)./w0) - del) <= del * 1e-2);
disp(N_bad_param);
% diary off;

%% Performance and results
% plotregression(targets, control_net(inputs));
% plotregression(targets(:, tr.testInd), control_net(inputs(:, tr.testInd)));

allvars = whos;
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
save(filename, allvars(tosave).name)
evalin('base',['load ' filename]);

%%
function J = period_error(w, control_net, N_inp, N_neu, N_out, ...
    t, X0, WK_R, WK_L, INSECT, N_single, N_iters, N_per_iter, N_p, ...
    get_args, param_type, m, Weights, delta0, n, X_ref, hints)
%%
% control_net = update_weights(w, control_net, N_inp, N_neu, N_out);
control_net = setwb(control_net, w, hints);

idx_con = 1:(1+N_single);
err = zeros(n, 1);

for j=1:n
%     j = 1;
%     X = zeros(1+N_single, length(X0));
    dR = delta0(4:6, j);
    X_del0 = [X0(1:3)+delta0(1:3, j); reshape(reshape(X0(4:12),3,3)*expmhat(dR), 9, 1); ...
        X0(13:18) + delta0(7:12, j)];
    param = control_net((1 ./ Weights.PerturbVariables)' .* get_error(X_del0, X0));
%     varargin = get_args(t(idx_con(1)), zeros(N_p, 1));

    X = crgr_xR_control_mex(INSECT, WK_R, WK_L, t(idx_con), X_del0, zeros(N_p,1), param, N_p, m, N_per_iter);
    if ~all(isfinite(X), 'all')
        err(j) = 1/eps;
    end
%     for con=1:m
%         param_idx = (1+(con-1)*N_p):(con*N_p);
%         idx = idx_con((1+(con-1)*N_per_iter):(1+con*N_per_iter));
%         param_m = param(param_idx);
%
%         [~, X_temp] = ode45(@(t,X) eom_param(INSECT, WK_R, WK_L, t, X, param_m, param_type, varargin{:}), ...
%         t(idx), X_del0, odeset('AbsTol',1e-6,'RelTol',1e-6));
%         if size(X_temp,1) ~= (N_per_iter+1)
%             err(j) = 1/eps;
%             break;
%         else
%             X(idx,:) = X_temp;
%         end
%
%         X_del0 = X(idx(end),:)';
%         dang0 = param_type(param_m, t(idx(end)), varargin{:});
%         varargin = get_args(t(idx(end)), dang0);
%     end

    % err(j) = err(j) + sum((Weights.OutputVariables .* (X(end, :) - X0')).^2);
    % err(j) = err(j) + sum((Weights.OutputVariables .* (X(end, :) - X(1, :))).^2);

	% int_dx = cumsum([zeros(1, 3); (X(:, 1:3) - X_ref(1:(1+m*N_per_iter), 1:3))*(t(2)-t(1))], 1);
	% int_dx = int_dx(1:end-1, :);
	preds = (1+N_per_iter):N_per_iter:(1+m*N_per_iter);
	J = sum(Weights.PredictionHorizon(1:10)' .* sqrt(...
			sum((Weights.OutputVariables .* (X(preds, :) - X_ref(preds,:))).^2, 2))); % + ...
			% sum((Weights.OutputVariables(1:3) .* int_dx(preds, :)).^2, 2)));
	err(j) = err(j) + J;
end

J = sum(err);

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

function [J, dJ] = weight_fun(w, w0)
% c1 = 1; c2 = 1;
% J = c1 * norm(w - w0) + c2 * norm(control_net(zeros(N_inp, 1)));
J = 0.5 * (norm(w - w0))^2;

if nargout > 1
    dJ = w - w0;
end

end

function [J, dJ] = fisher_norm(w, w0, F)
J = 0.5 * (w-w0)' * F * (w-w0);

if nargout > 1
    dJ = F * (w-w0);
end

end

function [J, dJ] = fisher_norm_eig(w, w0, eig_vals, eig_vecs)
e_vecs_x_w = eig_vecs' * (w-w0);
J = 0.5 * eig_vals' * (e_vecs_x_w).^2;

if nargout > 1
	dJ = eig_vecs * (eig_vals .* e_vecs_x_w);
end

end

function [J] = zero_obj(w, control_net, hints, zero_inp)
control_net = setwb(control_net, w, hints);
J = 0.5 * norm(control_net(zero_inp))^2;
end

function [c, ceq] = fisher_norm_cons(w, w0, eig_vals, eig_vecs, eps_cons)
% function [c, ceq, gradc, gradceq] = fisher_norm_cons(w, w0, eig_vals, eig_vecs, eps_cons)
% if nargout > 2
% 	e_vecs_x_w = eig_vecs' * (w-w0);
% 	c = 0.5 * eig_vals' * (e_vecs_x_w).^2 - eps_cons;
% 	gradc = eig_vecs * (eig_vals .* e_vecs_x_w);
% 	gradceq = [];
% else
% 	c = 0.5 * eig_vals' * (eig_vecs' * (w-w0)).^2 - eps_cons;
% end

c = 0.5 * eig_vals' * (eig_vecs' * (w-w0)).^2 - eps_cons;
ceq = [];

end

% function [c, ceq, gradc, gradceq] = zero_cons(w, control_net, N_inp, N_neu, N_out, hints)
function [c, ceq] = zero_cons(w, control_net, N_inp, N_neu, N_out, hints)
c = [];

% control_net = update_weights(w, control_net, N_inp, N_neu, N_out);
control_net = setwb(control_net, w, hints);
ceq = control_net(zeros(N_inp, 1));
% zero_inp = zeros(N_inp, 1);
% zero_out = control_net(zero_inp);
% ceq = 0.5 * norm(zero_out)^2;

% if nargout > 2
% 	gradc = [];
% 	% gradceq = defaultderiv('de_dwb', control_net, zero_inp, zero_out);
% 	gradceq = defaultderiv('de_dwb', control_net, zero_inp, zero_out) * zero_out;
% end

end

function control_net = update_weights(w, control_net, N_inp, N_neu, N_out)
%%
control_net.IW{1} = reshape(w(1:N_neu*N_inp), N_neu, N_inp);
idx = N_neu*N_inp;
control_net.IW{2} = reshape(w((idx+1):(idx+N_out*N_inp)), N_out, N_inp);
idx = idx + N_out*N_inp;
control_net.LW{2,1} = reshape(w((idx+1):(idx+N_out*N_neu)), N_out, N_neu);
idx = idx + N_out*N_neu;
control_net.b{1} = reshape(w((idx+1):(idx+N_neu)), N_neu, 1);
idx = idx+N_neu;
control_net.b{2} = reshape(w((idx+1):(idx+N_out)), N_out, 1);
idx = idx+N_out;
end

function w = get_weights(control_net)
w = [reshape(control_net.IW{1}, numel(control_net.IW{1}), 1);
    reshape(control_net.IW{2}, numel(control_net.IW{2}), 1);
    reshape(control_net.LW{2,1}, numel(control_net.LW{2,1}), 1);
    control_net.b{1};
    control_net.b{2}];
end

function dXdt = eom_param(INSECT, WK_R, WK_L, t, X, param, param_type, varargin)
%%
dang = param_type(param, t, varargin{:});
[WK_R, WK_L] = get_WK(WK_R, WK_L, dang);
dXdt = eom_QS_xR(INSECT, WK_R, WK_L, t, X);
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
