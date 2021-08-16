%% Simple NN Model with algorithms : Multilayer Perceptron (2 layer linear activation)
evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
filename = 'iterative_learning';

load('sim_QS_xR_hover_control_opt_200_big', 't', 'X_ref0', ...
    'WK_R', 'WK_L', 'INSECT', 'N_single', 'N_iters', 'N_per_iter', ...
    'get_args', 'param_type', 'm', 'N_dang', 'Weights', 'N_sims', ...
    'X0_pert', 'opt_param', 'des', 'opt_complete');

X0_pert_flat = squeeze(X0_pert(:, 1, :));
opt_param_flat = reshape(permute(opt_param, [1, 3, 2]), 3*N_sims, 6*N_iters);
err = zeros(3*N_sims, 12);
for i=1:(3*N_sims)
    err(i, :) = get_error(X0_pert_flat(i, :)', des.X0)';
end
inputs = (1 ./ Weights.PerturbVariables)' .* err'; targets = opt_param_flat';
inputs = inputs(:, opt_complete); targets = targets(:, opt_complete);

% load('fwuav', 'inputs', 'targets');
N_zero = 1000; % 25, 100
N_data = size(inputs, 2);
inputs(:, (N_data+1):(N_data+N_zero)) = zeros(12, N_zero);
targets(:, (N_data+1):(N_data+N_zero)) = zeros(60, N_zero);

N_features = size(inputs, 1);
[N_outputs, N_data] = size(targets);
N_iters = 5;

%% NN Model
control_net = cascadeforwardnet([36]);
% control_net = cascadeforwardnet([18, 36]);
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
alpha = 0.1;
N_inp = N_features; [N_out, N_neu] = size(control_net.LW{2,1});
solver = @fmincon;
problem.solver = func2str(solver);
problem.options = optimoptions(solver);
problem.options.Algorithm = 'sqp'; % interior-point, sqp, active-set
problem.options.Display = 'iter';
problem.options.UseParallel = true;
problem.options.MaxFunctionEvaluations = 1e5;
problem.options.MaxIterations = 25; % 10 iters
problem.options.ConstraintTolerance = 1e-10;
problem.options.StepTolerance = 1e-10;
% problem.options.OptimalityTolerance = 1e-10;
problem.options.SpecifyObjectiveGradient = false;
problem.options.ObjectiveLimit = 0;

%% Algorithm
y_star = targets;
perf = zeros(N_iters, 1); perf_z = zeros(N_iters, 1);
cons = zeros(N_iters, 1); cons_z = zeros(N_iters, 1);
cs = cell(N_iters, 1); es = cell(N_iters, 1); os = cell(N_iters, 1);

tic;
% diary filename
% PRETRAINING
control_net = init(control_net);
% control_net = configure(control_net, inputs, targets);
control_net = train(control_net, inputs, targets, 'useParallel', 'yes'); %
y = control_net(inputs);
w_star = get_weights(control_net);
w_y = w_star;

for i=1:N_iters
    %% CONSTRAINTS / MASTER STEP
    control_net = init(control_net);
    control_net = train(control_net, inputs, (1-alpha)*y_star + alpha*y, 'useParallel', 'yes'); % 'useParallel', 'yes'
    w0 = get_weights(control_net);
%         w0 = (1-alpha)*w_star + alpha*w_y;

    N_w = length(w0);
    del = 5e-2;
%             problem.x0 = w0 + del * rand(N_w, 1) .* abs(w0);
    problem.x0 = w0;
    problem.lb = w0 - del * abs(w0);
    problem.ub = w0 + del * abs(w0);
    problem.nonlcon = @(w) zero_cons(w, control_net, N_inp, N_neu, N_out);

%%
    % problem.objective = @(w) weight_fun(w, w0);
%     problem.objective = @(w) stability_obj(w, control_net, N_inp, N_neu, N_out, t, ...
%         X_ref0', WK_R, WK_L, INSECT, N_single, N_iters, N_per_iter, N_dang,...
%         get_args, param_type, m, Weights);
%     problem.nonlcon = @(w) stability_error(w, control_net, N_inp, N_neu, N_out, t, ...
%         X_ref0', WK_R, WK_L, INSECT, N_single, N_iters, N_per_iter, N_dang,...
%         get_args, param_type, m, Weights);

%% Period error
    n = 12; epsilon = 1e-6;
    delta0 = epsilon * Weights.PerturbVariables .* eye(n);
%     n = 1; epsilon = 1e-6;
%     delta0 = epsilon * Weights.PerturbVariables';
    problem.objective = @(w) period_error(w, control_net, N_inp, N_neu, N_out, ...
        t, X_ref0', WK_R, WK_L, INSECT, N_single, N_iters, N_per_iter, N_dang, ...
        get_args, param_type, m, Weights, delta0, n);

    [w, fval, exitflag, output] = solver(problem);
	es{i} = exitflag; os{i} = output;
    control_net = update_weights(w, control_net, N_inp, N_neu, N_out);
    control_net_s = control_net;

%%
    % w0 = get_weights(control_net);
    % N_w = length(w0);
    % del = 5e-2;
    % problem.x0 = w0;
    % problem.lb = w0 - del * abs(w0);
    % problem.ub = w0 + del * abs(w0);
    % problem.nonlcon = @(w) zero_cons(w, control_net, N_inp, N_neu, N_out);
    % problem.objective = @(w) weight_fun(w, w0);
	% problem.options.MaxIterations = 25;
	% problem.options.SpecifyObjectiveGradient = true;

    % [w, fval, exitflag, output] = solver(problem);
    % control_net = update_weights(w, control_net, N_inp, N_neu, N_out);
    % control_net_z = control_net;

    z = control_net(inputs);
    cons_z(i) = norm(control_net(zeros(N_features,1)));
    perf_z(i) = perform(control_net, y_star, z);

    %% TRAINING / LEARNING
    control_net = init(control_net);
    control_net = train(control_net, inputs, z, 'useParallel', 'yes'); %
	cs{i} = control_net;
    y = control_net(inputs);
    w_y = get_weights(control_net);
    cons(i) = norm(control_net(zeros(N_features,1)));
    perf(i) = perform(control_net, y_star, y);
end

time_taken = toc;
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
    get_args, param_type, m, Weights, delta0, n)
%%
control_net = update_weights(w, control_net, N_inp, N_neu, N_out);

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
    err(j) = err(j) + sum((Weights.OutputVariables .* (X(end, :) - X(1, :))).^2);
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

function J = stability_obj(varargin)
    [~, J] = stability_error(varargin{:});
end

function [c, ceq] = stability_error(w, control_net, N_inp, N_neu, N_out, ...
    t, X0, WK_R, WK_L, INSECT, N_single, N_iters, N_per_iter, N_p, ...
    get_args, param_type, m, Weights)
%%
c = [];
control_net = update_weights(w, control_net, N_inp, N_neu, N_out);
param = control_net(zeros(N_inp, 1));

X = zeros(1+N_single, 18);

idx_con = 1:(1+N_single);
varargin = get_args(t(idx_con(1)), zeros(N_p, 1));
for con=1:m
    param_idx = (1+(con-1)*N_p):(con*N_p);
    idx = idx_con((1+(con-1)*N_per_iter):(1+con*N_per_iter));
    param_m = param(param_idx);

    [~, X_temp] = ode45(@(t,X) eom_param(INSECT, WK_R, WK_L, t, X, param_m, param_type, varargin{:}), ...
    t(idx), X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
    if size(X_temp,1) ~= (N_per_iter+1)
        ceq = 1e10;
        return;
    else
        X(idx,:) = X_temp;
    end

    X0 = X(idx(end),:)';
    dang0 = param_type(param_m, t(idx(end)), varargin{:});
    varargin = get_args(t(idx(end)), dang0);
end

ceq = sum((Weights.OutputVariables .* (X(end, :) - X(1, :))).^2);
% ceq = X(end, :) - X0';

end

function [J, dJ] = weight_fun(w, w0)
% c1 = 1; c2 = 1;
% J = c1 * norm(w - w0) + c2 * norm(control_net(zeros(N_inp, 1)));
J = 0.5 * (norm(w - w0))^2;

if nargout > 1
    dJ = w - w0;
end

end

function [c, ceq] = zero_cons(w, control_net, N_inp, N_neu, N_out)
c = [];

control_net = update_weights(w, control_net, N_inp, N_neu, N_out);
ceq = control_net(zeros(N_inp, 1));
% ceq = norm(control_net(zeros(N_inp, 1)));
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
