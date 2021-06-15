%% Simple NN Model with algorithms : Multilayer Perceptron (2 layer linear activation)
evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
filename = 'iterative_learning';

load('sim_QS_xR_hover_control_opt_mc', 'inputs', 'targets', 't', 'X_ref0', ...
    'WK_R', 'WK_L', 'INSECT', 'N_single', 'N_iters', 'N_per_iter', ...
    'get_args', 'param_type', 'm', 'N_dang', 'Weights');
% load('fwuav', 'inputs', 'targets');
N_features = size(inputs, 1);
[N_outputs, N_data] = size(targets);
N_iters = 1;

%% NN Model
MLP_net = fitnet([N_features*2]);
MLP_net = configure(MLP_net, inputs, targets);
MLP_net.trainFcn = 'trainbr';% train{lm,br,bfg}
MLP_net.trainParam.epochs = 30; % 30
%  MLP_net.performParam.regularization = 0.1; % if not trainbr
% MLP_net.performParam.normalization = 'standard';
MLP_net.inputs{1}.processFcns = {'mapminmax', 'processpca'};
%             MLP_net.inputs{1}.processFcns = {'processpca'};
%             MLP_net.outputs{end}.processFcns = {};
MLP_net.layers{1:(end-1)}.transferFcn = 'poslin'; % leakyrelu, tansig, poslin, purelin
MLP_net.divideParam.trainRatio = 80/100;
MLP_net.divideParam.valRatio = 10/100;
MLP_net.divideParam.testRatio = 10/100;

%% Setup
y_star = targets;
% ys = zeros(N_iters, N_data);
ys = [];
perf = zeros(N_iters, 1);
cons = zeros(N_iters, 1);

alpha = 0.9;
N_inp = N_features; [N_out, N_neu] = size(MLP_net.LW{2,1});
solver = @fmincon;
problem.solver = func2str(solver);
problem.options = optimoptions(solver);
problem.options.Algorithm = 'interior-point'; % interior-point, sqp, active-set
problem.options.Display = 'iter';
problem.options.UseParallel = true;
problem.options.MaxFunctionEvaluations = 5e5; % 1e6 takes 30 mins
problem.options.MaxIterations = 5; % 1 min = 9 iters
%     problem.options.ConstraintTolerance = 1e-10; problem.options.OptimalityTolerance = 1e-10;

%% Algorithm
tic;

% PRETRAINING
% MLP_net = init(MLP_net);
% [MLP_net, tr] = train(MLP_net, inputs, targets, 'useParallel', 'yes'); %
% y = MLP_net(inputs);
y = y_star;
w_star = get_weights(MLP_net);
w_y = w_star;

for i=1:N_iters
    %% CONSTRAINTS / MASTER STEP
    MLP_net = init(MLP_net);
    MLP_net = train(MLP_net, inputs, (1-alpha)*y_star + alpha*y, 'useParallel', 'yes'); % 'useParallel', 'yes'
    w0 = get_weights(MLP_net);
%         w0 = (1-alpha)*w_star + alpha*w_y;

    N_w = length(w0);
    del = 5e-2;
%             problem.x0 = w0 + del * rand(N_w, 1) .* abs(w0);
    problem.x0 = w0;
%         problem.lb = w0 - del * abs(w0);
%         problem.ub = w0 + del * abs(w0);
    problem.objective = @(w) weight_fun(w, w0);
%     problem.nonlcon = @(w) zero_cons(w, MLP_net, N_inp, N_neu, N_out);
    problem.nonlcon = @(w) stability_error(w, MLP_net, N_inp, N_neu, N_out, t, ...
        X_ref0', WK_R, WK_L, INSECT, N_single, N_iters, N_per_iter, N_dang,...
        get_args, param_type, m, Weights);
    [w, fval, exitflag, output] = solver(problem); 
    MLP_net = update_weights(w, MLP_net, N_inp, N_neu, N_out);
    z = MLP_net(inputs);
    
    %% TRAINING / LEARNING
    MLP_net = init(MLP_net);
    [MLP_net, tr] = train(MLP_net, inputs, z, 'useParallel', 'yes'); %
    y = MLP_net(inputs);
    w_y = get_weights(MLP_net);
%     ys(i, :) = y;
%     ys = [ys ; y];
    cons(i) = norm(MLP_net(zeros(N_features,1)));
    perf(i) = perform(MLP_net, y_star, y);
end

time_taken = toc;

%% Performance and results
% tstPerform = perform(MLP_net, targets(:, tr.testInd), MLP_net(inputs(:, tr.testInd)));
plotregression(targets, MLP_net(inputs));
plotregression(targets(:, tr.testInd), MLP_net(inputs(:, tr.testInd)));

allvars = whos;
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
save(filename, allvars(tosave).name)
evalin('base',['load ' filename]);

%%
function [c, ceq] = stability_error(w, MLP_net, N_inp, N_neu, N_out, ...
    t, X0, WK_R, WK_L, INSECT, N_single, N_iters, N_per_iter, N_p, ...
    get_args, param_type, m, Weights)
%%
c = [];
MLP_net = update_weights(w, MLP_net, N_inp, N_neu, N_out);
param = MLP_net(zeros(N_inp, 1));

X = zeros(1+N_single, 18);

idx_con = 1:(1+N_single);
varargin = get_args(t(idx_con(1)), zeros(N_p, 1));
for con=1:m
    param_idx = (1+(con-1)*N_p):(con*N_p);
    idx = idx_con((1+(con-1)*N_per_iter):(1+con*N_per_iter));
    param_m = param(param_idx);

    [~, X(idx,:)] = ode45(@(t,X) eom_param(INSECT, WK_R, WK_L, t, X, param_m, param_type, varargin{:}), ...
    t(idx), X0, odeset('AbsTol',1e-6,'RelTol',1e-6));

    X0 = X(idx(end),:)';
    dang0 = param_type(param_m, t(idx(end)), varargin{:});
    varargin = get_args(t(idx(end)), dang0);
end

ceq = sum((Weights.OutputVariables .* (X(end, :) - X0')).^2);
% ceq = X(end, :) - X0';

end

function J = weight_fun(w, w0)
% c1 = 1; c2 = 1;
% J = c1 * norm(w - w0) + c2 * norm(MLP_net(zeros(N_inp, 1)));
J = norm(w - w0);
end

function [c, ceq] = zero_cons(w, MLP_net, N_inp, N_neu, N_out)
c = [];

MLP_net = update_weights(w, MLP_net, N_inp, N_neu, N_out);
ceq = MLP_net(zeros(N_inp, 1));
% ceq = norm(MLP_net(zeros(N_inp, 1)));
end

function MLP_net = update_weights(w, MLP_net, N_inp, N_neu, N_out)
%%
MLP_net.IW{1} = reshape(w(1:N_neu*N_inp), N_neu, N_inp);
idx = N_neu*N_inp;
% MLP_net.IW{2} = reshape(w((idx+1):(idx+N_out*N_inp)), N_out, N_inp);
% idx = idx + N_out*N_inp;
MLP_net.LW{2,1} = reshape(w((idx+1):(idx+N_out*N_neu)), N_out, N_neu);
idx = idx + N_out*N_neu;
MLP_net.b{1} = reshape(w((idx+1):(idx+N_neu)), N_neu, 1);
idx = idx+N_neu;
MLP_net.b{2} = reshape(w((idx+1):(idx+N_out)), N_out, 1);
idx = idx+N_out;
end

function w = get_weights(MLP_net)
w = [reshape(MLP_net.IW{1}, numel(MLP_net.IW{1}), 1);
%     reshape(MLP_net.IW{2}, numel(MLP_net.IW{2}), 1);
    reshape(MLP_net.LW{2,1}, numel(MLP_net.LW{2,1}), 1);
    MLP_net.b{1};
    MLP_net.b{2}];
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
