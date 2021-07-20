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
N_iters = 10;

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
perf = zeros(N_iters, 1);
perf_z = zeros(N_iters, 1);
cons = zeros(N_iters, 1);
cons_z = zeros(N_iters, 1);

alpha = 0.9;
N_inp = N_features; [N_out, N_neu] = size(MLP_net.LW{2,1});
solver = @fmincon;
problem.solver = func2str(solver);
problem.options = optimoptions(solver);
problem.options.Algorithm = 'interior-point'; % interior-point, sqp, active-set
problem.options.Display = 'iter';
problem.options.UseParallel = true;
problem.options.MaxFunctionEvaluations = 5e5; % 1e6 takes 30 mins
problem.options.MaxIterations = 10; % 5 iters
problem.options.ConstraintTolerance = 1e-10;
% problem.options.SpecifyObjectiveGradient = true;
problem.options.ObjectiveLimit = 0;
% problem.options.OptimalityTolerance = 1e-10;

%% Algorithm
tic;

% PRETRAINING
MLP_net = init(MLP_net);
% MLP_net = configure(MLP_net, inputs, targets);
[MLP_net, tr] = train(MLP_net, inputs, targets, 'useParallel', 'yes'); %
y = MLP_net(inputs);
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
    problem.lb = w0 - del * abs(w0);
    problem.ub = w0 + del * abs(w0);
%     problem.objective = @(w) weight_fun(w, w0);
    problem.nonlcon = @(w) zero_cons(w, MLP_net, N_inp, N_neu, N_out);
    
%     problem.objective = @(w) stability_obj(w, MLP_net, N_inp, N_neu, N_out, t, ...
%         X_ref0', WK_R, WK_L, INSECT, N_single, N_iters, N_per_iter, N_dang,...
%         get_args, param_type, m, Weights);
%     problem.nonlcon = @(w) stability_error(w, MLP_net, N_inp, N_neu, N_out, t, ...
%         X_ref0', WK_R, WK_L, INSECT, N_single, N_iters, N_per_iter, N_dang,...
%         get_args, param_type, m, Weights);

    % FLOQUET
%     N_period = 1;
%     N_single = 100;
%     N = 1 + N_period*N_single;
%     T=N_period/WK_R.f;
%     t=linspace(0,T,N);
%     X0 = X_ref0';
%     epsilon = 1e-6; % 1e-10 stable, 1e-12 unstable
%     n = 12; n_vars = 18; % for (x,R) control with @eom_hover_xR_control
%     problem.objective = @(w) floquet_func(w, w0, MLP_net, N_inp, N_neu, N_out, ...
%                     n, n_vars, INSECT, WK_R, X0, N, t, epsilon, ...
%                     Weights, X0, N_dang, m, N_single, N_period);
                
    n = 12; epsilon = 1e-6;
    delta0 = epsilon * Weights.PerturbVariables .* eye(n);
%     n = 1; epsilon = 1e-6;
%     delta0 = epsilon * Weights.PerturbVariables';
    problem.objective = @(w) period_error(w, MLP_net, N_inp, N_neu, N_out, ...
        t, X_ref0', WK_R, WK_L, INSECT, N_single, N_iters, N_per_iter, N_dang, ...
        get_args, param_type, m, Weights, delta0, n);

    [w, fval, exitflag, output] = solver(problem); 
    MLP_net = update_weights(w, MLP_net, N_inp, N_neu, N_out);
    z = MLP_net(inputs);
    cons_z(i) = norm(MLP_net(zeros(N_features,1)));
    
    %% TRAINING / LEARNING
    MLP_net = init(MLP_net);
    [MLP_net, tr] = train(MLP_net, inputs, z, 'useParallel', 'yes'); %
    y = MLP_net(inputs);
    w_y = get_weights(MLP_net);
%     ys(i, :) = y;
    cons(i) = norm(MLP_net(zeros(N_features,1)));
    perf(i) = perform(MLP_net, y_star, y);
    perf_z(i) = perform(MLP_net, z, y);
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
function J = period_error(w, MLP_net, N_inp, N_neu, N_out, ...
    t, X0, WK_R, WK_L, INSECT, N_single, N_iters, N_per_iter, N_p, ...
    get_args, param_type, m, Weights, delta0, n)
%%
MLP_net = update_weights(w, MLP_net, N_inp, N_neu, N_out);

idx_con = 1:(1+N_single);
err = zeros(n, 1);

for j=1:n
%     j = 1;
%     X = zeros(1+N_single, length(X0));
    dR = delta0(4:6, j);
    X_del0 = [X0(1:3)+delta0(1:3, j); reshape(reshape(X0(4:12),3,3)*expmhat(dR), 9, 1); ...
        X0(13:18) + delta0(7:12, j)];
    param = MLP_net((1 ./ Weights.PerturbVariables)' .* get_error(X_del0, X0));
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

    err(j) = err(j) + sum((Weights.OutputVariables .* (X(end, :) - X(1, :))).^2);
end

J = sum(err);

end

function J = floquet_func(w, w0, control_net, N_inp, N_neu, N_out, ...
    n, n_vars, INSECT, WK, X0, N, t, epsilon, ...
    Weights, des_X0, N_dang, m, N_single, N_period)
%% Objective function
control_net = update_weights(w, control_net, N_inp, N_neu, N_out);

[delta_mat] = sim_pert(@eom_hover_xR_control, n, n_vars, INSECT, WK, X0, N, t, epsilon, ...
    control_net, Weights, des_X0, N_dang, m, N_single, N_period);

% J = sum((Weights.PerturbVariables' .* delta_mat(:, :, end) / epsilon).^2, 'all');

B = delta_mat(:, :, 1) \ delta_mat(:, :, 1+N_single);
[~, rhos] = eig(B);
mus = log(abs(diag(rhos))) * WK.f;

idx = (mus >= 0);
% if any(idx)
%     J = log(eps + sum(exp(mus(idx)) - 1));
% else
%     J = log(eps) + max(mus)*1e3;
% end
J = exp(-(max(mus(~idx))^2)) + 1e6 * log(1 + sum(exp(mus(idx)) - 1));

if isnan(J)
    J = 1/eps;
end

end

function [delta_mat] = sim_pert(eom, n, n_vars, INSECT, WK, X0, N, t, epsilon, varargin)
%%
control_net = varargin{1};
Weights = varargin{2};
des_X0 = varargin{3};
N_p = varargin{4};
m = varargin{5};
N_single = varargin{6};
N_period = varargin{7};
N_per_iter = N_single / m;

delta0 = epsilon * Weights.PerturbVariables .* eye(n);
X0 = [X0(1:n_vars); reshape(delta0, n^2, 1);];
t0 = t(1); dang0 = zeros(N_p,1); param_mat = zeros(m*N_p, n); dang0_mat = zeros(N_p, n);
X = zeros(N, length(X0));

for period=1:N_period
    idx_con = (1+(period-1)*N_single):(1+period*N_single);
    delta_mat0 = reshape(X0((n_vars+1):(n^2+n_vars)), n, n);
    param = sim(control_net, (1 ./ Weights.PerturbVariables)' .* get_error(X0, des_X0));
    for j=1:n
        dR = delta_mat0(4:6, j);
        X_del0 = [X0(1:3)+delta_mat0(1:3, j); reshape(reshape(X0(4:12),3,3)*expmhat(dR), 9, 1); ...
            X0(13:18) + delta_mat0(7:12, j)];
        param_mat(:, j) = sim(control_net, (1 ./ Weights.PerturbVariables)' .* get_error(X_del0, des_X0));
    end
    for con=1:m
        param_idx = (1+(con-1)*N_p):(con*N_p);
        idx = idx_con((1+(con-1)*N_per_iter):(1+con*N_per_iter));
        param_m = param(param_idx);
        param_m_mat = param_mat(param_idx, :);
        
        [~, X_temp] = ode45(@(t,X) eom(n, n_vars, INSECT, WK, WK, t, X, ...
            t0, dang0, param_m, dang0_mat, param_m_mat), t(idx), X0);
        if size(X_temp,1) ~= (N_per_iter+1)
            delta_mat = 1e10 * delta0;
            return;
        else
            X(idx,:) = X_temp;
        end
        
        X0 = X(idx(end),:)';
        dang0 = linear_func(param_m, t(idx(end)), t0, dang0);
        for j=1:n
            dang0_mat(:,j) = linear_func(param_m_mat(:,j), t(idx(end)), t0, dang0_mat(:,j));
        end
        t0 = t(idx(end));
    end
end

delta_mat = reshape(X(:,n_vars+1:(n^2+n_vars))', n, n, N);

end

function [X_dot] = eom_hover_xR_control(n, n_vars, INSECT, WK_R_des, WK_L_des, t, X, varargin)
%%
t0 = varargin{1};
dang0 = varargin{2};
param = varargin{3};
dang0_mat = varargin{4};
param_mat = varargin{5};

x=X(1:3);
R=reshape(X(4:12),3,3);
x_dot=X(13:15);
W=X(16:18);

delta_mat=reshape(X(n_vars+1:(n^2+n_vars)), n, n);
X = X(1:n_vars);

%% Ideal uncontrolled values
[Euler_R, Euler_R_dot, Euler_R_ddot] = wing_kinematics(t,WK_R_des);
[Euler_L, Euler_L_dot, Euler_L_ddot] = wing_kinematics(t,WK_L_des);
[Q_R, Q_L, W_R, W_L, W_R_dot, W_L_dot] = wing_attitude(WK_R_des.beta, ...
    Euler_R, Euler_L, Euler_R_dot, Euler_L_dot, Euler_R_ddot, Euler_L_ddot);
[Q_A, W_A, W_A_dot] = abdomen_attitude(t, WK_R_des.f, WK_R_des);

%% Ideal controlled trajectory
dang = linear_func(param, t, t0, dang0);
[WK_R_con, WK_L_con] = get_WK(WK_R_des, WK_L_des, dang);
X_dot = eom_QS_xR(INSECT, WK_R_con, WK_L_con, t, X(1:18));
xi_dot = X_dot(13:18);
JJ_xR = get_cross_terms(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);

%% Controlled trajectory with perturbation
delta_mat_dot = zeros(n,n);
parfor j=1:n
    param = param_mat(:, j);
    dang = linear_func(param, t, t0, dang0_mat(:, j));
    [WK_R_del, WK_L_del] = get_WK(WK_R_des, WK_L_des, dang);
    
    dx = delta_mat(1:3, j); dR = delta_mat(4:6, j);
    dx_dot = delta_mat(7:9, j); dW = delta_mat(10:12, j);
    X_del = [x+dx; reshape(R*expmhat(dR), 9, 1); x_dot+dx_dot; W+dW];
    
    X_dot_del = eom_QS_xR(INSECT, WK_R_del, WK_L_del, t, X_del(1:18));
    xi_dot_del = X_dot_del(13:18);
    JJ_xR_del = get_cross_terms(INSECT, R*expmhat(dR), Q_R, Q_L, Q_A, ...
        x_dot+dx_dot, W+dW, W_R, W_L, W_A);

    dxi_ddot = JJ_xR \ ((JJ_xR_del*xi_dot_del(1:6) - JJ_xR*xi_dot(1:6)) -...
                        (JJ_xR_del - JJ_xR)*xi_dot(1:6));
    %
    d_mat = zeros(n, 1);
    d_mat(1:3) = dx_dot;
    d_mat(4:6) = -hat(W)*dR + dW;
    d_mat(7:12) = dxi_ddot(1:6);
    delta_mat_dot(:, j) = d_mat;
end

X_dot = [X_dot(1:n_vars); reshape(delta_mat_dot, n^2, 1);];

end

function [JJ_xR] = get_cross_terms(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A)
    JJ = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
    [JJ_11, ~, JJ_21, ~] = inertia_sub_decompose_6_9(JJ);
    C=[zeros(3,9);
        -Q_R -Q_L -Q_A];
    JJ_xR = (JJ_11 - C*JJ_21);
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
% J = c1 * norm(w - w0) + c2 * norm(MLP_net(zeros(N_inp, 1)));
J = 0.5 * (norm(w - w0))^2;

if nargout > 1
    dJ = w - w0;
end

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
