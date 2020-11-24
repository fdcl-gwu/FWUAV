function sim_QS_xR_hover_control_opt
% simulate the thorax (x) trajectory along with a controller for
% given thorax attiude, wing kinematics, abdomen attitude.

evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
des = load('sim_QS_xR_hover.mat', 'INSECT', 't', 'x', 'x_dot','R','W',...
    'x0', 'x_dot0', 'R0', 'W0', 'WK', 'X0');

filename = 'sim_QS_xR_hover_control_opt';
INSECT = des.INSECT;
WK = des.WK;
des.t = linspace(0,3/WK.f,6001)';
[~, X_ref]=ode45(@(t,X) eom_QS_xR(INSECT, WK, WK, t,X), des.t, des.X0', odeset('AbsTol',1e-6,'RelTol',1e-6));
des.x = X_ref(:, 1:3)'; des.x_dot = X_ref(:, 13:15)';
des.R = reshape(X_ref(:, 4:12)', 3, 3, []); des.W = X_ref(:, 16:18)';

des.x_fit = cell(3, 1); des.x_dot_fit = cell(3, 1);
des.R_fit = cell(3, 3); des.W_fit = cell(3, 1);
interp = 'fourier8'; % 'fourier8', 'cubicinterp'
for i=1:3
    des.x_fit{i} = fit(des.t, des.x(i, :)', interp);
    des.x_dot_fit{i} = fit(des.t, des.x_dot(i, :)', interp);
    des.W_fit{i} = fit(des.t, des.W(i, :)', interp);
    for j=1:3
        des.R_fit{i,j} = fit(des.t, squeeze(des.R(i, j, :)), interp);
    end
end

% load('sim_QS_xR_hover.mat', 'solutions');
% WK_arr = solutions(1).X; % Controllable : (2), Uncontrollable : (1)
% [WK, des.x_dot0, des.R0, des.W0] = get_WK0(WK, WK_arr);
% des.X0=[des.x0; reshape(des.R0,9,1); des.x_dot0; des.W0];

%% Optimization
rng default;
N_periods = 3; % 1, 3
N_iters = 10; % 2, 4, 10 per period
p = 2*N_iters; % Prediction Horizon multiplier of iter
m = N_iters; % Control Horizon multiplier of iter
N_single = 100; % per time period
N_per_iter = N_single / N_iters;
N = N_single*(N_periods-1) + 1 + N_per_iter*p;
err_bound = 1e-4; % Convergence criterion is a f(N, N_periods, WK)
T = (N_periods - 1 + p/N_iters)/WK.f;
t = linspace(0,T,N);
bound_param = 0.1; % Parameter bound; Use 0.25?

X_ref = zeros(N, 18);
X_ref0 = des.X0';
idx = 1:(1+N_single);
[~, X_ref(idx, :)]=ode45(@(t,X) eom_QS_xR(INSECT, WK, WK, t,X), ...
    t(idx), X_ref0, odeset('AbsTol',1e-6,'RelTol',1e-6));
idx = (1+N_single):N;
X_ref(idx, :) = X_ref(mod(idx-1, N_single)+1, :);

idx_con = 1:(1+N_single*N_periods);
des.x_fit_t = X_ref(idx_con, 1:3)'; des.x_dot_fit_t = X_ref(idx_con, 13:15)';
des.R_fit_t = reshape(X_ref(idx_con, 4:12)', 3, 3, []); des.W_fit_t = X_ref(idx_con, 16:18)';

% Weights.OutputVariables = [5*ones(1,3), 1*ones(1,9), 2.5*ones(1,3), 1*ones(1,3)]; % 5/15, 1/2.5, 2.5/5, 1
Weights.OutputVariables = [15*ones(1,3), 1.5*ones(1,9), 5*ones(1,3), 1*ones(1,3)]; % 5/15, 1/2.5, 2.5/5, 1
% Can try
% Weights.PredictionHorizon = logspace(-1, 1, p); % p elements
% Weights.PredictionHorizon = logspace(0, 1, p); % p elements
Weights.PredictionHorizon = logspace(-0.5, 1, p); % p elements
Weights.PredictionHorizon = Weights.PredictionHorizon / sum(Weights.PredictionHorizon);
% To multiply the perturbations
dx_max = 0.2*max(max(abs(X_ref(:, 1:3)), [], 1)); dtheta_max = 2.86*pi/180;
dx_dot_max = 0.05 * max(max(abs(X_ref(:, 13:15)), [], 1)); domega_max = 0.05 * max(max(abs(X_ref(:, 16:18)), [], 1));
Weights.PerturbVariables = [dx_max*ones(1,3), dtheta_max*ones(1,3), dx_dot_max*ones(1,3), domega_max*ones(1,3)];

WK_R = WK; WK_L = WK;

% FMINCON
N_dang = 6;
dang0 = zeros(N_dang, 1);

% param_type = @constant_func;
% get_args = @constant_args;
% lb = -bound_param * ones(N_dang, 1);
% ub = bound_param * ones(N_dang, 1);

param_type = @linear_func;
get_args = @linear_args;
lb = -2*bound_param*WK.f*N_iters * ones(N_dang, 1)/4;
ub = 2*bound_param*WK.f*N_iters * ones(N_dang, 1)/4;

problem.x0 = repmat(dang0, p, 1);
problem.lb = repmat(lb, p, 1);
problem.ub = repmat(ub, p, 1);
problem.x0 = 1e-2 * rand(N_dang*p, 1);
% Equality constraints work only for special values of p, m
I_66 = eye(6); O_66 = zeros(6);
problem.Aeq = [repmat(I_66, 1, N_iters), repmat(O_66, 1, N_iters);
               repmat(O_66, 1, N_iters), repmat(I_66, 1, N_iters)];
problem.beq = zeros(12, 1);
solver = @fmincon; % fmincon, patternsearch
problem.solver = func2str(solver);
problem.options = optimoptions(solver);
switch problem.solver
    case 'fmincon'
        problem.options.Algorithm = 'sqp'; % interior-point, sqp, fastest : active-set
%         problem.options.ScaleProblem = true;
    case 'patternsearch'
        problem.options.AccelerateMesh = true;
        problem.options.PollMethod = 'GSSPositiveBasisNp1';
        problem.options.UseCompletePoll = true;
%         problem.options.SearchFcn = {@searchneldermead,1,optimset('MaxFunEvals',100)};
end
% if p*N_dang*100 > 4000
%     problem.options.MaxFunctionEvaluations = 10000;
% end
problem.options.Display = 'iter';
problem.options.UseParallel = true;

%% Simulation
simulation_type = 'optimized'; % 'single', 'monte_carlo', 'optimized'
switch simulation_type
    case 'single'
%%
    e1 = [1 0 0]'; e2 = [0 1 0]'; e3 = [0 0 1]';
    des_X0 = X_ref(1, :);
    load('sim_QS_xR_hover_control_opt_mc.mat', 'dX');
    dX0 = dX(268, :);
    X0 = [des_X0(1:3)+dX0(1:3),...
            reshape(reshape(des_X0(4:12), 3, 3)*expmhat(dX0(6)*e3)*expmhat(dX0(5)*e2)*expmhat(dX0(4)*e1),1,9),...
            des_X0(13:18) + dX0(7:12)]';
    load('sim_QS_xR_hover_control_opt_net', 'control_net');
%     control_net = 'none';
    tic;
    [cost, opt_param, opt_fval, opt_iter, opt_firstorder, dang,...
        X, x, x_dot, R, Q_R, Q_L, Q_A, theta_A, W, W_R, W_R_dot, W_L, ...
        W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, f_tau, tau, Euler_R, ...
        Euler_R_dot] = simulate_opt_control(t, X0, X_ref, WK_R, WK_L, INSECT, ...
        N, N_periods, N_single, N_iters, N_per_iter, problem, solver, Weights, ...
        get_args, param_type, p, m, control_net);
    time_taken = toc;

    plot(0:1:N_periods*N_iters, cost);

    case 'monte_carlo'
%%
    load_mc_data = true; % To add similar data to the old values
%     % Longitudinal velocity perturbations only
%     Weights.PerturbVariables([1:6,8,10:12]) = 0;
%     Weights.PerturbVariables = 20 * Weights.PerturbVariables;
    N_sims = 100; % 48 take 1 day
%     problem.options.Display = 'none';
    old_seed = 'shuffle'; % default, shuffle
    if load_mc_data
        mc_old = load('sim_QS_xR_hover_control_opt_mc.mat', 'X0_pert', 'dX', ...
            'cost', 'opt_param', 'opt_fval', 'opt_iter', 'opt_firstorder', ...
            'opt_time', 'new_seed', 'N_sims', 'time_taken');
        old_seed = mc_old.new_seed;
%         arr_idx = 1:(3*mc_old.N_sims);
%         dX_idx = intersect(1 + mod(arr_idx(mc_old.opt_iter < 7)-1, mc_old.N_sims), 1 + mod(arr_idx(mc_old.cost(:,end) > mc_old.cost(:,1))-1, mc_old.N_sims));
%         dX0 = mc_old.dX(dX_idx, :);
%         N_sims = length(dX_idx);
    end
    t_init = tic;
    load('sim_QS_xR_hover_control_opt_mc', 'control_net', 'tr', 'inputs', 'targets');
%     control_net = 'none';
    [X0_pert, dX, cost, opt_param, opt_fval, opt_iter, opt_firstorder, ...
    opt_time, new_seed, N_sims] = monte_carlo(N_sims, t, X_ref, WK_R, WK_L, INSECT, ...
            N, N_periods, N_single, N_iters, N_per_iter, problem, solver, Weights, ...
            get_args, param_type, p, m, old_seed, control_net);
    time_taken = toc(t_init);
    if load_mc_data
        time_taken = time_taken + mc_old.time_taken;
        mc_new = struct('N_sims', N_sims, 'X0_pert', X0_pert, 'dX', dX, ...
            'cost', cost, 'opt_param', opt_param, 'opt_fval', opt_fval, ...
            'opt_iter', opt_iter, 'opt_firstorder', opt_firstorder, 'opt_time', opt_time);
        N_sims = mc_new.N_sims + mc_old.N_sims;
        X0_pert = zeros(N_sims*N_periods, size(X0_pert, 2), size(X0_pert, 3));
        dX = [mc_old.dX; dX];
        cost = zeros(N_sims*N_periods, size(cost, 2));
        opt_param = zeros(N_sims*N_periods, size(opt_param, 2), size(opt_param, 3));
        opt_fval = zeros(N_sims*N_periods, size(opt_fval, 2));
        opt_iter = zeros(N_sims*N_periods, size(opt_iter, 2));
        opt_firstorder = zeros(N_sims*N_periods, size(opt_firstorder, 2));
        opt_time = [mc_old.opt_time; opt_time];
        for mc_var = {'X0_pert', 'cost', 'opt_param', 'opt_fval', 'opt_iter', 'opt_firstorder'}
            for i=1:N_periods    
                eval([mc_var{1} '((1+(i-1)*N_sims):((i-1)*N_sims + mc_old.N_sims), :, :) = mc_old.' ...
                    mc_var{1} '((1+(i-1)*mc_old.N_sims):(i*mc_old.N_sims), :, :);']);
                eval([mc_var{1} '((1+(i-1)*N_sims + mc_old.N_sims):(i*N_sims), :, :) = mc_new.' ...
                    mc_var{1} '((1+(i-1)*mc_new.N_sims):(i*mc_new.N_sims), :, :);']);
%                 eval([mc_var{1} '((i-1)*N_sims+dX_idx, :, :) = mc_new.' mc_var{1} '((1+(i-1)*mc_new.N_sims):(i*mc_new.N_sims), :, :);']);
            end
        end
        clear('mc_old');
    end
    if all(cost(:, end) < cost(:, 1))
        fprintf('Cost has decreased in each period\n');
    else
        fprintf('Cost has increased for some initial condition\n');
    end
    
    case 'optimized'
%%
    mc = load('sim_QS_xR_hover_control_opt_mc.mat', 'N_sims', 'X0_pert', 'cost', 'opt_param');
    use_control_net = 'period'; % 'none', 'iter', 'period'
    mc.X0_pert_flat = squeeze(mc.X0_pert(:, 1, :)); mc.opt_param_flat = reshape(permute(mc.opt_param, [1, 3, 2]), 3*mc.N_sims, 6*N_iters);
    N_periods = 3;
    control_net = 0;
    % Regression techniques:
    % 1. mvregress (simple linear) : doesn't work
    % 2. newgrnn (simple neural) : not good enough??
    % 4. Normalization : data_norm = (1 ./ Weights.PerturbVariables)' .* data_x;
    % 5. PCA : [coeff, score, latent] = pca(data_norm');
    % TODO : Preprocessing Data, Shallow NN Time
    % Series, cvpartition (k-fold), fitensemble, Bagging/boosting
    switch use_control_net
        case 'iter'
        mc.opt_dang = zeros(size(mc.opt_param)); t_iter = (1/WK.f)/N_iters;
        mc.err = zeros(3*mc.N_sims, N_iters, 12);
        for i=1:(3*mc.N_sims)
            param = squeeze(mc.opt_param(i, :, :));
            mc.opt_dang(i, :, :) = (cumsum(param, 1) - 0.5 * param) * t_iter;
            for j=1:N_iters
                mc.err(i, j, :) = get_error(squeeze(mc.X0_pert(i, j, :)), des.X0);
            end
        end
        control_net = newgrnn((1 ./ Weights.PerturbVariables)' .* reshape(mc.err, 3*mc.N_sims*N_iters, [])', reshape(mc.opt_dang, 3*mc.N_sims*N_iters, 6)', 0.1);
        case 'period'
        mc.err = zeros(3*mc.N_sims, 12);
        for i=1:(3*mc.N_sims)
            mc.err(i, :) = get_error(mc.X0_pert_flat(i, :)', des.X0)';
        end
        inputs = (1 ./ Weights.PerturbVariables)' .* mc.err'; targets = mc.opt_param_flat';
        train_new_net = true;
        if train_new_net
            control_net = cascadeforwardnet([36]); % 24, 36
%             control_net = fitnet([6, 36]); % 24, 36
%             control_net.inputConnect(:) = 1;
% %             control_net.layerConnect(3,1) = 1;
            control_net.biasConnect = [0, 0]';
            control_net = configure(control_net, inputs, targets);
            control_net.trainFcn = 'trainbr';% train{lm,br,bfg}
            control_net.trainParam.epochs = 30; % 30
            control_net.performParam.normalization = 'standard';
            control_net.inputs{1}.processFcns = {'mapminmax', 'processpca'};
    %         control_net.inputWeights{1,1}.weightFcn
    %         control_net.layerWeights{2,1}.weightFcn
    %         control_net.layers{1}.transferFcn = 'poslin';
    %         control_net.layers{1}.netInputFcn = 'netprod';
            control_net.divideParam.trainRatio = 80/100;
            control_net.divideParam.valRatio = 10/100;
            control_net.divideParam.testRatio = 10/100;
    %         control_net.trainParam.mu = 1; control_net.trainParam.mu_dec = 0.8; control_net.trainParam.mu_inc = 1.5;
            control_net = init(control_net);
            [control_net, tr] = train(control_net, inputs, targets);
            tstPerform = perform(control_net, targets(:, tr.testInd), control_net(inputs(:, tr.testInd)));
        else
            load('sim_QS_xR_hover_control_opt.mat', 'control_net', 'tr');
%             control_net = newgrnn(inputs, targets, 0.5);
        end

%         beta = mvregress(mc.err, mc.opt_param_flat, 'algorithm', 'cwls');
%         control_net = newgrnn(mc.X0_pert_flat', mc.opt_param_flat', 0.1);
    end
    arr_idx = 1:(3*mc.N_sims);
    inc_idx = 1 + mod(arr_idx(mc.cost(:,end) > mc.cost(:,1)) - 1, mc.N_sims);
    opt_idx = inc_idx(1); % 3 for mc
    t = t(1:(1+N_periods*N_single));
    idx_con = 1:(1+N_single*N_periods);
    X0 = squeeze(mc.X0_pert(opt_idx, 1, :));
    if true
        rng('shuffle');
        dx = 2*rand(1,3)-1; dx = rand(1) * dx / norm(dx);
        dtheta = 2*rand(1,3)-1; dtheta = rand(1) * dtheta / norm(dtheta);
        dx_dot = 2*rand(1,3)-1; dx_dot = rand(1) * dx_dot / norm(dx_dot);
        domega = 2*rand(1,3)-1; domega = rand(1) * domega / norm(domega);
        dXi = Weights.PerturbVariables .* [dx, dtheta, dx_dot, domega];
        e1 = [1 0 0]'; e2 = [0 1 0]'; e3 = [0 0 1]';
        des_X0 = X_ref(1, :);
        des_R0 = reshape(des_X0(4:12), 3, 3);
        X0 = [des_X0(1:3)+dXi(1:3),...
            reshape(des_R0*expmhat(dXi(6)*e3)*expmhat(dXi(5)*e2)*expmhat(dXi(4)*e1),1,9),...
            des_X0(13:18) + dXi(7:12)]';
    end
    cost = zeros(1, 1+N_periods*N_iters);
    opt_param = zeros(N_dang*m, N_periods*N_iters/m);
    for i=1:N_periods
        cost((1+(i-1)*N_iters):(1+i*N_iters)) = mc.cost(opt_idx + (i-1)*mc.N_sims, :);
        if strcmp(use_control_net, 'none')
            opt_param(:, i) = reshape(squeeze(mc.opt_param(opt_idx + (i-1)*mc.N_sims, :, :))', 60, 1);
        end
    end
    tic;
    [cost, opt_param, dang, X, x, x_dot, R, Q_R, Q_L, Q_A, theta_A, W, W_R, W_R_dot, W_L, ...
    W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, f_tau, tau, Euler_R, ...
    Euler_R_dot] = simulate_control(t, X0, X_ref, WK_R, WK_L, INSECT, ...
    N_periods, N_single, N_iters, N_per_iter, problem, Weights, ...
    get_args, param_type, p, m, opt_param, use_control_net, control_net, des.X0);
    time_taken = toc;

    plot(0:1:N_periods*N_iters, cost);
end

%%
t = t(idx_con);
N = 1+N_periods*N_single;
% Get a list of all variables
allvars = whos;
% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics object
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
% Pass these variable names to save
save(filename, allvars(tosave).name)
evalin('base',['load ' filename]);

end

function [X0_pert, dX, cost, opt_param, opt_fval, opt_iter, opt_firstorder, ...
    opt_time, new_seed, N_sims] = monte_carlo(N_sims, t, X_ref, WK_R, WK_L, INSECT, ...
        N, N_periods, N_single, N_iters, N_per_iter, problem, solver, Weights, ...
        get_args, param_type, p, m, old_seed, control_net, varargin)
%%
e1 = [1 0 0]'; e2 = [0 1 0]'; e3 = [0 0 1]';
opt_time = zeros(N_sims, 1);
dX = zeros(N_sims, 12);
cost = zeros(N_sims, N_periods, 1+N_iters);
N_p = length(problem.x0)/p;
opt_param = zeros(N_sims, N_periods, N_iters, N_p); opt_fval = zeros(N_sims, N_periods, 1);
opt_iter = zeros(N_sims, N_periods, 1); opt_firstorder = zeros(N_sims, N_periods, 1);
X0_pert = zeros(N_sims, N_periods, 1+N_iters, 18);
idx_cost = zeros(N_periods, 1+N_iters); idx_X0 = zeros(N_periods, 1+N_iters);
for period=1:3
    idx_cost(period,:) = (1+(period-1)*N_iters):(1+period*N_iters);
    idx_X0(period,:,:) = (1+(period-1)*N_single):N_per_iter:(1+period*N_single);
end
des_X0 = X_ref(1, :);
des_R0 = reshape(des_X0(4:12), 3, 3);

% par_pool = gcp;
% nWorkers = par_pool.NumWorkers;
% ticBytes(par_pool);
% parfor n=1:nWorkers
%     rng(n);
% end
pause(1);

rng(old_seed);
try
    dX = varargin{1};
catch
    for i = 1:N_sims
        dx = 2*rand(1,3)-1; dx = rand(1) * dx / norm(dx);
        dtheta = 2*rand(1,3)-1; dtheta = rand(1) * dtheta / norm(dtheta);
        dx_dot = 2*rand(1,3)-1; dx_dot = rand(1) * dx_dot / norm(dx_dot);
        domega = 2*rand(1,3)-1; domega = rand(1) * domega / norm(domega);
        dX(i,:) = Weights.PerturbVariables .* [dx, dtheta, dx_dot, domega];
    end
end
new_seed = rng;

stopbutton = uicontrol('style', 'toggle', 'String', 'STOP');
set(stopbutton, 'Value', 0);
drawnow();%give time for callback
for i = 1:N_sims
    % perturbation
    dXi = dX(i,:);
    X0 = [des_X0(1:3)+dXi(1:3),...
        reshape(des_R0*expmhat(dXi(6)*e3)*expmhat(dXi(5)*e2)*expmhat(dXi(4)*e1),1,9),...
        des_X0(13:18) + dXi(7:12)]';
    tic;
    [cost_arr, opt_param_arr, opt_fval_arr, opt_iter_arr, opt_firstorder_arr, ~, X] = ...
        simulate_opt_control(t, X0, X_ref, WK_R, WK_L, INSECT, ...
        N, N_periods, N_single, N_iters, N_per_iter, problem, solver, Weights, ...
        get_args, param_type, p, m, control_net);
    opt_time(i) = toc;
    cost(i,:,:) = cost_arr(idx_cost);
    opt_param(i,:,:,:) = permute(reshape(opt_param_arr, N_p, N_iters, N_periods), [3, 2, 1]);
    opt_fval(i,:,:) = opt_fval_arr';
    opt_iter(i,:,:) = opt_iter_arr'; opt_firstorder(i,:,:) = opt_firstorder_arr';
    X0_pert(i,:,:,:) = reshape(X(idx_X0, :), N_periods, 1+N_iters, 18);
    fprintf('Completed simulation %d\n', i);
    drawnow();
    quitthis = get(stopbutton, 'Value');
    if quitthis
        N_sims = i;
        break;
    end
end

dX = dX(1:N_sims, :); opt_time = opt_time(1:N_sims, :);
% row (i, i+N_sims, i+2*N_sims ...) correspond to same dX
cost = reshape(cost(1:N_sims,:,:), N_sims*N_periods, 1+N_iters);
opt_param = reshape(opt_param(1:N_sims,:,:,:), N_sims*N_periods, N_iters, N_p);
opt_fval = reshape(opt_fval(1:N_sims,:,:), N_sims*N_periods, 1);
opt_iter = reshape(opt_iter(1:N_sims,:,:), N_sims*N_periods, 1);
opt_firstorder = reshape(opt_firstorder(1:N_sims,:,:), N_sims*N_periods, 1);
X0_pert = reshape(X0_pert(1:N_sims,:,:,:), N_sims*N_periods, 1+N_iters, 18);

% tocBytes(par_pool);

end

function [cost, opt_param, opt_fval, opt_iter, opt_firstorder, dang, ...
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
        opt_complete = false; opt_count = 0;
        while ~opt_complete
            [param, fval, exitflag, output] = solver(problem);
            opt_count = opt_count + 1;
%             if ((output.iterations >= 20 || period > 1) && output.firstorderopt < 1e3) || opt_count >= 3
            if (((period == 1 && output.iterations >= 15) || (period > 1 && output.iterations >= 7)) && (exitflag ~= -2)) || opt_count >= 3 % period > 1 || 
                opt_complete = true;
                fprintf('Exit flag is %d\n', exitflag);
            else
                if period == 1
                    problem.x0 = problem.x0 + 1e-2 * (2*rand(N_p*p, 1)-1) .* problem.ub;
                else
                    problem.x0((1+N_p*(p-m)):end) = problem.x0((1+N_p*(p-m)):end) + ...
                        5e-3 * (2*rand(N_p*m, 1)-1) .* problem.ub((1+N_p*(p-m)):end);
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
        for con=1:m
            param_idx = (1+(con-1)*N_p):(con*N_p);
            param_m = param(param_idx);
            idx = idx_con((1+(con-1)*N_per_iter):(1+con*N_per_iter));

            [~, X(idx,:)] = ode45(@(t,X) eom_param(INSECT, WK_R, WK_L, t, X, param_m, param_type, varargin{:}), ...
            t(idx), X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
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

function [cost, opt_param, dang, X, x, x_dot, R, Q_R, Q_L, Q_A, theta_A, W, W_R, W_R_dot, W_L, ...
    W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, f_tau, tau, Euler_R, ...
    Euler_R_dot] = simulate_control(t, X0, X_ref, WK_R, WK_L, INSECT, ...
    N_periods, N_single, N_iters, N_per_iter, problem, Weights, ...
    get_args, param_type, p, m, opt_param, use_control_net, control_net, des_X0)
%%
X = zeros(1+N_single*N_periods, 18);
N_p = length(problem.x0) / p;
dang = zeros(N_p, 1+N_single*N_periods);
cost = zeros(1, 1+N_periods*N_iters);
cost(1) = sqrt(sum((Weights.OutputVariables .* (X0' - X_ref(1, :))).^2));
if ~strcmp(use_control_net, 'none')
    opt_param = zeros(N_p*m, N_periods*N_iters/m);
end

for period=1:N_periods
    for iter=1:m:N_iters
        idx_con = (1+(period-1)*N_single+(iter-1)*N_per_iter):(1+(period-1)*N_single+(iter+m-1)*N_per_iter);
        varargin = get_args(t(idx_con(1)), dang(:,idx_con(1)));
        if strcmp(use_control_net, 'none')
            param = opt_param(:, 1+((period-1)*N_iters+(iter-1))/m);
        else
            if strcmp(use_control_net, 'period')
            param = sim(control_net, (1 ./ Weights.PerturbVariables)' .* get_error(X0, des_X0));
            opt_param(:,1+((period-1)*N_iters+(iter-1))/m) = param;
            end
        end
        for con=1:m
            param_idx = (1+(con-1)*N_p):(con*N_p);
            idx = idx_con((1+(con-1)*N_per_iter):(1+con*N_per_iter));
            if strcmp(use_control_net, 'iter')
                opt_dang = sim(control_net, (1 ./ Weights.PerturbVariables)' .* get_error(X0, des_X0));
                param_m = 2*(opt_dang - dang(:, idx(1))) / (t(idx(end)) - t(idx(1)));
                opt_param(param_idx, 1+((period-1)*N_iters+(iter-1))/m) = param_m;
            else
                param_m = param(param_idx);
            end

            [~, X(idx,:)] = ode45(@(t,X) eom_param(INSECT, WK_R, WK_L, t, X, param_m, param_type, varargin{:}), ...
            t(idx), X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
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
    end
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
J = 0;
X = zeros(size(X_ref));
int_dx = zeros(size(X_ref(:, 1:3))); dt = t(2) - t(1);
for pred=1:p
    param_idx = (1+(pred-1)*N_p):(pred*N_p);
    param_p = param(param_idx);
    idx = (1+(pred-1)*N_per_iter):(1+pred*N_per_iter);

    [~, X_temp] = ode45(@(t,X) eom_param(INSECT, WK_R, WK_L, t, X, param_p, param_type, varargin{:}), ...
    t(idx), X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
    X_len = size(X_temp,1);
    if X_len ~= (N_per_iter+1)
%         J = J + sqrt(sum((Weights.OutputVariables .* (X_temp(end,:) - X_ref(idx(1)+X_len-1,:))).^2));
%         J = J*(exp(N_per_iter+1-X_len)-1);
        J = 1e10;
        return;
    else
        X(idx,:) = X_temp;
    end
    int_temp = cumsum([int_dx(idx(1), :); (X(idx, 1:3) - X_ref(idx, 1:3))*dt], 1);
    int_dx(idx, :) = int_temp(1:end-1, :);
    X0 = X(idx(end),:)';
    dang0 = param_type(param_p, t(idx(end)), varargin{:});
    varargin = get_args(t(idx(end)), dang0);
    J = J + Weights.PredictionHorizon(pred) * ...
        sqrt(sum((Weights.OutputVariables .* (X(idx(end),:) - X_ref(idx(end),:))).^2)...
        + sum((Weights.OutputVariables(1:3) .* int_dx(idx(end), :)).^2));
end
% J = sqrt(mean(sum((Weights.OutputVariables .* (X - X_ref)).^2, 2))) + ...
%     100*sqrt(sum((Weights.OutputVariables .* (X(end,:) - X_ref(end,:))).^2));

% [~, X] = ode45(@(t,X) eom_param(INSECT, WK_R, WK_L, t, X, param, param_type, varargin{:}), ...
%     t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6)); % [t(1), t(end)] instead of [t]

end

function [c, ceq] = param_cons(param, p, N_iters)
%%
c = [];
param = reshape(param, 6, p);
ceq(1:6) = sum(param(:,1:N_iters), 2);
ceq(7:12) = sum(param(:,(N_iters+1):p), 2);
end

function dXdt = eom_param(INSECT, WK_R, WK_L, t, X, param, param_type, varargin)
%%
dang = param_type(param, t, varargin{:});
[WK_R, WK_L] = get_WK(WK_R, WK_L, dang);
dXdt = eom_QS_xR(INSECT, WK_R, WK_L, t, X);
end

function dang = constant_func(param, t, varargin)
dang = param;
end

function dang = linear_func(param, t, varargin)
t0 = varargin{1};
dang0 = varargin{2};
dang = dang0 + param*(t - t0);
end

function args = constant_args(t0, dang0)
args = {};
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

function MPC
%%
nx = 19;
ny = 19;
nu = 6;
nlobj = nlmpc(nx, ny, nu);

nlobj.Model.IsContinuousTime = true;
nlobj.Model.StateFcn = @(X, U, INSECT, WK_R, WK_L, X_T, Weights) eom_control(X, U, INSECT, WK_R, WK_L);
nlobj.Model.NumberOfParameters = 5;

% Ts = (T / N_periods) / 8; p = 51; m = 10;
% Ts = t(2) - t(1); p = N; m = 1;
p = 100; m = 1;
% m = [p/2, p/2];
% m = [p];
Ts = T/p;
nlobj.Ts = Ts; % sample time
nlobj.PredictionHorizon = p;
nlobj.ControlHorizon = m;

nlobj.Optimization.CustomCostFcn = @(X,U,e,data,INSECT,WK_R,WK_L,X_T,Weights) cost_fun(X,U,e,data,INSECT,WK_R,WK_L,X_T,Weights);
nlobj.Optimization.ReplaceStandardCost = true;
options = optimoptions(@fmincon,'Algorithm','interior-point','Display','none',...
    'SpecifyConstraintGradient', true, 'SpecifyObjectiveGradient', true,...
    'UseParallel', true);
nlobj.Optimization.SolverOptions = options;
nlobj.Weights.OutputVariables = [1, 1, 1, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.5, 0.5, 0.5, 0];
% nlobj.Weights.ManipulatedVariables = 0.05 * ones(1, nu);
% nlobj.Weights.ManipulatedVariablesRate = 1e10 * ones(1, nu);
Weights = nlobj.Weights;
Weights.Time = logspace(-1, 1, p)';

for i = 1:nu
    nlobj.MV(i).Min = -bound_param;
    nlobj.MV(i).Max = bound_param;
end

t0 = 0;
X0 = [X0; t0];
U0 = zeros(nu, 1);
X_ref = [des.x_fit_t' squeeze(reshape(des.R_fit_t, 1, 9, N))' des.x_dot_fit_t' des.W_fit_t' t'];
X_T = X_ref(end, :);
params = {INSECT, WK, WK, X_T, Weights};
validateFcns(nlobj, X0, U0, [], params);

nloptions = nlmpcmoveopt;
nloptions.MVTarget = zeros(1, 6);
nloptions.Parameters = params;
mv = nloptions.MVTarget;

tic;
% [uk,nloptions,info] = nlmpcmove(nlobj,X0,U0,X_ref(1:p, :),[],nloptions);
[uk,nloptions,info] = nlmpcmove(nlobj,X0,U0,X_ref(linspace(round(1+N/p), N, p), :),[],nloptions);
toc; % ~ 2153 seconds, p = 101, m = 1

dang = info.MVopt(1, :)';
end
