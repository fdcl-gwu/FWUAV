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

Weights.OutputVariables = [15*ones(1,3), 1*ones(1,9), 5*ones(1,3), 1*ones(1,3)]; % 15/25, 1/2.5, 5, 1
Weights.PredictionHorizon = logspace(-1, 1, p); % p elements
Weights.PredictionHorizon = Weights.PredictionHorizon / sum(Weights.PredictionHorizon);
% To multiply the perturbations
dx_max = 0.2*max(max(abs(X_ref(:, 1:3)), [], 1)); dtheta_max = 5*pi/180;
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
lb = -2*bound_param*WK.f*N_iters * ones(N_dang, 1);
ub = 2*bound_param*WK.f*N_iters * ones(N_dang, 1);

problem.x0 = repmat(dang0, p, 1);
problem.lb = repmat(lb, p, 1);
problem.ub = repmat(ub, p, 1);
problem.solver = 'fmincon';
problem.options = optimoptions(@fmincon);
problem.options.Algorithm = 'sqp'; % interior-point, sqp, fastest : active-set
if p*N_dang*100 > 4000
    problem.options.MaxFunctionEvaluations = 4000;
end
problem.options.Display = 'iter';
problem.options.UseParallel = true;

%% Simulation
simulation_type = 'single'; % 'single', 'monte_carlo'
switch simulation_type
    case 'single'
    load('sim_QS_xR_hover_control_opt_mcdata', 'X0_pert', 'cost', 'dX');
    idx_unstable = cost(:, end) > cost(:, 1);
    % X0 = X0_pert(idx_unstable, :);
    % X0 = X0(end, :)';
    dX0 = dX(idx_unstable, :);
    dX0 = dX0(end, :);
    e1 = [1 0 0]'; e2 = [0 1 0]'; e3 = [0 0 1]';
    des_X0 = X_ref(1, :);
    dX0 = Weights.PerturbVariables;
    X0 = [des_X0(1:3)+dX0(1:3),...
            reshape(reshape(des_X0(4:12), 3, 3)*expmhat(dX0(6)*e3)*expmhat(dX0(5)*e2)*expmhat(dX0(4)*e1),1,9),...
            des_X0(13:18) + dX0(7:12)]';
    tic;
    [cost, opt_param, dang, X, x, x_dot, R, Q_R, Q_L, Q_A, theta_A, W, W_R, W_R_dot, W_L, ...
        W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, f_tau, tau, Euler_R, ...
        Euler_R_dot] = simulate_control(t, X0, X_ref, WK_R, WK_L, INSECT, ...
        N, N_periods, N_single, N_iters, N_per_iter, problem, Weights, ...
        get_args, param_type, p, m);
    time_taken = toc;

    plot(0:1:N_periods*N_iters, cost);

    case 'monte_carlo'
    N_sims = 4; % 1 takes 4200s
    problem.options.Display = 'none';
    tic;
    [X0_pert, dX, cost, opt_param] = monte_carlo(N_sims, t, X_ref, WK_R, WK_L, INSECT, ...
            N, N_periods, N_single, N_iters, N_per_iter, problem, Weights, ...
            get_args, param_type, p, m);
    time_taken = toc;
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

function [X0_pert, dX, cost, opt_param] = monte_carlo(N_sims, t, X_ref, WK_R, WK_L, INSECT, ...
        N, N_periods, N_single, N_iters, N_per_iter, problem, Weights, ...
        get_args, param_type, p, m)
%%
e1 = [1 0 0]'; e2 = [0 1 0]'; e3 = [0 0 1]';
dX = zeros(N_sims, 12);
cost = zeros(N_sims, N_periods, 1+N_iters);
opt_param = zeros(N_sims, N_periods, length(problem.x0)*m/p);
X0_pert = zeros(N_sims, N_periods, 18);
idx_cost = zeros(N_periods, 1+N_iters);
for period=1:3
    idx_cost(period,1:(1+N_iters)) = (1+(period-1)*N_iters):(1+period*N_iters);
end
idx_X0 = 1:N_single:(N_single*N_periods);
des_X0 = X_ref(1, :);
des_R0 = reshape(des_X0(4:12), 3, 3);

par_pool = gcp;
nWorkers = par_pool.NumWorkers;
ticBytes(par_pool);
parfor n=1:nWorkers
    rng(n);
end
pause(1);

for i = 1:N_sims
    dX(i,:) = Weights.PerturbVariables .* randn(1, 12);
end
parfor i = 1:N_sims
    % perturbation
    dXi = dX(i,:);
    X0 = [des_X0(1:3)+dXi(1:3),...
        reshape(des_R0*expmhat(dXi(6)*e3)*expmhat(dXi(5)*e2)*expmhat(dXi(4)*e1),1,9),...
        des_X0(13:18) + dXi(7:12)]';
    [cost_arr, opt_param_arr, ~, X] = simulate_control(t, X0, X_ref, WK_R, WK_L, INSECT, ...
        N, N_periods, N_single, N_iters, N_per_iter, problem, Weights, ...
        get_args, param_type, p, m);
    cost(i,:,:) = cost_arr(idx_cost);
    opt_param(i,:,:) = opt_param_arr';
    X0_pert(i,:,:) = X(idx_X0, :);
    fprintf('Completed simulation %d\n', i);
end

cost = reshape(cost, N_sims*N_periods, 1+N_iters);
opt_param = reshape(opt_param, N_sims*N_periods, length(problem.x0)*m/p);
X0_pert = reshape(X0_pert, N_sims*N_periods, 18);

tocBytes(par_pool);

end

function [cost, opt_param, dang, X, x, x_dot, R, Q_R, Q_L, Q_A, theta_A, W, W_R, W_R_dot, W_L, ...
    W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, f_tau, tau, Euler_R, ...
    Euler_R_dot] = simulate_control(t, X0, X_ref, WK_R, WK_L, INSECT, ...
    N, N_periods, N_single, N_iters, N_per_iter, problem, Weights, ...
    get_args, param_type, p, m)
%%
cost = zeros(1, 1+N_periods*N_iters);
cost(1) = sqrt(sum((Weights.OutputVariables .* (X0' - X_ref(1, :))).^2));
X = zeros(1+N_single*N_periods, 18);
N_p = length(problem.x0) / p;
dang = zeros(N_p, 1+N_single*N_periods);
opt_param = zeros(N_p*m, N_periods*N_iters/m);

for period=1:N_periods
    for iter=1:m:N_iters
%         cost0 = cost(end);
%         problem.nonlcon = @(param) cost_cons(...
%             param, t(idx), X0, INSECT, WK_R, WK_L, X_ref(idx,:), Weights, param_type, cost0, varargin{:});
        idx_opt = (1+(period-1)*N_single+(iter-1)*N_per_iter):(1+(period-1)*N_single+(iter+p-1)*N_per_iter);
        dang0 = dang(:,idx_opt(1));
        varargin = get_args(t(idx_opt(1)), dang(:,idx_opt(1)));
        problem.objective = @(param) cost_fun(param, t(idx_opt), X0, ...
            INSECT, WK_R, WK_L, X_ref(idx_opt,:), Weights, param_type, ...
            get_args, p, N_p, N_per_iter, varargin{:});
        % Equality constraints work only for special values of p, m
        I_66 = eye(6); O_66 = zeros(6);
        problem.Aeq = [repmat(I_66, 1, N_iters), repmat(O_66, 1, N_iters);
                       repmat(O_66, 1, N_iters), repmat(I_66, 1, N_iters)];
        problem.beq = zeros(12, 1);
%         problem.nonlcon = @(param) param_cons(param, p, N_iters, dang0, idx_opt);
        [param, fval, exitflag] = fmincon(problem);
        fprintf('Exit flag is %d\n', exitflag);
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
            cost(1+(N_periods-1)*N_iters+(iter-1)+con) = sqrt(...
                sum((Weights.OutputVariables .* (X(idx(end), :) - X_ref(idx(end), :))).^2));
        end
        if p > 4
            problem.options.MaxFunctionEvaluations = 3000;
        end
    end
    fprintf('Period number %d completed\n', period);
end

x=X(:,1:3)';
x_dot=X(:,13:15)';

end

function J = cost_fun(param, t, X0, INSECT, WK_R, WK_L, X_ref, Weights, param_type, ...
    get_args, p, N_p, N_per_iter, varargin)
%%
J = 0;
for pred=1:p
    param_idx = (1+(pred-1)*N_p):(pred*N_p);
    param_p = param(param_idx);
    idx = (1+(pred-1)*N_per_iter):(1+pred*N_per_iter);

%     if length(idx) == 2
%         [~, X_temp] = ode45(@(t,X) eom_param(INSECT, WK_R, WK_L, t, X, param_p, param_type, varargin{:}), ...
%         t(idx), X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
%         X(idx,:) = [X_temp(1,:); X_temp(end,:)];
    [~, X_temp] = ode45(@(t,X) eom_param(INSECT, WK_R, WK_L, t, X, param_p, param_type, varargin{:}), ...
    t(idx), X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
    if size(X_temp, 1) ~= (N_per_iter+1)
        J = 1e10;
        return;
    else
        X(idx,:) = X_temp;
    end
    X0 = X(idx(end),:)';
    dang0 = param_type(param_p, t(idx(end)), varargin{:});
    varargin = get_args(t(idx(end)), dang0);
    J = J + Weights.PredictionHorizon(pred) * ...
        sqrt(sum((Weights.OutputVariables .* (X(idx(end),:) - X_ref(idx(end),:))).^2));
end
% J = sqrt(sum((Weights.OutputVariables .* (X(end,:) - X_ref(end,:))).^2));

% [~, X] = ode45(@(t,X) eom_param(INSECT, WK_R, WK_L, t, X, param, param_type, varargin{:}), ...
%     t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6)); % [t(1), t(end)] instead of [t]
% J = sum(sum((Weights.OutputVariables .* (X-X_ref)).^2));

end

function [c, ceq] = param_cons(param, p, N_iters, dang0, idx_opt)
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
