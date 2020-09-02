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

% Values obtained from parametric study
load('parametric_study_xR.mat', 'params');
des.params = params;
% Observing parametric study plots
% des.params.mat_aero(1,[3,5,6],:) = 0;
% des.params.mat_aero(2,[1,2,4,6],:) = 0;
% des.params.mat_aero(3,[3,5,6],:) = 0;
% des.params.mat_aero(3,[1,2,4],1) = 1; % A big number?
% des.params.mat_aero(4,[1,2,4,5,6],:) = 0;
% des.params.mat_aero(5,[3,5,6],:) = 0;
% des.params.mat_aero(6,[1,2,4],:) = 0;

% [f_a, M_a] vs [dphi_ms, dtheta_0s, dphi_mk, dphi_0s, dtheta_0k, dpsi_0k]
if rank(des.params.mat_aero(:,:,1), 1e-10) < 6
    error('Parameters are not able to generate aerodynamic forces and moments');
end

load('sim_QS_xR_hover.mat', 'solutions');
WK_arr = solutions(1).X; % Controllable : (2), Uncontrollable : (1)
[WK, des.x_dot0, des.R0, des.W0] = get_WK0(WK, WK_arr);
des.X0=[des.x0; reshape(des.R0,9,1); des.x_dot0; des.W0];

% des = rmfield(des, {'WK', 'INSECT', 't', 'x', 'x_dot', 'R', 'W', 'f_tau'}); % 'x_fit', 'x_dot_fit'

%% Gains
pol = poly([-7.8 + 19i, -7.8 - 19i, -0.003]);
if ~all(real(roots(pol)) < 0)
    error('The chosen gains are not suitable');
end
% Optimized gs = [427.1529   15.6076  13.4983];
% gains.Kp_pos = pol(3); gains.Kd_pos = pol(2); gains.Ki_pos = pol(4);
gains.Kp_pos = 4; gains.Kd_pos = 2; gains.Ki_pos = 1;
gains.KR = 4; gains.KOm = 2; gains.KI = 1; gains.cI = 1e-1;

%% Single simulation
rng default;
eps = 1e-3;
% dx0 = rand(3,1)*eps;
dx0 = zeros(3, 1);
dx_dot0 = zeros(3, 1);
x0 = des.x0 + dx0;
x_dot0 = des.x_dot0 + dx_dot0;
int_d_x0 = zeros(3, 1);
R0 = des.R0*expmhat(rand(1)*1e-1*[1,0,0]);
W0 = des.W0 + [10; 100; 0]*eps;
int_att0 = zeros(3,1);
X0 = [x0; reshape(R0, 9, 1); x_dot0; W0;] + [1e-2*ones(3,1); zeros(9,1); 0.05*ones(3,1); 1*ones(3,1)].*rand(18,1);
% wt = 0;

%% Optimization
N_periods = 1; % 1, 3, 15
N_iters = 2; % 2, 4 per period
max_f_iter = 500; % 100, 500 per optimization
N_single = 100; % per time period
N_per_iter = N_single / N_iters;
N = N_single*N_periods + 1;
err_bound = 1e-4; % Convergence criterion is a f(N, N_periods, WK)
T = N_periods/WK.f;
t = linspace(0,T,N);
bound_param = 0.1; % Parameter bound; Use 0.25?

X_ref = zeros(N, 18);
X_ref0 = des.X0';
for period=1:N_periods
    idx = (1+(period-1)*N_single):(1+period*N_single);
    if period == 1
        [~, X_ref(idx, :)]=ode45(@(t,X) eom_QS_xR(INSECT, WK, WK, t,X), t(idx), X_ref0, odeset('AbsTol',1e-6,'RelTol',1e-6));
        idx_prev = idx;
    else
        X_ref(idx, :) = X_ref(idx_prev, :);
        idx_prev = idx;
    end
end
des.x_fit_t = X_ref(:, 1:3)'; des.x_dot_fit_t = X_ref(:, 13:15)';
des.R_fit_t = reshape(X_ref(:, 4:12)', 3, 3, N); des.W_fit_t = X_ref(:, 16:18)';

Weights.PerturbVariables = [0.025*ones(1,3), 0.25*ones(1,3), 0.1*ones(1,3), 1*ones(1,3)]; % To multiply the perturbations
Weights.OutputVariables = [25*ones(1,3), 2.5*ones(1,9), 5*ones(1,3), 1*ones(1,3)]; % 10, 2.5, 2.5, 2
% Weights.Time = logspace(-1, 1, N_single+1)';

WK_R = WK; WK_L = WK;

% FMINCON
options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter',...
    'MaxFunctionEvaluations',max_f_iter,'UseParallel',true); %'PlotFcn',@optimplotfval,
N_dang = 6;
dang0 = zeros(N_dang, 1);

param_type = @constant_func;
get_args = @constant_args;
lb = -bound_param * ones(1, N_dang);
ub = bound_param * ones(1, N_dang);

% param_type = @linear_func;
% get_args = @linear_args;
% lb = -2*bound_param*WK.f*N_iters * ones(1, N_dang);
% ub = 2*bound_param*WK.f*N_iters * ones(1, N_dang);

problem.x0 = dang0;
problem.lb = lb;
problem.ub = ub;
problem.solver = 'fmincon';
problem.options = options;

%% Single simulation
load('sim_QS_xR_hover_control_opt_mcdata', 'X0_pert', 'cost');
X0_unstable = X0_pert(cost(:, 3) > cost(:, 1), :);
X0 = X0_unstable(end, :)';
problem.options.MaxFunctionEvaluations = 1000;
problem.options.Algorithm = 'active-set'; % interior-point, sqp, fastest : active-set
N_iters = 1; % 2, 4 per period
N_per_iter = N_single / N_iters;
tic;
[cost, dang, X, x, x_dot, R, Q_R, Q_L, Q_A, theta_A, W, W_R, W_R_dot, W_L, ...
    W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, f_tau, tau, Euler_R, ...
    Euler_R_dot] = simulate_control(t, X0, X_ref, WK_R, WK_L, INSECT, ...
    N, N_periods, N_single, N_iters, N_per_iter, problem, Weights, ...
    get_args, param_type);
toc;

plot(0:1:N_periods*N_iters, cost);

%% Monte carlo
% N_sims = 512; % 512 takes 36478s, 16 ; 71.2s per sim
% problem.options.Display = 'none';
% [X0_pert, dX, cost, dang] = monte_carlo(N_sims, t, X_ref, WK_R, WK_L, INSECT, ...
%         N, N_periods, N_single, N_iters, N_per_iter, problem, Weights, ...
%         get_args, param_type);

%%
% Get a list of all variables
allvars = whos;
% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics object
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
% Pass these variable names to save
save(filename, allvars(tosave).name)
evalin('base',['load ' filename]);

end

function [X0_pert, dX, cost, dang] = monte_carlo(N_sims, t, X_ref, WK_R, WK_L, INSECT, ...
        N, N_periods, N_single, N_iters, N_per_iter, problem, Weights, ...
        get_args, param_type)
%%
e1 = [1 0 0]'; e2 = [0 1 0]'; e3 = [0 0 1]';
X0_pert = zeros(N_sims, 18);
dX = zeros(N_sims, 12);
cost = zeros(N_sims, N_periods*N_iters+1);
dang = zeros(N_sims, length(problem.x0), N_periods*N_iters);
idx_dang = [1, 1+N_per_iter];
des_X0 = X_ref(1, :);
des_R0 = reshape(des_X0(4:12), 3, 3);

tic;
par_pool = gcp;
nWorkers = par_pool.NumWorkers;
% ticBytes(par_pool);
parfor n=1:nWorkers
    rng(n);
end
pause(1);

parfor i = 1:N_sims
    % perturbation
    dXi = Weights.PerturbVariables .* randn(1, 12);
    X0 = [des_X0(1:3)+dXi(1:3),...
        reshape(des_R0*expmhat(dXi(6)*e3)*expmhat(dXi(5)*e2)*expmhat(dXi(4)*e1),1,9),...
        des_X0(13:18) + dXi(7:12)]';
    dX(i,:) = dXi; X0_pert(i,:) = X0';
    [cost(i,:), dang_arr] = simulate_control(t, X0, X_ref, WK_R, WK_L, INSECT, ...
        N, N_periods, N_single, N_iters, N_per_iter, problem, Weights, ...
        get_args, param_type);
    dang(i,:,:) = dang_arr(:, idx_dang);
    fprintf('Completed simulation %d\n', i);
end

% tocBytes(par_pool);
toc;

end

function [cost, dang, X, x, x_dot, R, Q_R, Q_L, Q_A, theta_A, W, W_R, W_R_dot, W_L, ...
    W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, f_tau, tau, Euler_R, ...
    Euler_R_dot] = simulate_control(t, X0, X_ref, WK_R, WK_L, INSECT, ...
    N, N_periods, N_single, N_iters, N_per_iter, problem, Weights, ...
    get_args, param_type)
%%
cost = [sum((Weights.OutputVariables .* (X0' - X_ref(1,:))).^2)];
X = zeros(N, 18);
dang = zeros(length(problem.x0), N);

for period=1:N_periods
%     idx = (1+(period-1)*N_single):(1+period*N_single);
    for iter=1:N_iters
        idx = (1+(period-1)*N_single+(iter-1)*N_per_iter):(1+(period-1)*N_single+iter*N_per_iter);
        varargin = get_args(t, dang, idx);
        problem.objective = @(param) cost_fun(param, t(idx), X0, ...
            INSECT, WK_R, WK_L, X_ref(idx,:), Weights, param_type, varargin{:});
%         cost0 = cost(end);
%         problem.nonlcon = @(param) cost_cons(...
%             param, t(idx), X0, INSECT, WK_R, WK_L, X_ref(idx,:), Weights, param_type, cost0, varargin{:});
        [param, fval, exitflag] = fmincon(problem);
% %           MULTISTART
%         problem.options.Display = 'none';
%         ptmatrix(1, :) = problem.x0';
%         N_points = 10;
%         ptmatrix(2:N_points, :) = problem.lb + rand(N_points-1, length(problem.x0)) .* (problem.ub - problem.lb);
%         tpoints = CustomStartPointSet(ptmatrix);
%         ms = MultiStart('Display','iter','PlotFcn',@gsplotbestf,'MaxTime',0.5*3600);
%         [param, fval, exitflag, output, solutions] = run(ms, problem, tpoints);
        
        fprintf('Exit flag is %d\n', exitflag);
        [~, X(idx,:)] = ode45(@(t,X) eom_param(INSECT, WK_R, WK_L, t, X, param, param_type, varargin{:}), ...
            t(idx), X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
        cost = [cost, fval];
        for k=idx
            dang(:,k) = param_type(param, t(k), varargin{:});
            [WK_R_new, WK_L_new] = get_WK(WK_R, WK_L, dang);
            [X_dot(:,k), R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k),...
                theta_A(k), W(:,k), W_R(:,k), W_R_dot(:,k), W_L(:,k),...
                W_L_dot(:,k), W_A(:,k), W_A_dot(:,k), F_R(:,k), F_L(:,k), M_R(:,k),...
                M_L(:,k), f_a(:,k), f_g(:,k), f_tau(:,k), tau(:,k), ...
                Euler_R(:,k), Euler_R_dot(:,k)] = ...
                eom_QS_xR(INSECT, WK_R_new, WK_L_new, t(k), X(k, :)');
        end
        X0 = X(idx(end),:)';
    end
%     dang0 = dang(:, idx(end));
    fprintf('Period number %d completed\n', period);
end

x=X(:,1:3)';
x_dot=X(:,13:15)';

end

function J = cost_fun(param, t, X0, INSECT, WK_R, WK_L, X_ref, Weights, param_type, varargin)
[~, X] = ode45(@(t,X) eom_param(INSECT, WK_R, WK_L, t, X, param, param_type, varargin{:}), ...
    [t(1), t(end)], X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
% [~, X] = ode45(@(t,X) eom_QS_xR(INSECT, WK_R, WK_L, t, X), ...
%     t, X0(1:18), odeset('AbsTol',1e-6,'RelTol',1e-6));

% J = sum(((X(end,:) - X_ref(end,:))).^2);
J = sum((Weights.OutputVariables .* (X(end,:) - X_ref(end,:))).^2);
% J = sum(sum((Weights.OutputVariables .* (X-X_ref) .* Weights.Time).^2));
% J = sum(sum(((X-X_ref)).^2));

end

function [c, ceq] = cost_cons(param, t, X0, INSECT, WK_R, WK_L, X_ref, Weights, param_type, cost0, varargin)
[~, X] = ode45(@(t,X) eom_param(INSECT, WK_R, WK_L, t, X, param, param_type, varargin{:}), ...
    [t(1), t(end)], X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
J = sum((Weights.OutputVariables .* (X(end,:) - X_ref(end,:))).^2);
c(1) = J - cost0;
ceq = [];
end

function dXdt = eom_param(INSECT, WK_R, WK_L, t, X, param, param_type, varargin)
dang = param_type(param, t, varargin{:});
[WK_R, WK_L] = get_WK(WK_R, WK_L, dang);
dXdt = eom_QS_xR(INSECT, WK_R, WK_L, t, X);
end

function dang = constant_func(param, t, varargin)
dang = param;
end

function args = constant_args(t, dang, idx)
args = {};
end

function dang = linear_func(param, t, varargin)
t0 = varargin{1};
dang0 = varargin{2};
dang = dang0 + param*(t - t0);
end

function args = linear_args(t, dang, idx)
args{1} = t(idx(1));
args{2} = dang(:, idx(1));
end

function [WKR_new, WKL_new] = get_WK(WKR_old, WKL_old, dang)
    WKR_new = WKR_old; WKL_new = WKL_old;
    WKR_new.phi_m = WKR_new.phi_m + dang(1) + dang(3);
    WKL_new.phi_m = WKL_new.phi_m + dang(1) - dang(3);
    WKR_new.theta_0 = WKR_new.theta_0 + dang(2) + dang(5);
    WKL_new.theta_0 = WKL_new.theta_0 + dang(2) - dang(5);
    WKR_new.phi_0 = WKR_new.phi_0 + dang(4);
    WKL_new.phi_0 = WKL_new.phi_0 + dang(4);
    WKR_new.psi_m = WKR_new.psi_m + dang(6);
    WKL_new.psi_m = WKL_new.psi_m - dang(6);
%     WKR_new.psi_m = WKR_new.psi_m + dang(4) + dang(6);
%     WKL_new.psi_m = WKL_new.psi_m + dang(4) - dang(6);
    
%     WKR_new.theta_A_m = WKR_new.theta_A_m + dang(7);
%     WKL_new.theta_A_m = WKL_new.theta_A_m + dang(7);
end

function MPC
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
