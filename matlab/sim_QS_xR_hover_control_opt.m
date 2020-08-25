function sim_QS_xR_hover_control_opt
% simulate the thorax (x) trajectory along with a controller for
% given thorax attiude, wing kinematics, abdomen attitude.

evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
des = load('sim_QS_xR_hover.mat', 'INSECT', 't', 'x', 'x_dot','R','W',...
    'x0', 'x_dot0', 'R0', 'W0', 'WK');

filename = 'sim_QS_xR_hover_control_opt';
INSECT = des.INSECT;
WK = des.WK;
des.x_fit = cell(3, 1); des.x_dot_fit = cell(3, 1);
des.R_fit = cell(3, 3); des.W_fit = cell(3, 1);
% 'fourier8', 'cubicinterp'
for i=1:3
    des.x_fit{i} = fit(des.t, des.x(i, :)', 'fourier8');
    des.x_dot_fit{i} = fit(des.t, des.x_dot(i, :)', 'fourier8');
    des.W_fit{i} = fit(des.t, des.W(i, :)', 'fourier8');
    for j=1:3
        des.R_fit{i,j} = fit(des.t, squeeze(des.R(i, j, :)), 'fourier8');
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

N_period = 10; %15
N_single = 100;
N = N_single*N_period + 1;
err_bound = 1e-4; % Convergence criterion is a f(N, N_period, WK)
T = N_period/WK.f;
t = linspace(0,T,N);
bound_param = 0.1; % Parameter bound; Use 0.25?

des.x_fit_t = zeros(3, N); des.x_dot_fit_t = zeros(3, N);
des.R_fit_t = zeros(3, 3, N); des.W_fit_t = zeros(3, N);
for i=1:3
    des.x_fit_t(i, :) = des.x_fit{i}(t);
    des.x_dot_fit_t(i, :) = des.x_dot_fit{i}(t);
    des.W_fit_t(i, :) = des.W_fit{i}(t);
    for j=1:3
        des.R_fit_t(i, j, :) = des.R_fit{i,j}(t);
    end
end
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
X0 = [x0; reshape(R0, 9, 1); x_dot0; W0;];% + eps*rand(18,1);
% wt = 0;

%% Optimization
X_ref = [des.x_fit_t' squeeze(reshape(des.R_fit_t, 1, 9, N))' des.x_dot_fit_t' des.W_fit_t'];

Weights.OutputVariables = [0.25*ones(1,3), 0.25*ones(1,9), 1*ones(1,3), 2*ones(1,3)];
Weights.Time = logspace(-1, 1, N_single+1)';
% nlobj.Weights.ManipulatedVariables = 0.05 * ones(1, nu);
% nlobj.Weights.ManipulatedVariablesRate = 1e10 * ones(1, nu);

A = []; b = []; Aeq = []; beq = [];
N_dang = 6;
dang0 = zeros(N_dang,1);
dang = zeros(N_dang, N);
X = zeros(N, 18);
WK_R = WK; WK_L = WK;
lb = -bound_param * ones(1, N_dang);
ub = bound_param * ones(1, N_dang);
% nonlcon = @(WK_arr) traj_condition(WK_arr, WK, INSECT, N, x0, final_pos);

tic;
rng default; % For reproducibility

% FMINCON
options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter',...
    'MaxFunctionEvaluations',100,'UseParallel',true); %'PlotFcn',@optimplotfval,

for period=1:N_period
    idx = (1+(period-1)*N_single):(1+period*N_single);
    [dang_val, fval, exitflag, output] = fmincon(@(dang) cost_fun(dang, t(idx), X0, INSECT, WK_R, WK_L, X_ref(idx,:), Weights),...
    dang0,A,b,Aeq,beq,lb,ub,[],options);
    [WK_R_new, WK_L_new] = get_WK(WK_R, WK_L, dang_val);
    dang(:,idx) = repmat(dang_val, 1, N_single+1);
    [~, X(idx,:)] = ode45(@(t,X) eom_QS_xR(INSECT, WK_R_new, WK_L_new, t, X), ...
        t(idx), X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
    for k=idx
        [X_dot(:,k), R(:,:,k), Q_R(:,:,k), Q_L(:,:,k), Q_A(:,:,k),...
            theta_A(k), W(:,k), W_R(:,k), W_R_dot(:,k), W_L(:,k),...
            W_L_dot(:,k), W_A(:,k), W_A_dot(:,k), F_R(:,k), F_L(:,k), M_R(:,k),...
            M_L(:,k), f_a(:,k), f_g(:,k), f_tau(:,k), tau(:,k), ...
            Euler_R(:,k), Euler_R_dot(:,k)] = ...
            eom_QS_xR(INSECT, WK_R_new, WK_L_new, t(k), X(k, :)');
    end
    X0 = X(idx(end),:)';
%     dang0 = dang(:, idx(end));
    fprintf('Period number %d completed\n', period);
end

fprintf('Optimization has been completed\n');
disp(output);
toc;

x=X(:,1:3)';
x_dot=X(:,13:15)';

% tic;
% [err_pos, N_conv, x, x_dot, int_d_x, R, Q_R, Q_L, Q_A, theta_A, W, W_R, W_R_dot, W_L, ...
%     W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, f_tau, tau, Euler_R, ...
%     Euler_R_dot, err_xR, dang] =  simulate_control(gains, WK, INSECT, des, X0, N, ...
%     N_single, N_period, t, wt, bound_param, err_bound);
% toc;

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

function J = cost_fun(dang, t, X0, INSECT, WK_R, WK_L, X_ref, Weights)
[WK_R, WK_L] = get_WK(WK_R, WK_L, dang);

[~, X] = ode45(@(t,X) eom_QS_xR(INSECT, WK_R, WK_L, t, X), ...
    t, X0(1:18), odeset('AbsTol',1e-6,'RelTol',1e-6));

% J = sum(((X(end,:) - X_ref(end,:))).^2);
J = sum((Weights.OutputVariables .* (X(end,:) - X_ref(end,:))).^2);
% J = sum(sum((Weights.OutputVariables .* (X-X_ref) .* Weights.Time).^2));
% J = sum(sum(((X-X_ref)).^2));


end

function [WKR_new, WKL_new] = get_WK(WKR_old, WKL_old, dang)
    WKR_new = WKR_old; WKL_new = WKL_old;
    WKR_new.phi_m = WKR_new.phi_m + dang(1) + dang(3);
    WKL_new.phi_m = WKL_new.phi_m + dang(1) - dang(3);
    WKR_new.theta_0 = WKR_new.theta_0 + dang(2) + dang(5);
    WKL_new.theta_0 = WKL_new.theta_0 + dang(2) - dang(5);
%     WKR_new.phi_0 = WKR_new.phi_0 + dang(4);
%     WKL_new.phi_0 = WKL_new.phi_0 + dang(4);
%     WKR_new.psi_m = WKR_new.psi_m + dang(6);
%     WKL_new.psi_m = WKL_new.psi_m - dang(6);
    WKR_new.psi_m = WKR_new.psi_m + dang(4) + dang(6);
    WKL_new.psi_m = WKL_new.psi_m + dang(4) - dang(6);
    
%     WKR_new.theta_A_m = WKR_new.theta_A_m + dang(7);
%     WKL_new.theta_A_m = WKL_new.theta_A_m + dang(7);
end
