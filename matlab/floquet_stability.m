function floquet_stability
% Study the stability of a periodic trajecoty using floquet theory.

evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
filename='floquet_stability';

%% Linearized dynamics
load('sim_QS_xR_hover.mat', 'INSECT', 'WK', 'X0');
bool_sort_mus=false; % false if n ~= 6

N = 1001;
N_period = 2;
T=N_period/WK.f;
ix_d = (N-1)/N_period;
t=linspace(0,T,N);
dt = t(2) - t(1);
epsilon = 1e-8;

% n is #perturbation states; n_vars is #actual states
% n = 3; n_vars = 6; % for nominal hover with @eom_hover_vel if position is not periodic
% [delta_mat, F_linear] = sim_pert(@eom_hover_vel, n, n_vars, INSECT, WK, X0, N, t, epsilon);
% n = 6; n_vars = 6; % for nominal hover with @eom_hover
% [delta_mat, F_linear] = sim_pert(@eom_hover, n, n_vars, INSECT, WK, X0, N, t, epsilon);

% load('sim_QS_x_hover_control.mat', 'des', 'gains', 'wt', 'bound_param');
% int_d_x0 = zeros(3,1);
% X0 = [X0; int_d_x0;];
% n = 9; n_vars = 9; % for controlled hover with @eom_hover_control
% [delta_mat, F_linear] = sim_pert(@eom_hover_control, n, n_vars, INSECT, WK, X0, N, t, epsilon, des, gains, wt, bound_param);

% n = 12; n_vars = 18; % for (x,R) stability with @eom_hover_xR; n=8 for just pitch, n=12 for all modes;
% [delta_mat, F_linear] = sim_pert(@eom_hover_xR, n, n_vars, INSECT, WK, X0, N, t, epsilon);
load('sim_QS_xR_hover_control.mat', 'des', 'gains', 'wt', 'bound_param');
% gains.KR = 40; gains.KOm = 10;

load('sim_QS_xR_hover.mat', 'solutions');
WK_arr = solutions(2).X; % Controllable : (2), Uncontrollable : (1)
[WK, des.x_dot0, des.R0, des.W0] = get_WK0(WK, WK_arr);
des.X0=[des.x0; reshape(des.R0,9,1); des.x_dot0; des.W0];

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
X0 = des.X0;

X0 = [X0; zeros(3,1); zeros(3,1)];
n = 18; n_vars = 24; % for (x,R) control with @eom_hover_xR_control
tic;
[delta_mat, F_linear] = sim_pert(@eom_hover_xR_control, n, n_vars, INSECT, WK, X0, N, t, epsilon, des, gains, wt, bound_param);
toc;
% [delta_mat, F_linear, gains] = optimized_gains(@eom_hover_xR_control, n, n_vars, INSECT, WK, X0, N, t, epsilon, des, gains, wt, bound_param);

B = zeros(n, n, 1+ix_d);
start_ix = max(1, round((N_period-2)/N_period * N));
for i=start_ix:(start_ix+ix_d)
    j = i-start_ix+1;
    d_F = (F_linear(:, :, i) - F_linear(:, :, i+ix_d))./ F_linear(:, :, i);
    %
    max_F = min(max(F_linear(:,:,i), [], 'all'), max(F_linear(:,:,i+ix_d), [], 'all'));
    d_F(max(abs(F_linear(:, :, i+ix_d)/max_F), abs(F_linear(:, :, i)/max_F)) < 1e-4) = 0; % Insignificant values
    idx_small = max(abs(F_linear(:, :, i+ix_d)/max_F), abs(F_linear(:, :, i)/max_F)) < 1e-1;
    d_F(d_F(idx_small) < 0.25) = 0; % Small values
    %
    if(~all(abs(d_F) < 5e-2, 'all'))
        disp(i)
        disp(d_F)
    end
    B(:, :, j) = delta_mat(:, :, i) \ delta_mat(:, :, i+ix_d);
end
B = B(:, :, round(end/2));
[e_vecs, rhos] = eig(B);
mus = log(abs(diag(rhos))) * WK.f;

char_soln_mat = zeros(n, n, N);
per_val_mat = zeros(n, n, N);

for i=1:N
    char_soln_mat(:, :, i) = delta_mat(:, :, i) * e_vecs;
    Y_0 = diag(exp(mus*t(i)));
    per_val_mat(:, :, i) = char_soln_mat(:, :, i) / Y_0;
end

%%
if bool_sort_mus
    [mus_s, idx_mus] = sort(mus,'descend');
    e_vecs_s = e_vecs(:,idx_mus);
    mus_sort = mus_s;
    e_vecs_sort = e_vecs_s;
    for i=4:6
        if all(abs(e_vecs_s([1,3,4,6],i) - 0) < 1e-10)
            lat_idx = i;
        end
    end
    mus_sort(end) = mus_s(lat_idx); mus_sort(lat_idx) = mus_s(end);
    e_vecs_sort(:,end) = e_vecs_s(:,lat_idx); e_vecs_sort(:,lat_idx) = e_vecs_s(:,end);
    [~, idx_sort] = ismember(mus_sort, mus); % Sort according to modes
end

% Get a list of all variables
allvars = whos;
% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics object
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
% Pass these variable names to save
save(filename, allvars(tosave).name)
evalin('base',['load ' filename]);
end

function [delta_mat, F_linear, gains] = optimized_gains(eom, n, n_vars, INSECT, WK, X0, N, t, epsilon, varargin)
%% The optimization algorithm
A = []; b = []; Aeq = []; beq = [];
% Initial value of gains_arr = [Kp_pos, Kd_pos, Ki_pos, KR, KOm, KI, cI]
gains0 = varargin{2};
gains_arr0 = struct2array(gains0);
lb = [0, 0, 0, 0, 0, 0, 0];
ub = [1000, 500, 100, 100, 100, 50, 10];
nonlcon = @(gains_arr) gains_condition(gains_arr);

tic;
rng default; % For reproducibility

% FMINCON
options = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter',...
    'MaxFunctionEvaluations',30,'PlotFcn',@optimplotfval,'UseParallel',true);
[gains_arr, fval, exitflag, output] = fmincon(@(gains_arr) objective_func(gains_arr, eom, n, n_vars, INSECT, WK, X0, N, t, epsilon, varargin{:}),...
    gains_arr0,A,b,Aeq,beq,lb,ub,nonlcon,options);

% % MULTISTART, PARTICLESWARM
% ptmatrix(1, :) = gains_arr0;
% N_points = 5;
% ptmatrix(2:N_points, :) = lb + rand(N_points-1, length(gains_arr0)) .* (ub - lb);
% tpoints = CustomStartPointSet(ptmatrix);
% ms = MultiStart('Display','iter','PlotFcn',@gsplotbestf,'MaxTime',1*3600);
% options = optimoptions(@fmincon,'Algorithm','sqp',...
%     'UseParallel',false,'MaxIterations',2000,'MaxFunctionEvaluations',6000);%'ConstraintTolerance',1e-5,'StepTolerance',1e-8,'OptimalityTolerance',1e-5);
% problem = createOptimProblem('fmincon','objective',@(gains_arr) objective_func(gains_arr, eom, n, n_vars, INSECT, WK, X0, N, t, epsilon, varargin{:}),...
%     'x0',gains_arr0,'lb',lb,'ub',ub,'nonlcon',nonlcon,'options',options);
% [gains_arr, fval, exitflag, output, solutions] = run(ms, problem, tpoints);

fprintf('Optimization has been completed\n');
disp(output);
toc;

%%
gains = struct('Kp_pos', gains_arr(1), 'Kd_pos', gains_arr(2), 'Ki_pos', gains_arr(3),...
        'KR', gains_arr(4), 'KOm', gains_arr(5), 'KI', gains_arr(6) , 'cI', gains_arr(7));
varargin{2} = gains;
[delta_mat, F_linear] = sim_pert(eom, n, n_vars, INSECT, WK, X0, N, t, epsilon, varargin{:});

end

function [c, ceq] = gains_condition(gains_arr)
    c = real(roots([1, gains_arr(2), gains_arr(1), gains_arr(3)]));
    ceq = [];
end

function [J] = objective_func(gains_arr, eom, n, n_vars, INSECT, WK, X0, N, t, epsilon, varargin)
%% Objective function
gains = struct('Kp_pos', gains_arr(1), 'Kd_pos', gains_arr(2), 'Ki_pos', gains_arr(3),...
        'KR', gains_arr(4), 'KOm', gains_arr(5), 'KI', gains_arr(6) , 'cI', gains_arr(7));
varargin{2} = gains;
delta_mat = sim_pert(eom, n, n_vars, INSECT, WK, X0, N, t, epsilon, varargin{:});
B = delta_mat(:, :, round(N/4)) \ delta_mat(:, :, round(N/4) + 2*round(N/4));
[~, rhos] = eig(B);
mus = log(abs(diag(rhos))) * WK.f;
J = sum(exp(mus(mus>0)) - 1);

if isnan(J)
    J = 1/eps;
end

end

function [delta_mat, F_linear] = sim_pert(eom, n, n_vars, INSECT, WK, X0, N, t, epsilon, varargin)
%%
delta0 = eye(n)*epsilon;
X0 = [X0(1:n_vars); reshape(delta0, n^2, 1);];

[t,X] = ode45(@(t,X) eom(n, n_vars, INSECT, WK, WK, t, X, varargin), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
F_linear = zeros(n, n, N);
parfor k=1:N
    [~, F_linear(:, :, k)] = eom(n, n_vars, INSECT, WK, WK, t(k), X(k, :)', varargin);
end

delta_mat = reshape(X(:,n_vars+1:(n^2+n_vars))', n, n, N);

end

function [X_dot, F_linear]= eom_hover_vel(n, n_vars, INSECT, WK_R, WK_L, t, X, varargin)
%% Stability analysis for only velocity perturbations
x=X(1:3);
x_dot=X(4:6);
delta_mat=reshape(X(n_vars+1:(n^2+n_vars)), n, n);

X = X(1:6);
[X_dot, R, Q_R, Q_L, Q_A, theta_B, theta_A, W, W_dot, W_R, ...
    W_R_dot, W_L, W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, ...
    f_tau, tau] = eom_QS_x(INSECT, WK_R, WK_L, t, X);
JJ = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
JJ_11 = inertia_sub_decompose_3_12(JJ);

[d_L_R, d_L_L, d_D_R, d_D_L]=wing_QS_aerodynamics_linearized(INSECT, W_R, W_L, W_R_dot, W_L_dot, x_dot, R, W, Q_R, Q_L);
d_F_R = d_L_R + d_D_R;
d_F_L = d_L_L + d_D_L;
F_linear = JJ_11\ R * (Q_R*d_F_R + Q_L*d_F_L);
X_dot=[X_dot; reshape(F_linear*delta_mat, n^2, 1);];

end

function [X_dot, F_linear]= eom_hover(n, n_vars, INSECT, WK_R, WK_L, t, X, varargin)
%% Stability analysis for nominal hover trajectory
x=X(1:3);
x_dot=X(4:6);
delta_mat=reshape(X(n_vars+1:(n^2+n_vars)), n, n);

X = X(1:6);
[X_dot, R, Q_R, Q_L, Q_A, theta_B, theta_A, W, W_dot, W_R, ...
    W_R_dot, W_L, W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, ...
    f_tau, tau] = eom_QS_x(INSECT, WK_R, WK_L, t, X);
JJ = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
JJ_11 = inertia_sub_decompose_3_12(JJ);

F_linear = zeros(n, n);
[d_L_R, d_L_L, d_D_R, d_D_L]=wing_QS_aerodynamics_linearized(INSECT, W_R, W_L, W_R_dot, W_L_dot, x_dot, R, W, Q_R, Q_L);
d_F_R = d_L_R + d_D_R;
d_F_L = d_L_L + d_D_L;
F_linear(1:3, 4:6) = eye(3);
F_linear(4:6, 4:6) = JJ_11\ R * (Q_R*d_F_R + Q_L*d_F_L);
X_dot=[X_dot; reshape(F_linear*delta_mat, n^2, 1);];

end

function [X_dot, F_linear]= eom_hover_control(n, n_vars, INSECT, WK_R_des, WK_L_des, t, X, varargin)
des = varargin{1}{1};
gains = varargin{1}{2};
wt = varargin{1}{3};
bound_param = varargin{1}{4};
%% Stability analysis for controlled hover trajectory
x=X(1:3);
x_dot=X(4:6);
int_d_x=X(7:9);
delta_mat = reshape(X(n_vars+1:(n^2+n_vars)), n, n);
X=X(1:9);

%% Ideal uncontrolled values
[Euler_R, Euler_R_dot, Euler_R_ddot] = wing_kinematics(t,WK_R_des);
[Euler_L, Euler_L_dot, Euler_L_ddot] = wing_kinematics(t,WK_L_des);
[Q_R, Q_L, W_R, W_L, W_R_dot, W_L_dot] = wing_attitude(WK_R_des.beta, ...
    Euler_R, Euler_L, Euler_R_dot, Euler_L_dot, Euler_R_ddot, Euler_L_ddot);
[R, W] = body_attitude(t, WK_R_des.f, WK_R_des); % body
[Q_A, W_A, W_A_dot] = abdomen_attitude(t, WK_R_des.f, WK_R_des); % abdomen
[JJ_A, KK_A] = inertia_wing_sub(INSECT.m_A, INSECT.mu_A, INSECT.nu_A, INSECT.J_A, R, Q_A, X(4:6), W, W_A);
f_abd_des = -(JJ_A(1:3, 7:9)*W_A_dot + KK_A(1:3, 7:9)*W_A);

[L_R, L_L, D_R, D_L]=wing_QS_aerodynamics(INSECT, ...
    W_R, W_L, W_R_dot, W_L_dot, x_dot, R, W, Q_R, Q_L);
F_R = L_R + D_R;
F_L = L_L + D_L;
f_a_des = R*Q_R*F_R + R*Q_L*F_L;

%% Ideal controlled trajectory
d_x = zeros(3, 1); d_x_dot = zeros(3, 1);
for k=1:3
    d_x(k) = des.x_fit{k}(t) - x(k);
    d_x_dot(k) = des.x_dot_fit{k}(t) - x_dot(k);
end
pos_err = INSECT.m*(gains.Kp_pos * d_x + gains.Kd_pos * d_x_dot + gains.Ki_pos * int_d_x);
[~, WK_R, WK_L] = position_controller(pos_err, WK_R_des, WK_L_des, des, f_a_des, f_abd_des, wt, bound_param);

[X_dot, R, ~, ~, Q_A, ~, ~, W, ~, ~, ~, ~, ~, W_A, W_A_dot, ~, ~, ~, ~,...
    f_a] = eom_QS_x(INSECT, WK_R, WK_L, t, X(1:6));
[JJ_A, KK_A] = inertia_wing_sub(INSECT.m_A, INSECT.mu_A, INSECT.nu_A, INSECT.J_A, R, Q_A, X(4:6), W, W_A);
f_abd = -(JJ_A(1:3, 7:9)*W_A_dot + KK_A(1:3, 7:9)*W_A);
f_total_1 = f_a(1:3) + f_abd;

%% Controlled trajectory with perturbation
delta_mat_dot = zeros(n,n);
parfor j=1:size(delta_mat, 2)
    X_del = X + delta_mat(:, j);
    x=X_del(1:3);
    x_dot=X_del(4:6);
    d_x_del = zeros(3,1); d_x_dot_del = zeros(3,1);
    for k=1:3
        d_x_del(k) = des.x_fit{k}(t) - x(k);
        d_x_dot_del(k) = des.x_dot_fit{k}(t) - x_dot(k);
    end
    pos_err = INSECT.m*(gains.Kp_pos * d_x_del + gains.Kd_pos * d_x_dot_del + gains.Ki_pos * int_d_x);
    [~, WK_R, WK_L] = position_controller(pos_err, WK_R_des, WK_L_des, des, f_a_des, f_abd_des, wt, bound_param);

    [~, R, ~, ~, Q_A, ~, ~, W, ~, ~, ~, ~, ~, W_A, W_A_dot, ~, ~, ~, ~, ...
        f_a] = eom_QS_x(INSECT, WK_R, WK_L, t, X_del(1:6));
    [JJ_A, KK_A] = inertia_wing_sub(INSECT.m_A, INSECT.mu_A, INSECT.nu_A, INSECT.J_A, R, Q_A, X_del(4:6), W, W_A);
    f_abd = -(JJ_A(1:3, 7:9)*W_A_dot + KK_A(1:3, 7:9)*W_A);
    
    f_total_2 = f_a(1:3) + f_abd;
    del_f = f_total_2 - f_total_1;
    delta_mat_dot(:,j) = [delta_mat(4:6,j); del_f / INSECT.m; delta_mat(1:3,j);];
end

F_linear = delta_mat_dot / delta_mat;
X_dot = [X_dot; d_x; reshape(delta_mat_dot, n^2, 1);];

end

function [dang, WK_R, WK_L] = position_controller(pos_err, WK_R_des, WK_L_des, des, f_a_des, f_abd_des, wt, bound_param)
%%
    mul_phi = [des.params.df_a_1_by_dphi_m(round((3 - sign(f_a_des(1)))/2)); 0; des.params.df_a_3_by_dphi_m];
    mul_theta = [des.params.df_a_1_by_dtheta_0(round((3 - sign(f_a_des(1)))/2)); 0; des.params.df_a_3_by_dtheta_0];
    mul_theta_A = [des.params.df_a_1_by_dtheta_A_m(round((3 - sign(f_abd_des(1)))/2)); 0; ...
        des.params.df_a_3_by_dtheta_A_m(round((3 - sign(f_abd_des(3)))/2))];
    temp_A = [mul_phi(1), mul_theta(1), mul_theta_A(1);
              mul_phi(3), mul_theta(3), mul_theta_A(3);
              wt * mul_phi(1), wt * mul_theta(1), -(1-wt) * mul_theta_A(1)];
    rhs = [pos_err(1); pos_err(3); 0];
    dang = zeros(5, 1);
    dang(1:3) = temp_A \ (rhs);
    %
%     if wt == 0
%         dang(1:2) = temp_A(1:2,1:2) \ rhs(1:2);
%     else
%         % Minimum norm solution
%         temp_A = temp_A(1:2,:);
%         dang(1:3) = temp_A' * ((temp_A*temp_A') \ rhs(1:2));
%     end
    %
    dang(4) = pos_err(2) / des.params.df_a_2_by_dphi_m(2);
    idx = abs(dang) > bound_param;
    dang(idx) = bound_param * sign(dang(idx));
    dphi_m_R = dang(1) + dang(4);
    dphi_m_L = dang(1) - dang(4);
    dtheta_0 = dang(2);
    dtheta_A_m = dang(3);
    dpsi_m = dang(5);

    WK_R = WK_R_des; WK_L = WK_L_des;
    WK_R.phi_m = WK_R_des.phi_m + dphi_m_R;
    WK_L.phi_m = WK_L_des.phi_m + dphi_m_L;
    WK_R.theta_0 = WK_R_des.theta_0 + dtheta_0;
    WK_L.theta_0 = WK_L_des.theta_0 + dtheta_0;
    WK_R.psi_m = WK_R_des.psi_m + dpsi_m;
    WK_L.psi_m = WK_L_des.psi_m - dpsi_m;
    WK_R.theta_A_m = WK_R_des.theta_A_m + dtheta_A_m;
    WK_L.theta_A_m = WK_L_des.theta_A_m + dtheta_A_m;
end

function [X_dot, F_linear]= eom_hover_xR(n, n_vars, INSECT, WK_R, WK_L, t, X, varargin)
%% Stability analysis for nominal hover trajectory
x=X(1:3);
x_dot=X(13:15);
delta_mat=reshape(X(n_vars+1:(n^2+n_vars)), n, n);

X = X(1:n_vars);
[X_dot, R, Q_R, Q_L, Q_A, ~, W, W_R, ...
    W_R_dot, W_L, W_L_dot, W_A, W_A_dot] = eom_QS_xR(INSECT, WK_R, WK_L, t, X);
xi_1_dot = X_dot(13:18);

delta_mat_dot = zeros(n, n);
xi_dot = [X_dot(13:18); W_R_dot; W_L_dot; W_A_dot;];
JJ = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);

parfor j=1:size(delta_mat, 2)
    dx = delta_mat(1:3, j); dR = delta_mat(4:6, j); dx_dot = delta_mat(7:9, j); dW = delta_mat(10:12, j);
    X_del = [X(1:3)+dx; reshape(R*expmhat(dR), 9, 1); X(13:18) + [dx_dot; dW]];
%     dx = delta_mat(1:3, j); dR = delta_mat(4, j); dx_dot = delta_mat(5:7, j); dW = delta_mat(8, j);
%     X_del = [X(1:3)+dx; reshape(R*expmhat([0, dR, 0]), 9, 1); X(13:18) + [dx_dot; 0; dW; 0]];
    %
    [X_dot_del, R_del, Q_R_del, Q_L_del, Q_A_del, ~, W_del, W_R_del, ...
        W_R_dot_del, W_L_del, W_L_dot_del, W_A_del, W_A_dot_del] = eom_QS_xR(INSECT, WK_R, WK_L, t, X_del);
    x_dot_del = X_del(13:15);
    xi_dot_del = [X_dot_del(13:18); W_R_dot_del; W_L_dot_del; W_A_dot_del;];
    JJ_del = inertia(INSECT, R_del, Q_R_del, Q_L_del, Q_A_del, x_dot_del, W_del, W_R_del, W_L_del, W_A_del);
    dxi_ddot = ((JJ_del*xi_dot_del - JJ*xi_dot) - (JJ_del - JJ)*xi_dot);
    dxi_ddot = JJ(1:6,1:6) \ dxi_ddot(1:6);
% %     dxi_ddot = JJ \ dxi_ddot;
% %     dxi_ddot = dxi_ddot(1:6);
    %
    d_mat = zeros(n, 1);
    d_mat(1:3) = dx_dot;
    d_mat(4:6) = -hat(W)*dR + dW;
    d_mat(7:12) = dxi_ddot(1:6);
%     d_mat(1:3) = dx_dot;
%     d_mat(4) = dW;
%     d_mat(5:7) = dxi_ddot(1:3);
%     d_mat(8) = dxi_ddot(5);
    delta_mat_dot(:, j) = d_mat;
end
F_linear = delta_mat_dot / delta_mat;
X_dot = [X_dot; reshape(delta_mat_dot, n^2, 1);];

end

function [X_dot, F_linear] = eom_hover_xR_control(n, n_vars, INSECT, WK_R_des, WK_L_des, t, X, varargin)
%%
des = varargin{1}{1};
gains = varargin{1}{2};
wt = varargin{1}{3};
bound_param = varargin{1}{4};

x=X(1:3);
R=reshape(X(4:12),3,3);
x_dot=X(13:15);
W=X(16:18);
int_d_x=X(19:21);
int_att=X(22:24);

delta_mat=reshape(X(n_vars+1:(n^2+n_vars)), n, n);
X = X(1:n_vars);

%% Ideal uncontrolled values
d_x = zeros(3, 1); des_x_dot = zeros(3, 1);
des_R = zeros(3, 3); des_W = zeros(3, 1);
for k=1:3
    d_x(k) = des.x_fit{k}(t) - x(k);
    des_x_dot(k) = des.x_dot_fit{k}(t);
    des_W(k) = des.W_fit{k}(t);
    for j=1:3
        des_R(k, j) = des.R_fit{k,j}(t);
    end
end

[Euler_R, Euler_R_dot, Euler_R_ddot] = wing_kinematics(t,WK_R_des);
[Euler_L, Euler_L_dot, Euler_L_ddot] = wing_kinematics(t,WK_L_des);
[Q_R, Q_L, W_R, W_L, W_R_dot, W_L_dot] = wing_attitude(WK_R_des.beta, ...
    Euler_R, Euler_L, Euler_R_dot, Euler_L_dot, Euler_R_ddot, Euler_L_ddot);
[Q_A, W_A, W_A_dot] = abdomen_attitude(t, WK_R_des.f, WK_R_des); % abdomen

[L_R, L_L, D_R, D_L, M_R, M_L]=wing_QS_aerodynamics(INSECT, ...
    W_R, W_L, W_R_dot, W_L_dot, des_x_dot, des_R, des_W, Q_R, Q_L);
F_R=L_R+D_R;
F_L=L_L+D_L;
M_A=zeros(3,1);
f_a=[des_R*Q_R*F_R + des_R*Q_L*F_L;
    hat(INSECT.mu_R)*Q_R*F_R + hat(INSECT.mu_L)*Q_L*F_L;
    M_R; M_L; M_A];
f_a_1=f_a(1:6);
f_a_2=f_a(7:15);

[~, dU]=potential(INSECT,x,des_R,Q_R,Q_L,Q_A);
f_g=-dU;
f_g_1=f_g(1:6);
f_g_2=f_g(7:15);
C=[zeros(3,9); -Q_R -Q_L -Q_A];

tmp_f_des = f_a_1+f_g_1 - C*(f_a_2+f_g_2);

%% Ideal controlled trajectory
JJ_xR = get_cross_terms(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
[~, ~, ~, u_control, int_att_dot] = controller(t, des, gains, X, JJ_xR, ...
    WK_R_des, WK_L_des, wt, bound_param);

[X_dot, R, Q_R, Q_L, Q_A, theta_A, W] = ...
    eom_QS_xR_ideal(INSECT, WK_R_des, WK_L_des, t, X(1:18), u_control, tmp_f_des);
xi_dot = [X_dot(13:18)];

% f_total_1 = f_a(1:3);

%% Controlled trajectory with perturbation
delta_mat_dot = zeros(n,n);
parfor j=1:size(delta_mat, 2)
    dx = delta_mat(1:3, j); dR = delta_mat(4:6, j);
    dx_dot = delta_mat(7:9, j); dW = delta_mat(10:12, j);
    X_del = [X(1:3)+dx; reshape(R*expmhat(dR), 9, 1); X(13:24) + delta_mat(7:18, j)];
    
    JJ_xR_del = get_cross_terms(INSECT, R*expmhat(dR), Q_R, Q_L, Q_A, ...
        x_dot+dx_dot, W+dW, W_R, W_L, W_A);
    [~, ~, ~, u_control, int_att_dot_del] = controller(t, des, gains, X_del, JJ_xR, ...
        WK_R_des, WK_L_des, wt, bound_param);
    
    [X_dot_del] = ...
    eom_QS_xR_ideal(INSECT, WK_R_des, WK_L_des, t, X_del(1:18), u_control, tmp_f_des);

    xi_dot_del = [X_dot_del(13:18)];
%     dxi_ddot = ((JJ_del*xi_dot_del - JJ*xi_dot) - (JJ_del - JJ)*xi_dot);
%     dxi_ddot = JJ(1:6,1:6) \ dxi_ddot(1:6);

    dxi_ddot = JJ_xR \ ((JJ_xR_del*xi_dot_del(1:6) - JJ_xR*xi_dot(1:6)) - (JJ_xR_del - JJ_xR)*xi_dot(1:6));
    %
    d_mat = zeros(n, 1);
    d_mat(1:3) = dx_dot;
    d_mat(4:6) = -hat(W)*dR + dW;
    d_mat(7:12) = dxi_ddot(1:6);
    d_mat(13:15) = dx;
    d_mat(16:18) = int_att_dot_del - int_att_dot;
    delta_mat_dot(:, j) = d_mat;
end

F_linear = delta_mat_dot / delta_mat;
X_dot = [X_dot(1:18); d_x; int_att_dot; reshape(delta_mat_dot, n^2, 1);];

end

function [dang, WK_R, WK_L, u_control, int_att_dot] = controller(t, des, gains, X,...
    JJ_xR, WK_R_des, WK_L_des, wt, bound_param)
%%
    x=X(1:3);
    R=reshape(X(4:12),3,3);
    x_dot=X(13:15);
    W=X(16:18);
    int_d_x=X(19:21);
    int_att=X(22:24);
    
    d_x = zeros(3, 1); d_x_dot = zeros(3, 1);
    des_R = zeros(3, 3); des_W = zeros(3, 1);
    for k=1:3
        d_x(k) = des.x_fit{k}(t) - x(k);
        d_x_dot(k) = des.x_dot_fit{k}(t) - x_dot(k);
        des_W(k) = des.W_fit{k}(t);
        for j=1:3
            des_R(k, j) = des.R_fit{k,j}(t);
        end
    end
    e_R = 0.5*vee(des_R'*R - R'*des_R);
    e_Om = W - R'*des_R*des_W;
    int_att_dot = e_Om + gains.cI * e_R;

    pos_err = (gains.Kp_pos * d_x + gains.Kd_pos * d_x_dot + gains.Ki_pos * int_d_x);
    att_err = (- gains.KR*e_R - gains.KOm*e_Om - gains.KI*int_att);
    err_xR = JJ_xR*[pos_err; att_err];
    
    u_control = err_xR;% + tmp_mat - tmp_des;
%     u_control = err_xR - tmp_des;

    dang = zeros(6, 1);
    WK_R = WK_R_des; WK_L = WK_L_des;
%     WK_R.phi_m = WK_R_des.phi_m + dphi_m_R;
%     WK_L.phi_m = WK_L_des.phi_m + dphi_m_L;
%     WK_R.theta_0 = WK_R_des.theta_0 + dtheta_0;
%     WK_L.theta_0 = WK_L_des.theta_0 + dtheta_0;
%     WK_R.psi_m = WK_R_des.psi_m + dpsi_m;
%     WK_L.psi_m = WK_L_des.psi_m - dpsi_m;
%     WK_R.theta_A_m = WK_R_des.theta_A_m + dtheta_A_m;
%     WK_L.theta_A_m = WK_L_des.theta_A_m + dtheta_A_m;
end

function [JJ_xR] = get_cross_terms(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A)
    JJ = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
    [JJ_11, ~, JJ_21, ~] = inertia_sub_decompose_6_9(JJ);
    C=[zeros(3,9);
        -Q_R -Q_L -Q_A];
    JJ_xR = (JJ_11 - C*JJ_21);

%     LL = KK - 0.5*KK';
%     co_ad=blkdiag(zeros(3,3), -hat(W), -hat(W_R), -hat(W_L), -hat(W_A));
%     [LL_11, ~, LL_21, ~] = inertia_sub_decompose_6_9(LL);
%     [co_ad_11, ~, ~, co_ad_22] = inertia_sub_decompose_6_9(co_ad);
%     xi_1 = [x_dot; W];
%     tmp_mat = -(co_ad_11*JJ_11-C*co_ad_22*JJ_21)*xi_1 + (LL_11-C*LL_21)*xi_1;
end
