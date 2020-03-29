function floquet_stability
% Study the stability of a periodic trajecoty using floquet theory.

evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
filename='floquet_stability';

%% Linearized dynamics
load('sim_QS_x_hover.mat', 'INSECT', 'WK', 'X0');
% INSECT.scale=1e-2;
% INSECT.name='MONARCH';
% WK.ab_type='varying';
% WK.bo_type='varying';
% load('sim_QS_x_hover_hawkmoth.mat');
% INSECT.name='NOT_MONARCH';

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
% n = 12; n_vars = 6; % attitude stability with @eom_hover_attitude
% [delta_mat, F_linear] = sim_pert(@eom_hover_attitude, n, n_vars, INSECT, WK, X0, N, t, epsilon);

load('sim_QS_x_hover_control.mat', 'des', 'gains', 'wt', 'bound_param');
int_d_x0 = zeros(3,1);
X0 = [X0; int_d_x0;];
n = 9; n_vars = 9; % for controlled hover with @eom_hover_control
[delta_mat, F_linear] = sim_pert(@eom_hover_control, n, n_vars, INSECT, WK, X0, N, t, epsilon, des, gains, wt, bound_param);

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

%% Study for various insects
% stop_idx = N;
% N_sims = 10;
% conv_rate_osc = zeros(N_sims, 6);
% var_name_to_save = 'conv_rate_osc';
% for c_ix=4:6 % Column index for perturbation direction
% %     delta_g_mag = vecnorm(reshape(delta_mat(1:3, c_ix, :), 3, N));
% %     delta_xi_mag = vecnorm(reshape(delta_mat(4:6, c_ix, :), 3, N));
%     conv_rate_osc(:, c_ix) = repmat(mus(c_ix), [N_sims, 1]);
% end
% % save('sim_QS_x_hover_conv_rate', var_name_to_save, '-append');

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

function [X_dot, F_linear]= eom_hover_attitude(n, n_vars, INSECT, WK_R, WK_L, t, X, varargin)
%% Stability analysis for nominal hover trajectory
x=X(1:3);
x_dot=X(4:6);
delta_mat=reshape(X(n_vars+1:(n^2+n_vars)), n, n);

X = X(1:6);
[X_dot, R, Q_R, Q_L, Q_A, theta_B, theta_A, W, W_dot, W_R, ...
    W_R_dot, W_L, W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, ...
    f_tau, tau] = eom_QS_x(INSECT, WK_R, WK_L, t, X);

delta_mat_dot = zeros(n, n);
x_ddot = X_dot(4:6);
xi_dot = [x_ddot; W_dot; W_R_dot; W_L_dot; W_A_dot;];
[JJ, euler_rhs] = EL_equation_terms(INSECT, x, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A, W_R_dot, W_L_dot, tau, f_tau);
parfor j=1:size(delta_mat, 2)
    dx = delta_mat(1:3, j); dx_dot = delta_mat(4:6, j); dR = delta_mat(7:9, j); dW = delta_mat(10:12, j);
    [JJ_new, euler_rhs_new] = EL_equation_terms(INSECT, x+dx, R*expmhat(dR), Q_R, Q_L, Q_A, x_dot+dx_dot, W+dW, W_R, W_L, W_A, W_R_dot, W_L_dot, tau, f_tau);
    dx_ddot = ((euler_rhs_new - euler_rhs) - (JJ_new - JJ)*xi_dot);
    dx_ddot = JJ(1:6, 1:6) \ dx_ddot(1:6);
    d_mat = zeros(n, 1);
    d_mat(1:3) = dx_dot;
    d_mat(4:6) = dx_ddot(1:3);
    d_mat(7:9) = -hat(W)*dR + dW;
    d_mat(10:12) = dx_ddot(4:6);
    delta_mat_dot(:, j) = d_mat;
end
F_linear = delta_mat_dot / delta_mat;
X_dot = [X_dot; reshape(delta_mat_dot, n^2, 1);];

end

function [JJ, euler_rhs] = EL_equation_terms(INSECT, x, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A, W_R_dot, W_L_dot, tau, f_tau)
%%
[L_R, L_L, D_R, D_L, M_R, M_L, ...
    F_rot_R, F_rot_L, M_rot_R, M_rot_L]=wing_QS_aerodynamics(INSECT, W_R, W_L, W_R_dot, W_L_dot, x_dot, R, W, Q_R, Q_L);
F_R=L_R+D_R+F_rot_R;
F_L=L_L+D_L+F_rot_L;
M_R=M_R+M_rot_R;
M_L=M_L+M_rot_L;
M_A=zeros(3,1);

f_a=[R*Q_R*F_R + R*Q_L*F_L; hat(INSECT.mu_R)*Q_R*F_R + hat(INSECT.mu_L)*Q_L*F_L;
    M_R; M_L; M_A];
[~, dU]=potential(INSECT,x,R,Q_R,Q_L,Q_A);
f_g=-dU;

xi = [x_dot; W; W_R; W_L; W_A];
[JJ, KK] = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
LL = KK - 0.5*KK';
co_ad=blkdiag(zeros(3,3), -hat(W), -hat(W_R), -hat(W_L), -hat(W_A));

f_tau = [zeros(3,1); f_tau(4:6); Q_R'*tau(4:6); Q_L'*tau(7:9); Q_A'*tau(10:12);];

euler_rhs = co_ad*JJ*xi - LL*xi + f_a + f_g + f_tau;

end
