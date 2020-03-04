function floquet_stability
% Study the stability of a periodic trajecoty using floquet theory.

evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');

%% Linearized dynamics
load('sim_QS_x_hover_control.mat', 'gains');
load('sim_QS_x_hover.mat');
filename='floquet_stability';

% INSECT.scale=1e-2;
% INSECT.name='MONARCH';
% WK.ab_type='varying';
% WK.bo_type='varying';
% load('sim_QS_x_hover_hawkmoth.mat');
% INSECT.name='NOT_MONARCH';

N = 1001;
N_period = 10;
T=N_period/WK.f;
ix_d = (N-1)/N_period;
t=linspace(0,T,N);
dt = T/(N-1);
epsilon = 1e-8;

% n is the number of perturbation states 
% n = 3; % for nominal hover with @eom_hover_vel if position is not periodic
% [delta_mat, F_linear] = sim_pert(@eom_hover_vel, n, INSECT, WK, X0, N, t, epsilon);
% n = 6; % for nominal hover with @eom_hover
% [delta_mat, F_linear] = sim_pert(@eom_hover, n, INSECT, WK, X0, N, t, epsilon);
% n = 9; % for controlled hover with @eom_hover_control
% [delta_mat, F_linear] = sim_pert(@eom_hover_control, n, INSECT, WK, X0, N, t, epsilon, gains);
n = 12; % attitude stability with @eom_hover_attitude
[delta_mat, F_linear] = sim_pert(@eom_hover_attitude, n, INSECT, WK, X0, N, t, epsilon);

B = zeros(n, n, 1+ix_d);
start_ix = max(1, round((N_period-2)/N_period * N));
for i=start_ix:(start_ix+ix_d)
    j = i-start_ix+1;
    d_F = (F_linear(:, :, i) - F_linear(:, :, i+ix_d))./ F_linear(:, :, i);
    d_F(max(abs(F_linear(:, :, i+ix_d)), abs(F_linear(:, :, i))) < 1e-10) = 0;
    if(~all(abs(d_F) < 5e-2, [1, 2]))
        disp(d_F)
    end
    B(:, :, j) = delta_mat(:, :, i) \ delta_mat(:, :, i+ix_d);
end
B = B(:, :, round(end/2));
[e_vecs, rhos] = eig(B);
mus = log(diag(rhos)) * WK.f;

char_soln_mat = zeros(n, n, N);
per_val_mat = zeros(n, n, N);

for i=1:N
    char_soln_mat(:, :, i) = delta_mat(:, :, i) * e_vecs;
    Y_0 = diag(exp(mus*t(i)));
    per_val_mat(:, :, i) = char_soln_mat(:, :, i) / Y_0;
end

stop_idx = N;

%% Study for various insects
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

function [delta_mat, F_linear] = sim_pert(eom, n, INSECT, WK, X0, N, t, epsilon, varargin)
%%
delta0 = diag(rand(n, 1))*epsilon;
X0 = [X0(1:6); reshape(delta0, n^2, 1);];

[t,X]=ode45(@(t,X) eom(n, INSECT, WK, WK, t, X, varargin), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
F_linear=zeros(n, n, N);
for k=1:N
    [~, F_linear(:, :, k)] = eom(n, INSECT, WK, WK, t(k), X(k, :)', varargin);
end

delta_mat=reshape(X(:,7:(n^2+6))', n, n, N);

end

function [X_dot F_linear R Q_R Q_L Q_A theta_B theta_A W W_dot W_R W_R_dot W_L W_L_dot W_A W_A_dot F_R F_L M_R M_L f_a f_g f_tau tau]= eom_hover_vel(n, INSECT, WK_R, WK_L, t, X, varargin)
%% Stability analysis for only velocity perturbations
x=X(1:3);
x_dot=X(4:6);
delta_mat=reshape(X(7:(n^2+6)), n, n);

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

function [X_dot F_linear R Q_R Q_L Q_A theta_B theta_A W W_dot W_R W_R_dot W_L W_L_dot W_A W_A_dot F_R F_L M_R M_L f_a f_g f_tau tau]= eom_hover(n, INSECT, WK_R, WK_L, t, X, varargin)
%% Stability analysis for nominal hover trajectory
x=X(1:3);
x_dot=X(4:6);
delta_mat=reshape(X(7:(n^2+6)), n, n);

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

function [X_dot F_linear R Q_R Q_L Q_A theta_B theta_A W W_dot W_R W_R_dot W_L W_L_dot W_A W_A_dot F_R F_L M_R M_L f_a f_g f_tau tau]= eom_hover_control(n, INSECT, WK_R, WK_L, t, X, varargin)
%% Stability analysis for controlled hover trajectory
x=X(1:3);
x_dot=X(4:6);
delta_mat=reshape(X(7:(n^2+6)), n, n);
gains = varargin{1}{1};

X = X(1:6);
[X_dot, R, Q_R, Q_L, Q_A, theta_B, theta_A, W, W_dot, W_R, ...
    W_R_dot, W_L, W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, ...
    f_tau, tau] = eom_QS_x(INSECT, WK_R, WK_L, t, X);
JJ = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
JJ_11 = inertia_sub_decompose_3_12(JJ);

F_linear = zeros(n, n);
[d_L_R, d_L_L, d_D_R, d_D_L]=wing_QS_aerodynamics_linearized(INSECT, W_R, W_L, W_R_dot, W_L_dot, x_dot, R, W, Q_R, Q_L);
d_F_R=d_L_R+d_D_R;
d_F_L=d_L_L+d_D_L;
F_linear(1:3, 4:6) = eye(3);
F_linear(4:6, 7:9) = eye(3);
F_linear(7:9, 1:3) = -gains.Ki_pos * eye(3);
F_linear(7:9, 4:6) = -gains.Kp_pos * eye(3);
F_linear(7:9, 7:9) = JJ_11\ R * (Q_R*d_F_R + Q_L*d_F_L) - gains.Kd_pos * eye(3);
X_dot=[X_dot; reshape(F_linear*delta_mat, n^2, 1);];

end

function [X_dot F_linear R Q_R Q_L Q_A theta_B theta_A W W_dot W_R W_R_dot W_L W_L_dot W_A W_A_dot F_R F_L M_R M_L f_a f_g f_tau tau]= eom_hover_attitude(n, INSECT, WK_R, WK_L, t, X, varargin)
%% Stability analysis for nominal hover trajectory
x=X(1:3);
x_dot=X(4:6);
delta_mat=reshape(X(7:(n^2+6)), n, n);

X = X(1:6);
[X_dot, R, Q_R, Q_L, Q_A, theta_B, theta_A, W, W_dot, W_R, ...
    W_R_dot, W_L, W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, ...
    f_tau, tau] = eom_QS_x(INSECT, WK_R, WK_L, t, X);

delta_mat_dot = zeros(n, n);
x_ddot = X_dot(4:6);
xi_dot = [x_ddot; W_dot; W_R_dot; W_L_dot; W_A_dot;];
parfor j=1:size(delta_mat, 2)
    [JJ, euler_rhs] = EL_equation_terms(INSECT, x, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A, W_R_dot, W_L_dot, tau, f_tau);
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

% [JJ_11, JJ_12, JJ_21, JJ_22] = inertia_sub_decompose_3_12(JJ);
% [LL_11, LL_12, LL_21, LL_22] = inertia_sub_decompose_3_12(LL);
% [~, ~, ~, co_ad_22] = inertia_sub_decompose_3_12(co_ad);
% 
% xi_1_dot = JJ_11\( -JJ_12*xi_2_dot -LL_11*xi_1 - LL_12*xi_2 + f_a_1 + f_g_1);
% f_tau_2 = JJ_21*xi_1_dot + JJ_22*xi_2_dot - co_ad_22*(JJ_21*xi_1 + JJ_22*xi_2) ...
%     + LL_21*xi_1 + LL_22*xi_2 - f_a_2 - f_g_2;
% f_tau = [zeros(3,1); -tau(4:6)-tau(7:9)-tau(10:12); Q_R'*tau(4:6); Q_L'*tau(7:9); Q_A'*tau(10:12);];
f_tau = [zeros(3,1); f_tau(4:6); Q_R'*tau(4:6); Q_L'*tau(7:9); Q_A'*tau(10:12);];

euler_rhs = co_ad*JJ*xi - LL*xi + f_a + f_g + f_tau;

end
