function floquet_stability_opt
% Study the stability of a periodic trajecoty using floquet theory.

evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
filename='floquet_stability_opt';

%% Linearized dynamics
bool_sort_mus=false; % false if n ~= 6

load('sim_QS_xR_hover_control_opt.mat', 'INSECT', 'WK', 'des', ...
    'control_net', 'Weights', 'N_dang', 'm', 'X0_pert');
N_period = 2;
N_single = 500;
N = 1 + N_period*N_single;
N_char = 1 + (N_period-1)*N_single;
T=N_period/WK.f;
t=linspace(0,T,N);
dt = t(2) - t(1);
des_X0 = des.X0;
X0 = des_X0;
epsilon = 1e-6; % 1e-10 stable, 1e-12 unstable
del = 1e0;
% epsilon = 1e0; % 1e-10 stable, 1e-12 unstable
% del = 0;
rng(1);
dx = 2*rand(1,3)-1; dx = rand(1) * dx / norm(dx);
dtheta = 2*rand(1,3)-1; dtheta = rand(1) * dtheta / norm(dtheta);
dx_dot = 2*rand(1,3)-1; dx_dot = rand(1) * dx_dot / norm(dx_dot);
domega = 2*rand(1,3)-1; domega = rand(1) * domega / norm(domega);
dXi = del * Weights.PerturbVariables .* [dx, dtheta, dx_dot, domega];
e1 = [1 0 0]'; e2 = [0 1 0]'; e3 = [0 0 1]';
des_R0 = reshape(des_X0(4:12), 3, 3);
X0 = [des_X0(1:3)'+dXi(1:3),...
    reshape(des_R0*expmhat(dXi(6)*e3)*expmhat(dXi(5)*e2)*expmhat(dXi(4)*e1),1,9),...
    des_X0(13:18)' + dXi(7:12)]';

n = 12; n_vars = 18; % for (x,R) control with @eom_hover_xR_control
tic;
[delta_mat, F_linear] = sim_pert(@eom_hover_xR_control, n, n_vars, INSECT, WK, X0, N, t, epsilon, ...
    control_net, Weights, des_X0, N_dang, m, N_single, N_period);
% [delta_mat, F_linear] = sim_pert_non(@eom_hover_xR_control_non, n, n_vars, INSECT, WK, X0, N, t, epsilon, ...
%     control_net, Weights, des_X0, N_dang, m, N_single, N_period);
time_taken = toc;

delta_mag = zeros(n, N);
X_ref = zeros(N, 18);
idx = 1:(1+N_single);
[~, X_ref(idx, :)]=ode45(@(t,X) eom_QS_xR(INSECT, WK, WK, t,X), ...
    t(idx), des_X0, odeset('AbsTol',1e-6,'RelTol',1e-6));
idx = (1+N_single):N;
X_ref(idx, :) = X_ref(mod(idx-1, N_single)+1, :);

for k=1:N    
    [~, ~, ~, ~, ~,~, ~, ~, ~, ~,~, ~, ~, ~, ~, ~,~, ~, ~, ~, ~, ~, ...
        Euler_R_dot(:,k)]= eom_QS_xR(INSECT, WK, WK, t(k), X_ref(k,:)');
    delta_mag(:, k) = vecnorm((1./Weights.PerturbVariables)' .* delta_mat(:, :, k), 2, 1);
%     delta_mag(:, k) = max(abs((1./Weights.PerturbVariables)' .* delta_mat(:, :, k)), [], 1);
end

B = zeros(n, n, N_char);
mu_vals = zeros(n, N_char);
delta_ratio = zeros(n, N_char);
start_ix = N - N_char - N_single + 1;
for i=start_ix:(start_ix+N_char-1)
    j = i-start_ix+1;
    d_F = (F_linear(:, :, i) - F_linear(:, :, i+N_single))./ F_linear(:, :, i);
    %
    max_F = min(max(F_linear(:,:,i), [], 'all'), max(F_linear(:,:,i+N_single), [], 'all'));
    d_F(max(abs(F_linear(:, :, i+N_single)/max_F), abs(F_linear(:, :, i)/max_F)) < 1e-4) = 0; % Insignificant values
    idx_small = max(abs(F_linear(:, :, i+N_single)/max_F), abs(F_linear(:, :, i)/max_F)) < 1e-1;
    d_F(d_F(idx_small) < 0.25) = 0; % Small values
    %
    if(~all(abs(d_F) < 5e-2, 'all'))
        disp(i)
        disp(d_F)
    end
    B(:, :, j) = delta_mat(:, :, i) \ delta_mat(:, :, i+N_single);
    mu_vals(:,j) = log(abs(eig(B(:,:,j))))*WK.f;
    delta_ratio(:,j) = delta_mag(:, i) .\ delta_mag(:, i+N_single);
end
B = B(:, :, round(end/2));
B = delta_mat(:, :, 1) \ delta_mat(:, :, 1+N_single);
[e_vecs, rhos] = eig(B);
mus = log(abs(diag(rhos))) * WK.f;

h_char_exp = figure;
hold on;
for i=1:n
    plot(t(1:N_char)*WK.f, mu_vals(i, :));
end
patch_downstroke(h_char_exp,t(1:N_char)*WK.f,Euler_R_dot(:, 1:N_char));
plot(t(1:N_char)*WK.f, zeros(1,N_char), 'k--', 'LineWidth', 2);
axis('tight');
xlabel('$t/T$','interpreter','latex');
ylabel('$\mu$','interpreter','latex');

h_mag_rat = figure;
hold on;
for i=1:n
    plot(t(1:N_char)*WK.f, delta_ratio(i, :));
end
patch_downstroke(h_mag_rat,t(1:N_char)*WK.f,Euler_R_dot(:, 1:N_char));
plot(t(1:N_char)*WK.f, ones(1,N_char), 'k--', 'LineWidth', 2);
axis('tight');
xlabel('$t/T$','interpreter','latex');
ylabel('$\|\frac{\delta(t+T)}{\delta(t)}\|$','interpreter','latex');

h_mag = figure;
hold on;
for i=1:n
    plot(t(1:N)*WK.f, delta_mag(i, :)/delta_mag(i, 1));
end
patch_downstroke(h_mag,t(1:N)*WK.f,Euler_R_dot(:, 1:N));
plot(t(1:N)*WK.f, ones(1,N), 'k--', 'LineWidth', 2);
axis('tight');
xlabel('$t/T$','interpreter','latex');
ylabel('$\|\frac{\delta(t)}{\delta(0)}\|$','interpreter','latex');

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

function [delta_mat, F_linear] = sim_pert(eom, n, n_vars, INSECT, WK, X0, N, t, epsilon, varargin)
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
F_linear = zeros(n, n, N);

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
        
        [~, X(idx,:)] = ode45(@(t,X) eom(n, n_vars, INSECT, WK, WK, t, X, ...
            t0, dang0, param_m, dang0_mat, param_m_mat), t(idx), X0);
        parfor k=idx
            [~, F_linear(:,:,k)] = eom(n, n_vars, INSECT, WK, WK, t(k), X(k, :)', ...
                t0, dang0, param_m, dang0_mat, param_m_mat);
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

function [delta_mat, F_linear] = sim_pert_non(eom, n, n_vars, INSECT, WK, X0, N, t, epsilon, varargin)
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
X_mat0 = zeros(n_vars, n);
for j=1:n
    dR = delta0(4:6, j);    
    X_mat0(:, j) = [X0(1:3)+delta0(1:3, j); reshape(reshape(X0(4:12),3,3)*expmhat(dR), 9, 1); ...
        X0(13:18) + delta0(7:12, j)];
end
X0 = [X0(1:n_vars); reshape(X_mat0, n_vars*n, 1);];
t0 = t(1); dang0 = zeros(N_p,1); param_mat = zeros(m*N_p, n); dang0_mat = zeros(N_p, n);
X = zeros(N, length(X0));
delta_mat = zeros(n, n, N);
F_linear = zeros(n, n, N);

for period=1:N_period
    idx_con = (1+(period-1)*N_single):(1+period*N_single);
    X_mat0 = reshape(X0((n_vars+1):(n_vars*n+n_vars)), n_vars, n);
    param = sim(control_net, (1 ./ Weights.PerturbVariables)' .* get_error(X0, des_X0));
    for j=1:n
        X_del0 = X_mat0(:, j);
        param_mat(:, j) = sim(control_net, (1 ./ Weights.PerturbVariables)' .* get_error(X_del0, des_X0));
    end
    for con=1:m
        param_idx = (1+(con-1)*N_p):(con*N_p);
        idx = idx_con((1+(con-1)*N_per_iter):(1+con*N_per_iter));
        param_m = param(param_idx);
        param_m_mat = param_mat(param_idx, :);
        
        [~, X(idx,:)] = ode45(@(t,X) eom(n, n_vars, INSECT, WK, WK, t, X, ...
            t0, dang0, param_m, dang0_mat, param_m_mat), t(idx), X0);
        parfor k=idx
            X_mat = reshape(X(k,n_vars+1:(n*n_vars+n_vars))', n_vars, n);
            Xd = X(k,1:n_vars)';
            for j=1:n
                delta_mat(:,j,k) = get_error(X_mat(:, j), Xd);
            end
        end
        X0 = X(idx(end),:)';
        dang0 = linear_func(param_m, t(idx(end)), t0, dang0);
        for j=1:n
            dang0_mat(:,j) = linear_func(param_m_mat(:,j), t(idx(end)), t0, dang0_mat(:,j));
        end
        t0 = t(idx(end));
    end
end

end

function [X_dot, F_linear] = eom_hover_xR_control(n, n_vars, INSECT, WK_R_des, WK_L_des, t, X, varargin)
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
% int_d_x=X(19:21);

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
% [X_dot, R, Q_R, Q_L, Q_A, ~, W, W_R, ~, W_L, ~, W_A] = ...
X_dot = eom_QS_xR(INSECT, WK_R_con, WK_L_con, t, X(1:18));
xi_dot = X_dot(13:18);
JJ_xR = get_cross_terms(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
% JJ = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);

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
    
%     JJ_del = inertia(INSECT, R, Q_R, Q_L, Q_A, x_dot, W, W_R, W_L, W_A);
%     dxi_ddot = ((JJ_del*xi_dot_del - JJ*xi_dot) - (JJ_del - JJ)*xi_dot);
%     dxi_ddot = JJ(1:6,1:6) \ dxi_ddot(1:6);

    dxi_ddot = JJ_xR \ ((JJ_xR_del*xi_dot_del(1:6) - JJ_xR*xi_dot(1:6)) -...
                        (JJ_xR_del - JJ_xR)*xi_dot(1:6));
    %
    d_mat = zeros(n, 1);
    d_mat(1:3) = dx_dot;
    d_mat(4:6) = -hat(W)*dR + dW;
    d_mat(7:12) = dxi_ddot(1:6);
%     d_mat(13:15) = dx;
    delta_mat_dot(:, j) = d_mat;
end

F_linear = delta_mat_dot / delta_mat;
X_dot = [X_dot(1:n_vars); reshape(delta_mat_dot, n^2, 1);];

end

function [X_dot] = eom_hover_xR_control_non(n, n_vars, INSECT, WK_R_des, WK_L_des, t, X, varargin)
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
% int_d_x=X(19:21);

X_mat=reshape(X(n_vars+1:(n*n_vars+n_vars)), n_vars, n);
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
% [X_dot, R, Q_R, Q_L, Q_A, ~, W, W_R, ~, W_L, ~, W_A] = ...
X_dot = eom_QS_xR(INSECT, WK_R_con, WK_L_con, t, X(1:18));

%% Controlled trajectory with perturbation
X_mat_dot = zeros(n_vars,n);
parfor j=1:n
    param = param_mat(:, j);
    dang = linear_func(param, t, t0, dang0_mat(:, j));
    [WK_R_del, WK_L_del] = get_WK(WK_R_des, WK_L_des, dang);
    
    X_del = X_mat(:, j);
    X_mat_dot(:, j) = eom_QS_xR(INSECT, WK_R_del, WK_L_del, t, X_del(1:18));
end

X_dot = [X_dot(1:n_vars); reshape(X_mat_dot, n*n_vars, 1);];

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

function dang = linear_func(param, t, varargin)
t0 = varargin{1};
dang0 = varargin{2};
dang = dang0 + param*(t - t0);
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
