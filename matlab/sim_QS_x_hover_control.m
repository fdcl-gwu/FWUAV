function sim_QS_x_hover_control
% simulate the thorax (x) trajectory along with a controller for
% given thorax attiude, wing kinematics, abdomen attitude.

evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
des = load('sim_QS_x_hover.mat', 'INSECT', 't', 'x', 'x_dot',...
    'f_tau', 'x0', 'x_dot0', 'WK');

filename = 'sim_QS_x_hover_control';
INSECT = des.INSECT;
WK = des.WK;
des.x_fit = cell(3, 1); des.x_dot_fit = cell(3, 1);
% 'fourier8', 'cubicinterp'
for i=1:3
    des.x_fit{i} = fit(des.t, des.x(i, :)', 'fourier8');
    des.x_dot_fit{i} = fit(des.t, des.x_dot(i, :)', 'fourier8');
end

% Values obtained from parametric study
load('parametric_study.mat', 'params');
des.params = params;
% des.params.df_a_2_by_dpsi_m = [1.663e-3, -1.444e-3]; % dpsi_m_R > 0, dpsi_m_L < 0; can use either of these values
% des.params.df_r_1_by_dphi_m = [1.7e-4, -2.74e-4]; % dphi_m_R > 0, dphi_m_L > 0
% des.params.df_r_3_by_dphi_m = [1.97e-3, -4.31e-3]; % dphi_m_R > 0, dphi_m_L > 0
%     des.f_abd(:, i) = -(JJ_A(1:3, 7:9)*W_A_dot + KK_A(1:3, 7:9)*W_A);

N = 1001;
N_period = 15;
err_bound = 1e-4; % Convergence criterion is a f(N, N_period, WK)
N_single = round((N-1)/N_period);
T = N_period/WK.f;
t = linspace(0,T,N);
bound_param = 0.1; % Parameter bound; Use 0.25?

des.x_fit_t = zeros(3, N); des.x_dot_fit_t = zeros(3, N);
for j=1:3
    des.x_fit_t(j, :) = des.x_fit{j}(t);
    des.x_dot_fit_t(j, :) = des.x_dot_fit{j}(t);
end
des = rmfield(des, {'WK', 'INSECT', 't', 'x', 'x_dot', 'f_tau'}); % 'x_fit', 'x_dot_fit'

%% Gains
pol = poly([-7.8 + 19i, -7.8 - 19i, -0.003]);
if ~all(real(roots(pol)) < 0)
    error('The chosen gains are not suitable');
end
% Optimized gs = [427.1529   15.6076  13.4983];
gains.Kp_pos = pol(3); gains.Kd_pos = pol(2); gains.Ki_pos = pol(4);

%% Single simulation
rng default;
eps1 = 1e-1/2; eps2 = 1e-1/2;
dx0 = [rand(1); 0; rand(1);]*eps1;
dx0 = [-0.2037; 0; 0.2264];
dx_dot0 = zeros(3, 1)*eps2;
x0 = des.x0 + dx0;
% x0 = [-2.5; 0 ;-2.5];
% x0 = [0; 0.8; 0];
x_dot0 = des.x_dot0 + dx_dot0;
int_d_x0 = zeros(3, 1);
X0 = [x0; x_dot0; int_d_x0];
wt = 0.1;

[err_pos, N_conv, x, x_dot, int_d_x, R, Q_R, Q_L, Q_A, theta_B, theta_A, W, W_dot, W_R, W_R_dot, W_L, ...
    W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, f_tau, tau, Euler_R, ...
    Euler_R_dot, pos_err, dang] =  simulate_control(gains, WK, INSECT, des, X0, N, ...
    N_single, N_period, t, wt, bound_param, err_bound);

%% Monte Carlo
% N_sims = 10000;
% eps = 3e0;
% wts = [0, 0.1];
% [x_pert, err_pos, N_conv] = monte_carlo(N_sims, eps, wts, gains, WK, ...
%             INSECT, des, N, N_single, N_period, t, bound_param, err_bound);

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

function [x_pert, err_pos, N_conv] = monte_carlo(N_sims, eps, wts, gains, WK, ...
            INSECT, des, N, N_single, N_period, t, bound_param, err_bound)
%%
N_wts = length(wts);
err_pos = zeros(N_sims, N_wts);
N_conv = zeros(N_sims, N_wts);
x_pert = zeros(3, N_sims);
des_x0 = des.x0;
des_x_dot0 = des.x_dot0;
x_dot0 = des_x_dot0;
int_d_x0 = zeros(3, 1);

tic;
par_pool = gcp;
nWorkers = par_pool.NumWorkers;
ticBytes(par_pool);
parfor n=1:nWorkers
    rng(n);
end
pause(1);
parfor i = 1:N_sims
    % % x-z perturbation
    r = rand(1);
    theta = rand(1) * 2*pi;
    dx = [r*cos(theta); 0; r*sin(theta);]*eps;
    % % y perturbation
%     y = rand(1) - 0.5;
%     dx = [0; y; 0;]*eps;
    %
    x_pert(:, i) = dx;
    x0 = des_x0 + dx;
    X0 = [x0; x_dot0; int_d_x0;];
    for l = 1:N_wts
        wt = wts(l);
        [err_pos(i, l), N_conv(i, l)] =  simulate_control(gains, WK, ...
        INSECT, des, X0, N, N_single, N_period, t, wt, bound_param, err_bound);
%         fprintf('%d, %d\n', l, i);
    end
end
tocBytes(par_pool);
toc;

end

function [err_pos, N_conv, x, x_dot, int_d_x, R, Q_R, Q_L, Q_A, theta_B, theta_A, W, W_dot, W_R, W_R_dot, W_L, ...
    W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, f_tau, tau, Euler_R, ...
    Euler_R_dot, pos_err, dang] =  simulate_control(gains, WK, INSECT, des, X0, N, N_single, N_period, t, wt, bound_param, err_bound)
%%
X = zeros(N, 9);
X(1, :) = X0;
dt = t(2) - t(1);
% % Euler Method
for i=1:(N-1)
%     if mod(i-1, N_single/100) == 0
%         x0 = X(i, 1:3)';
%     end
    %
    [X_dot(:,i), R(:,:,i) Q_R(:,:,i) Q_L(:,:,i) Q_A(:,:,i) theta_B(i) theta_A(i) ...
        W(:,i) W_dot(:,i) W_R(:,i) W_R_dot(:,i) W_L(:,i) W_L_dot(:,i) W_A(:,i) ...
        W_A_dot(:,i) F_R(:,i) F_L(:,i) M_R(:,i) M_L(:,i) f_a(:,i) f_g(:,i) ...
        f_tau(:,i) tau(:,i) Euler_R(:,i) Euler_R_dot(:,i) pos_err(:, i) dang(:, i)]... 
        = eom_control(INSECT, WK, WK, t(i), X(i,:)', des, gains, i, X0(1:3), wt, bound_param);
    X(i+1, :) = X(i, :) + dt * X_dot(:, i)';
end
i = i + 1;
[X_dot(:,i), R(:,:,i) Q_R(:,:,i) Q_L(:,:,i) Q_A(:,:,i) theta_B(i) theta_A(i) ...
        W(:,i) W_dot(:,i) W_R(:,i) W_R_dot(:,i) W_L(:,i) W_L_dot(:,i) W_A(:,i) ...
        W_A_dot(:,i) F_R(:,i) F_L(:,i) M_R(:,i) M_L(:,i) f_a(:,i) f_g(:,i) ...
        f_tau(:,i) tau(:,i) Euler_R(:,i) Euler_R_dot(:,i) pos_err(:, i) dang(:, i)]... 
        = eom_control(INSECT, WK, WK, t(i), X(i,:)', des, gains, i, X0(1:3), wt, bound_param);
% [t, X]=ode45(@(t,X) eom_control(INSECT, WK, WK, t, X, des, gains, i, X0(1:3), wt), t, X0, ...
%         odeset('AbsTol',1e-6,'RelTol',1e-6));

t_single = (max(t)/N_period);
err = vecnorm(pos_err, 2, 1);
N_conv = Inf;
for j = 1:N_single:(N-N_single)
    err_pos = sum(err(j:(j+N_single))) * dt / t_single;
    if err_pos <= err_bound
        N_conv = (j-1)/N_single + 1;
        break;
    end
end
err_pos = sum(err(:, N-N_single:N)) * dt / t_single;
% err_pos = norm(pos_err(:, end));

x=X(:,1:3)';
x_dot=X(:,4:6)';
int_d_x=X(:,7:9)';

end

function [X_dot, R, Q_R, Q_L, Q_A, theta_B, theta_A, W, W_dot, W_R, W_R_dot, W_L, ...
    W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, f_tau, tau, Euler_R, ...
    Euler_R_dot, pos_err, dang]= eom_control(INSECT, WK_R, WK_L, t, X, des, gains, i, x0, wt, bound_param)
%% Dynamics along with the designed control
x=X(1:3);
x_dot=X(4:6);
int_d_x=X(7:9);
X = X(1:6);

%% Ideal values
[Euler_R, Euler_R_dot, Euler_R_ddot] = wing_kinematics(t,WK_R);
[Euler_L, Euler_L_dot, Euler_L_ddot] = wing_kinematics(t,WK_L);
[Q_R, Q_L, W_R, W_L, W_R_dot, W_L_dot] = wing_attitude(WK_R.beta, ...
    Euler_R, Euler_L, Euler_R_dot, Euler_L_dot, Euler_R_ddot, Euler_L_ddot);
[R, W] = body_attitude(t, WK_R.f, WK_R); % body
[Q_A, W_A, W_A_dot] = abdomen_attitude(t, WK_R.f, WK_R); % abdomen
[JJ_A, KK_A] = inertia_wing_sub(INSECT.m_A, INSECT.mu_A, INSECT.nu_A, INSECT.J_A, R, Q_A, X(4:6), W, W_A);
f_abd = -(JJ_A(1:3, 7:9)*W_A_dot + KK_A(1:3, 7:9)*W_A);
% f_abd(isnan(f_abd)) = 0;

[L_R, L_L, D_R, D_L]=wing_QS_aerodynamics(INSECT, ...
    W_R, W_L, W_R_dot, W_L_dot, x_dot, R, W, Q_R, Q_L);
F_R = L_R + D_R;
F_L = L_L + D_L;
f_a = R*Q_R*F_R + R*Q_L*F_L;
% f_a(isnan(f_a)) = 0;
f_total = f_a + f_abd; % Use this for gettings signs

%% Control design
% d_x = zeros(3, 1); d_x_dot = zeros(3, 1);
% for j=1:3
%     d_x(j) = des.x_fit{j}(t) - (x(j) - x0(j)); % make the position periodic
% end
d_x = des.x_fit_t(:, i) - x;
d_x_dot = des.x_dot_fit_t(:, i) - x_dot;
pos_err = INSECT.m*(gains.Kp_pos * d_x + gains.Kd_pos * d_x_dot + gains.Ki_pos * int_d_x);

% mul_phi = R * Q_R * [get_m(sign(F_R(1)), des.df_r_1_by_dphi_m); 0; get_m(sign(F_R(3)), des.df_r_3_by_dphi_m)];
% e2 = [0; 1; 0]; mul_theta = 2 * R * Q_R * hat(e2) * F_R;
mul_phi = [des.params.df_a_1_by_dphi_m(round((3 - sign(f_a(1)))/2)); 0; des.params.df_a_3_by_dphi_m];
mul_theta = [des.params.df_a_1_by_dtheta_0(round((3 - sign(f_a(1)))/2)); 0; des.params.df_a_3_by_dtheta_0];
mul_theta_A = [des.params.df_a_1_by_dtheta_A_m(round((3 - sign(f_abd(1)))/2)); 0; ...
    des.params.df_a_3_by_dtheta_A_m(round((3 - sign(f_abd(3)))/2))];
temp_A = [mul_phi(1), mul_theta(1), mul_theta_A(1);
          mul_phi(3), mul_theta(3), mul_theta_A(3);
          wt * mul_phi(1), wt * mul_theta(1), -(1-wt) * mul_theta_A(1)];
rhs = [pos_err(1); pos_err(3); 0];
dang = zeros(5, 1);
% dang(1:3) = temp_A \ (rhs);
%
if wt == 0
    dang(1:2) = temp_A(1:2,1:2) \ rhs(1:2);
else
    % Minimum norm solution
    temp_A = temp_A(1:2,:);
%     pos_err = tanh(pos_err) * 1e-2; %(norm(temp_A, 'fro') * bound_param);
%     rhs = [pos_err(1); pos_err(3); 0];
    dang(1:3) = temp_A' * ((temp_A*temp_A') \ rhs(1:2));
end
%
% dang(4) = pos_err(2) / des.params.df_a_2_by_dphi_m(round((3 - sign(pos_err(2)+des.f_a(2, i)))/2));
dang(4) = pos_err(2) / des.params.df_a_2_by_dphi_m(2);
% dang(5) = pos_err(2) / des.params.df_a_2_by_dpsi_m(1);

idx = abs(dang) > bound_param;
dang(idx) = bound_param * sign(dang(idx));
dphi_m_R = dang(1) + dang(4);
dphi_m_L = dang(1) - dang(4);
dtheta_0 = dang(2);
dtheta_A_m = dang(3);
dpsi_m = dang(5);

WK_R.phi_m = WK_R.phi_m + dphi_m_R;
WK_L.phi_m = WK_L.phi_m + dphi_m_L;
WK_R.theta_0 = WK_R.theta_0 + dtheta_0;
WK_L.theta_0 = WK_L.theta_0 + dtheta_0;
WK_R.psi_m = WK_R.psi_m + dpsi_m;
WK_L.psi_m = WK_L.psi_m - dpsi_m;
WK_R.theta_A_m = WK_R.theta_A_m + dtheta_A_m;
WK_L.theta_A_m = WK_L.theta_A_m + dtheta_A_m;

X = X(1:6);
[X_dot, R, Q_R, Q_L, Q_A, theta_B, theta_A, W, W_dot, W_R, ...
    W_R_dot, W_L, W_L_dot, W_A, W_A_dot, F_R, F_L, M_R, M_L, f_a, f_g, ...
    f_tau, tau, Euler_R, Euler_R_dot] = eom_QS_x(INSECT, WK_R, WK_L, t, X);

X_dot=[X_dot; d_x;];

end
