function sim_QS_x_hover_control
% simulate the thorax (x) trajectory along with a controller for
% given thorax attiude, wing kinematics, abdomen attitude.

evalin('base','clear all');
close all;
addpath('./modules', './sim_data', './plotting');
des=load('sim_QS_x_hover.mat',...
    'INSECT', 't', 'N', 'x', 'x_dot', 'R', 'Q_R', 'Q_L', 'W_R', 'W_L', 'f_tau',...
    'x0', 'x_dot0', 'Q_A', 'WK', 'x_ddot', 'f_a');

filename='sim_QS_x_hover_control_monte_carlo';
INSECT = des.INSECT;
WK = des.WK;
des.x_fit = cell(3, 1); des.x_dot_fit = cell(3, 1); des.f_a_fit = cell(3, 1);
% 'fourier8', 'cubicinterp'
for i=1:3
    des.x_fit{i} = fit(des.t, des.x(i, :)', 'fourier8');
    des.x_dot_fit{i} = fit(des.t, des.x_dot(i, :)', 'fourier8');
    des.f_a_fit{i} = fit(des.t, des.f_a(i, :)', 'fourier8');
end

% Values obtained from parametric study
des.df_a_1_by_dphi_m = 0.54e-3 / 0.1; % dphi_m_R > 0, dphi_m_L > 0
des.df_a_1_by_dtheta_m = -0.6e-3 / 0.1; % dtheta_m_R > 0, dtheta_m_L > 0
des.df_a_2_by_dpsi_m = 0.3e-3 / 0.1; % dpsi_m_R > 0, dpsi_m_L < 0
des.df_a_2_by_dphi_m = 3.5e-4 / 0.1; % dphi_m_R > 0, dphi_m_L < 0
des.df_a_3_by_dphi_m = 0.47e-3 / 0.1; % dphi_m_R > 0, dphi_m_L > 0
des.df_a_1_by_dtheta_A_m = 1.7e-4 / 0.1; % dtheta_A_m > 0
des.df_a_3_by_dtheta_A_m = 1.25e-4 / 0.1; % dtheta_A_m > 0
wt = 0.1; % weight, wt > 0.5 is unstable?

rng default;
eps1 = 1e-3; eps2 = 1e-1;
x0 = des.x0 + rand(3,1)*eps1;
x_dot0 = des.x_dot0 + rand(3,1)*eps2;
X0 = [x0; x_dot0;];
N = 1001; % N = 2001 causes some numerical divergence
N_period = 10;
N_single = round((N-1)/N_period);
T = N_period/WK.f;
t = linspace(0,T,N);

%% Optimization
% g0 = [gains.Kp_pos, gains.Kd_pos, gains.Ki_pos];
% lb = [-1000, -1000, -100];
% ub = [1000, 1000, 100];
% tic;
% % rng default; % For reproducibility
% options = optimoptions(@surrogateopt, 'InitialPoints', g0,...
%     'MaxFunctionEvaluations', 50000, 'MaxTime', 10000, ...
%     'PlotFcn', @surrogateoptplot);
% [gs, fval, exitflag, output] = surrogateopt(@(gs) obtain_err_gains(gs, WK, INSECT, des, X0, N, t), lb, ub, options);
% fprintf('Optimization has been completed\n');
% disp(fval);
% disp(output);
% toc;
% gs = [427.1529   15.6076  13.4983];

%% Simulation
pol = poly([-7.8 + 19i, -7.8 - 19i, -0.003]);
if ~all(real(roots(pol)) < 0)
    error('The chosen gains are not suitable');
end
gains.Kp_pos = pol(3); gains.Kd_pos = pol(2); gains.Ki_pos = pol(4);

% [t, X]=ode45(@(t,X) eom_control(INSECT, WK, WK, t, X, des, gains, i), t, X0, odeset('AbsTol',1e-6,'RelTol',1e-6));

% [err_pos x x_dot R Q_R Q_L Q_A theta_B theta_A W W_dot W_R W_R_dot W_L ...
%     W_L_dot W_A W_A_dot F_R F_L M_R M_L f_a f_g f_tau tau Euler_R ...
%     Euler_R_dot pos_err] =  simulate_control(gains, WK, INSECT, des, X0, N, N_single, N_period, t, wt);

%% Monte Carlo
N_sims = 1000;
epsilons1 = [1e-4, 1e-3, 1e-2, 1e-1]; N_eps1 = length(epsilons1);
epsilons2 = [1e-2, 1e-2, 1e-0]; N_eps2 = length(epsilons2);
wts = [0, 0.1, 0.2]; N_wts = length(wts);
err_pos = zeros(N_sims, N_eps1, N_wts);
x_pert = zeros(3, N_sims, N_eps1);
des_x0 = des.x0;
des_x_dot0 = des.x_dot0;

t1 = clock; t_max = 24*3600;
tic;
parfor i = 1:N_sims
    x_dot0 = des_x_dot0;
    for j = 1:N_eps1
        eps1 = epsilons1(j);
        dx = randn(2, 1)*eps1;
        dx(3) = dx(2); dx(2) = 0;
        x_pert(:, i, j) = dx;
        x0 = des_x0 + dx;
        X0 = [x0; x_dot0;];
        for l = 1:N_wts
            wt = wts(l);
            err_pos(i, j, l) =  simulate_control(gains, WK, ...
            INSECT, des, X0, N, N_single, N_period, t, wt);
        end
    end
%     t2 = clock;
%     if etime(t2, t1) > t_max
%         break;
%     end
end
toc;

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

function [err_pos x x_dot R Q_R Q_L Q_A theta_B theta_A W W_dot W_R W_R_dot W_L ...
    W_L_dot W_A W_A_dot F_R F_L M_R M_L f_a f_g f_tau tau Euler_R ...
    Euler_R_dot pos_err] =  simulate_control(gains, WK, INSECT, des, X0, N, N_single, N_period, t, wt)
%%
X = zeros(N, 9);
X(1, :) = [X0; zeros(3, 1)];
dt = t(2) - t(1);
% % Explicit RK4
for i=1:(N-1)
%     if mod(i-1, N_single) == 0
%         x0 = X(i, 1:3)';
%     end
    if i == 1
        f_a_im1 = zeros(6, 1);
    else
        f_a_im1(1:3, :) = f_a(1:3, i-1);
        [JJ_A, KK_A] = inertia_wing_sub(INSECT.m_A, INSECT.mu_A, INSECT.nu_A, INSECT.J_A, R(:, :, i-1), Q_A(:, :, i-1), X(i-1, 4:6)', W(:, i-1), W_A(:, i-1));
        f_a_im1(4:6, :) = -(JJ_A(1:3, 7:9)*W_A_dot(:, i-1) + KK_A(1:3, 7:9)*W_A(:, i-1));
    end
    %
    [X_dot(:,i), R(:,:,i) Q_R(:,:,i) Q_L(:,:,i) Q_A(:,:,i) theta_B(i) theta_A(i) ...
        W(:,i) W_dot(:,i) W_R(:,i) W_R_dot(:,i) W_L(:,i) W_L_dot(:,i) W_A(:,i) ...
        W_A_dot(:,i) F_R(:,i) F_L(:,i) M_R(:,i) M_L(:,i) f_a(:,i) f_g(:,i) ...
        f_tau(:,i) tau(:,i) Euler_R(:,i) Euler_R_dot(:,i) pos_err(:, i)]... 
        = eom_control(INSECT, WK, WK, t(i), X(i,:)', des, gains, i, X0(1:3), f_a_im1, wt);
    X(i+1, :) = X(i, :) + dt * X_dot(:, i)';
end
i = i + 1;
[X_dot(:,i), R(:,:,i) Q_R(:,:,i) Q_L(:,:,i) Q_A(:,:,i) theta_B(i) theta_A(i) ...
        W(:,i) W_dot(:,i) W_R(:,i) W_R_dot(:,i) W_L(:,i) W_L_dot(:,i) W_A(:,i) ...
        W_A_dot(:,i) F_R(:,i) F_L(:,i) M_R(:,i) M_L(:,i) f_a(:,i) f_g(:,i) ...
        f_tau(:,i) tau(:,i) Euler_R(:,i) Euler_R_dot(:,i) pos_err(:, i)]... 
        = eom_control(INSECT, WK, WK, t(i), X(i,:)', des, gains, i, X0(1:3), f_a_im1, wt);

err_pos = vecnorm(pos_err, 2, 1);
err_pos = trapz(t(N-N_single:N)', err_pos(N-N_single:N)) / max(t/N_period);

x=X(:,1:3)';
x_dot=X(:,4:6)';

end

function [X_dot R Q_R Q_L Q_A theta_B theta_A W W_dot W_R W_R_dot W_L ...
    W_L_dot W_A W_A_dot F_R F_L M_R M_L f_a f_g f_tau tau Euler_R ...
    Euler_R_dot pos_err]= eom_control(INSECT, WK_R, WK_L, t, X, des, gains, i, x0, f_a_im1, wt)
%% Dynamics along with the designed control

x=X(1:3);
x_dot=X(4:6);
int_d_x=X(7:9);

% Control design
d_x = zeros(3, 1); d_x_dot = zeros(3, 1);
for j=1:3
%     d_x(j) = des.x_fit{j}(t) - (x(j) - x0(j)); % Trying to make the
%     position periodic first while using this expression
    d_x(j) = des.x_fit{j}(t) - x(j);
    d_x_dot(j) = des.x_dot_fit{j}(t) - x_dot(j);
end
pos_err = INSECT.m*(gains.Kp_pos * d_x + gains.Kd_pos * d_x_dot + gains.Ki_pos * int_d_x);

% % 1st strategy
% dphi_m_R = sign(f_a_im1(3)) * pos_err(3) / des.df_a_3_by_dphi_m + sign(f_a_im1(2)) * pos_err(2) / des.df_a_2_by_dphi_m;
% dphi_m_L = sign(f_a_im1(3)) * pos_err(3) / des.df_a_3_by_dphi_m - sign(f_a_im1(2)) * pos_err(2) / des.df_a_2_by_dphi_m;
% dtheta_m = sign(f_a_im1(1)) * (pos_err(1) - ...
%     sign(f_a_im1(1)) * (dphi_m_R+dphi_m_L)/2 * des.df_a_1_by_dphi_m) / des.df_a_1_by_dtheta_m;
% % 2nd strategy
% temp_A = [sign_pos(f_a_im1(1)) * des.df_a_1_by_dphi_m, sign_pos(f_a_im1(4)) * des.df_a_1_by_dtheta_A_m ;
%           sign_pos(f_a_im1(3)) * des.df_a_3_by_dphi_m, sign_pos(f_a_im1(6)) * des.df_a_3_by_dtheta_A_m ;];
% if abs(det(temp_A)) < 1e-10
%     dang = zeros(2, 1);
% else
%     dang = temp_A \ [pos_err(1); pos_err(3)];
% end
% dphi_m = dang(1); dtheta_A_m = dang(2);
% dphi_m_R = dphi_m; dphi_m_L = dphi_m;
% dpsi_m = sign(f_a_im1(2)) * pos_err(2) / des.df_a_2_by_dpsi_m;
% dtheta_m = 0;
% % 3rd strategy
dtheta_A_m = sign(f_a_im1(6)) * wt * pos_err(3) / des.df_a_3_by_dtheta_A_m;
dphi_m_R = sign(f_a_im1(3)) * (1-wt) * pos_err(3) / des.df_a_3_by_dphi_m + ...
    sign(f_a_im1(2)) * pos_err(2) / des.df_a_2_by_dphi_m;
dphi_m_L = sign(f_a_im1(3)) * (1-wt) * pos_err(3) / des.df_a_3_by_dphi_m - ...
    sign(f_a_im1(2)) * pos_err(2) / des.df_a_2_by_dphi_m;
dtheta_m = sign(f_a_im1(1)) * (pos_err(1) - ...
    sign(f_a_im1(1)) * (dphi_m_R+dphi_m_L)/2 * des.df_a_1_by_dphi_m -...
    sign(f_a_im1(4)) * dtheta_A_m * des.df_a_1_by_dtheta_A_m) / des.df_a_1_by_dtheta_m;
dpsi_m = 0;

WK_R.phi_m = WK_R.phi_m + dphi_m_R;
WK_L.phi_m = WK_L.phi_m + dphi_m_L;
WK_R.theta_m = WK_R.theta_m + dtheta_m;
WK_L.theta_m = WK_L.theta_m + dtheta_m;
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
